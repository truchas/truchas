#include "f90_assert.fpp"

module fd_robin_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_cfunc1_class
  use simpl_mesh_type
  use complex_vector_func_containers
  use bndry_face_group_builder_type
  implicit none
  private

  type, extends(bndry_cfunc1), public :: fd_robin_bndry_func
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    integer :: ngroup = 0
    integer, allocatable :: xgroup(:), face_index(:), fedge(:,:)
    type(complex_vector_func_box), allocatable :: g(:)
    logical :: computed = .false.
    real(r8) :: tlast = -huge(1.0_r8)
    ! temporaries used during construction
    type(bndry_face_group_builder), allocatable :: builder
    type(complex_vector_func_list) :: glist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type

contains

  subroutine init(this, mesh)
    class(fd_robin_bndry_func), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh, omit_offp=.false., no_overlap=.true.)
  end subroutine

  subroutine add(this, g, setids, stat, errmsg)
    class(fd_robin_bndry_func), intent(inout) :: this
    class(complex_vector_func), allocatable, intent(inout) :: g
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%glist%append(g)
  end subroutine

  subroutine add_complete(this)

    use integer_map_type

    class(fd_robin_bndry_func), intent(inout) :: this

    integer :: i, j, n
    type(integer_map) :: edge_map

    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%face_index)
    deallocate(this%builder)

    !! Flatten function list to array
    call complex_vector_func_list_to_box_array(this%glist, this%g)

    !! Tag edges contained in the specified boundary faces.
    do j = 1, size(this%face_index)
      do i = 1, size(this%mesh%fedge,dim=1)
        call edge_map%set(this%mesh%fedge(i,this%face_index(j)), 0)
      end do
    end do

    n = edge_map%size()
    allocate(this%index(n), this%value(n))

    !! Generate the list of tagged edges and the inverse mapping of
    !! edges to tagged edges.
    n = 0
    do j = 1, this%mesh%nedge
      if (edge_map%contains(j)) then
        n = n + 1
        this%index(n) = j
        call edge_map%set(j, n)
      end if
    end do

    !! Extract the section of the mesh%fedge array for the specified boundary
    !! faces and re-index the values to point to the tagged edges.
    allocate(this%fedge(size(this%mesh%fedge,dim=1),size(this%face_index)))
    do j = 1, size(this%face_index)
      do i = 1, size(this%mesh%fedge,dim=1)
        this%fedge(i,j) = edge_map%val(this%mesh%fedge(i,this%face_index(j)))
      end do
    end do

  end subroutine add_complete


  subroutine compute(this, t)

    use simplex_geometry, only: tri_edge_normal
    use gauss_quad_tri75

    class(fd_robin_bndry_func), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, j, k, n
    real(r8) :: nxq(3,3)
    complex(r8) :: g(3), g_dot_nxq(7,3), s

    if (this%computed .and. t == this%tlast) return ! value already set for this T

    this%tlast = t
    this%computed = .true.

    this%value = 0.0_r8

    do n = 1, this%ngroup
      do i = this%xgroup(n), this%xgroup(n+1)-1 ! loop over faces in group n
        associate (x => this%mesh%x(:,this%mesh%fnode(:,this%face_index(i))))

          nxq = -tri_edge_normal(x)

          do k = 1, 7
            g = this%g(n)%eval([matmul(x, GQTRI75_phi(:,k)), t])
            g_dot_nxq(k,:) = matmul(g, nxq)
          end do

          ! contribution to face edge 1 (23)
          s = 0.0_r8
          do k = 1, 7
            s = s + GQTRI75_wgt(k)*&
                (GQTRI75_phi(2,k)*g_dot_nxq(k,3) - GQTRI75_phi(3,k)*g_dot_nxq(k,2))
          end do
          j = this%fedge(1,i)
          this%value(j) = this%value(j) + 0.5_r8 * s

          ! contribution to face edge 2 (13)
          s = 0.0_r8
          do k = 1, 7
            s = s + GQTRI75_wgt(k)*&
                (GQTRI75_phi(1,k)*g_dot_nxq(k,3) - GQTRI75_phi(3,k)*g_dot_nxq(k,1))
          end do
          j = this%fedge(2,i)
          this%value(j) = this%value(j) + 0.5_r8 * s

          ! contribution to face edge 3 (12)
          s = 0.0_r8
          do k = 1, 7
            s = s + GQTRI75_wgt(k)*&
                (GQTRI75_phi(1,k)*g_dot_nxq(k,2) - GQTRI75_phi(2,k)*g_dot_nxq(k,1))
          end do
          j = this%fedge(3,i)
          this%value(j) = this%value(j) + 0.5_r8 * s
        end associate
      end do
    end do

    !NB: Only the values for on-process edges are complete; off-process edge
    !    values should be ignored or overwritten, directly or indirectly.

  end subroutine compute

end module fd_robin_bndry_func_type
