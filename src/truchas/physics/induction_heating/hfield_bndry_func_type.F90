#include "f90_assert.fpp"

module hfield_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func1_class
  use simpl_mesh_type
  use scalar_func_containers
  use vector_func_containers
  use bndry_edge_group_builder_type
  implicit none
  private

  type, extends(bndry_func1), public :: hfield_bndry_func
    type(simpl_mesh), pointer :: mesh => null()
    integer :: ngroup = 0
    integer, allocatable :: xgroup(:)
    real(r8), allocatable :: gvalue(:)
    type(scalar_func_box), allocatable :: f(:)
    logical :: computed = .false.
    real(r8) :: tlast = -huge(1.0_r8)
    ! temporaries used during construction
    type(bndry_edge_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
    type(vector_func_list) :: glist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type

contains

  subroutine init(this, mesh)
    class(hfield_bndry_func), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh, no_overlap=.false.)
  end subroutine

  subroutine add(this, f, g, setids, stat, errmsg)
    class(hfield_bndry_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    class(vector_func), allocatable, intent(inout) :: g
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg) ! TODO: onP only?
    if (stat /= 0) return
    call this%flist%append(f)
    call this%glist%append(g)
  end subroutine

  subroutine add_complete(this, omit_edge_list)

    class(hfield_bndry_func), intent(inout) :: this
    integer, intent(in), optional :: omit_edge_list(:)

    integer :: j, n, stat
    integer, allocatable :: xgroup(:), index(:), emap(:)
    type(vector_func_box), allocatable :: g(:)
    real(r8) :: c, p(3)

    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, xgroup, index)
    call this%builder%get_edge_groups(this%ngroup, this%xgroup, this%index, stat, omit_edge_list)
    INSIST(stat == 0)
    deallocate(this%builder)

    !! Flatten function lists to arrays
    call scalar_func_list_to_box_array(this%flist, this%f)
    call vector_func_list_to_box_array(this%glist, g)

    allocate(this%value(size(this%index)), this%gvalue(size(this%index)), emap(this%mesh%nedge))

    this%gvalue = 0.0_r8
    do n = 1, this%ngroup
      !! EMAP maps on-process edge indices to their group index, or 0 if not in a group.
      emap = 0  !TODO: replace by more efficient algorithm/structure?
      do j = this%xgroup(n), this%xgroup(n+1)-1
        emap(this%index(j)) = j
      end do
      associate (face => index(xgroup(n):xgroup(n+1)-1))
        do j = 1, size(face) ! loop over faces in group n
          associate (edge => emap(this%mesh%fedge(:,face(j))))
            !! Compute integrals of G over face and assemble results to the
            !! edges. Orientation of the face relative to the outward boundary
            !! normal must be accounted for and is inferred from the fcell array.
            c = merge(1.0_r8, -1.0_r8, this%mesh%fcell(2,face(j)) == 0) / 6.0_r8
            call project_on_face_edges(this%mesh%x(:,this%mesh%fnode(:,face(j))), g(n), p)
            if (edge(1) > 0) this%gvalue(edge(1)) = this%gvalue(edge(1)) + c * (p(3) - p(2))
            if (edge(2) > 0) this%gvalue(edge(2)) = this%gvalue(edge(2)) - c * (p(1) - p(3))
            if (edge(3) > 0) this%gvalue(edge(3)) = this%gvalue(edge(3)) + c * (p(2) - p(1))
          end associate
        end do
      end associate
    end do

  end subroutine add_complete

  !! Length-weighted edge projections using trapezoid rule
  subroutine project_on_face_edges(x, g, p)
    real(r8), intent(in) :: x(:,:)
    type(vector_func_box), intent(in) :: g
    real(r8), intent(out):: p(:)
    integer :: k
    real(r8) :: gv(3,3)
    do k = 1, 3
      gv(:,k) = g%eval(x(:,k))
    end do
    p(1) = 0.5_r8 * dot_product(x(:,3)-x(:,2), gv(:,3)+gv(:,2))
    p(2) = 0.5_r8 * dot_product(x(:,1)-x(:,3), gv(:,1)+gv(:,3))
    p(3) = 0.5_r8 * dot_product(x(:,2)-x(:,1), gv(:,2)+gv(:,1))
  end subroutine

  subroutine compute(this, t)
    class(hfield_bndry_func), intent(inout) :: this
    real(r8), intent(in) :: t
    integer :: j, n
    real(r8) :: s
    if (this%computed .and. t == this%tlast) return ! value already set for this T
    do n = 1, this%ngroup
      s = this%f(n)%eval([t])
      do j = this%xgroup(n), this%xgroup(n+1)-1
        this%value(j) = s*this%gvalue(j)
      end do
    end do
    this%tlast = t
    this%computed = .true.
  end subroutine

end module hfield_bndry_func_type
