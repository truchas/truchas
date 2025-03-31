!!
!! BNDRY_NODE_CFUNC_TYPE
!!
!! This module defines an implementation of the base class BNDRY_CFUNC1 that
!! describes a time-dependent complex-valued function on a subset of the
!! boundary nodes of a SIMPL_MESH type mesh.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module bndry_node_cfunc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_cfunc1_class
  use simpl_mesh_type
  use complex_scalar_func_class
  use complex_scalar_func_containers
  use bndry_node_group_builder_type
  implicit none
  private

  type, extends(bndry_cfunc1), public :: bndry_node_cfunc
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: tlast = -huge(1.0_r8)
    logical :: computed = .false.
    integer :: ngroup
    integer, allocatable :: xgroup(:), faces(:), fnode(:,:)
    type(complex_scalar_func_box), allocatable :: f(:)
    integer :: compute_type
    ! temporaries used during construction
    type(bndry_node_group_builder), allocatable :: builder
    type(complex_scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type

contains

  subroutine init(this, mesh, compute_type, no_overlap)
    class(bndry_node_cfunc), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    integer, intent(in) :: compute_type
    logical, intent(in), optional :: no_overlap
    this%mesh => mesh
    this%compute_type = compute_type
    allocate(this%builder)
    call this%builder%init(mesh, no_overlap=no_overlap)
  end subroutine

  subroutine add(this, f, setids, stat, errmsg)
    class(bndry_node_cfunc), intent(inout) :: this
    class(complex_scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(f)
  end subroutine

  subroutine add_complete(this)
    class(bndry_node_cfunc), intent(inout) :: this
    integer :: stat
    ASSERT(allocated(this%builder))
    select case (this%compute_type)
    case (1)
      call this%builder%get_node_groups(this%ngroup, this%xgroup, this%index, stat)
      INSIST(stat == 0)
      deallocate(this%builder)
    case (2)
      call this%builder%get_face_groups(this%ngroup, this%xgroup, this%faces, this%fnode, this%index)
    case default
      INSIST(.false.)
    end select
    allocate(this%value(size(this%index)))
    call complex_scalar_func_list_to_box_array(this%flist, this%f)
  end subroutine

  subroutine compute(this, t)

    class(bndry_node_cfunc), intent(inout) :: this
    real(r8), intent(in) :: t

    !! Verify that THIS is in the correct state.
    ASSERT(allocated(this%index))

    if (this%computed .and. t == this%tlast) return  ! values already set for this T

    this%tlast = t
    this%computed = .true.

    select case (this%compute_type)
    case(1)

      block
        integer :: n, i
        do n = 1, this%ngroup
          do i = this%xgroup(n), this%xgroup(n+1)-1 ! loop over nodes in group n
            this%value(i) = this%f(n)%eval([this%mesh%x(:,this%index(i)), t])
          end do
        end do
      end block

    case (2)
    
      block
        use gauss_quad_tri75
        complex(r8) :: f(7), s
        integer i, j, k, l, n
        this%value = 0.0_r8
        do n = 1, this%ngroup
          do i = this%xgroup(n), this%xgroup(n+1)-1 ! loop over faces in group n
            associate (x => this%mesh%x(:,this%mesh%fnode(:,this%faces(i))))
              do k = 1, 7
                f(k) = this%f(n)%eval([matmul(x, GQTRI75_phi(:,k)), t])
              end do
            end associate
            do l = 1, size(this%fnode,dim=1)
              s = 0.0_r8
              do k = 1, 7
                s = s + GQTRI75_wgt(k)*f(k)*GQTRI75_phi(l,k)
              end do
              j = this%fnode(l,i)
              this%value(j) = this%value(j) + s*this%mesh%area(this%faces(i))
            end do
          end do
        end do
      end block

    end select

  end subroutine compute

end module bndry_node_cfunc_type
