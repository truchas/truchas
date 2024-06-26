!!
!! MTC_BNDRY_FUNC_TYPE
!!
!! This module defines an extension of the base class BNDRY_FUNC3 that
!! implements the MTC boundary condition flux function on a subset of the
!! boundary faces of a mesh type that extends the UNSTR_BASE_MESH class.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module mtc_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func3_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use bndry_face_group_builder_type
  implicit none
  private

  type, extends(bndry_func3), public :: mtc_bndry_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: f(:), g(:)
    ! temporaries used during construction
    type(bndry_face_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist, glist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
    procedure :: compute_value => compute
    procedure :: compute_deriv2 => compute
  end type mtc_bndry_func

contains

  subroutine init(this, mesh)
    class(mtc_bndry_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, mtc, tamb, setids, stat, errmsg)
    class(mtc_bndry_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: mtc, tamb
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(mtc)
    call this%glist%append(tamb)
  end subroutine add

  subroutine add_complete(this)
    class(mtc_bndry_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index)), this%deriv2(size(this%index)))
    call scalar_func_list_to_box_array(this%flist, this%f)
    call scalar_func_list_to_box_array(this%glist, this%g)
  end subroutine add_complete

  subroutine compute(this, t, var1, var2)
    class(mtc_bndry_func), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: var1(:), var2(:)
    integer :: n, j
    real(r8) :: args(-2:size(this%mesh%x,dim=1)), c
    ASSERT(allocated(this%index))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv2 => this%deriv2(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          args(-2) = var1(index(j))
          args(-1) = var2(index(j))
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:3) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
          end associate
          c = this%mesh%area(index(j)) * this%f(n)%eval(args)
          value(j) = c * (var2(index(j)) - this%g(n)%eval(args))
          deriv2(j) = c
        end do
      end associate
    end do
  end subroutine compute

end module mtc_bndry_func_type
