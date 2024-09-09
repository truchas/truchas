!!
!! BNDRY_FACE_CFUNC_TYPE
!!
!! This module defines an implementation of the base class BNDRY_CFUNC1 that
!! describes a time-dependent complex-valued function on a subset of the
!! boundary faces of a SIMPL_MESH type mesh.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module bndry_face_cfunc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_cfunc1_class
  use simpl_mesh_type
  use complex_scalar_func_class
  use complex_scalar_func_containers
  use bndry_face_group_builder_type
  implicit none
  private

  type, extends(bndry_cfunc1), public :: bndry_face_cfunc
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: tlast = -huge(1.0_r8)
    logical :: computed = .false.
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(complex_scalar_func_box), allocatable :: farray(:)
    ! temporaries used during construction
    type(bndry_face_group_builder), allocatable :: builder
    type(complex_scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type

contains

  subroutine init(this, mesh, bndry_only, omit_offp)
    class(bndry_face_cfunc), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    logical, intent(in), optional :: bndry_only, omit_offp
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh, bndry_only, omit_offp)
  end subroutine

  subroutine add(this, f, setids, stat, errmsg)
    class(bndry_face_cfunc), intent(inout) :: this
    class(complex_scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(f)
  end subroutine

  subroutine add_complete(this)
    class(bndry_face_cfunc), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index)))
    call complex_scalar_func_list_to_box_array(this%flist, this%farray)
  end subroutine

  subroutine compute(this, t)

    class(bndry_face_cfunc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: n

    !! Verify that THIS is in the correct state.
    ASSERT(allocated(this%index))

    if (this%computed .and. t == this%tlast) return  ! values already set for this T

    this%tlast = t
    this%computed = .true.

    do n = 1, this%ngroup
      associate(value => this%value(this%xgroup(n):this%xgroup(n+1)-1))
        value = this%farray(n)%eval([t])
      end associate
    end do

  end subroutine compute

end module bndry_face_cfunc_type
