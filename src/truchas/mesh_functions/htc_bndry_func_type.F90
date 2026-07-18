!!
!! HTC_BNDRY_FUNC_TYPE
!!
!! This module defines an extension of the base class BNDRY_FIELD_FUNC that
!! implements the HTC boundary condition flux function on a subset of the
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

module htc_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_field_func_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use bndry_face_group_builder_type
  implicit none
  private

  type, extends(bndry_field_func), public :: htc_bndry_func
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
    procedure :: compute_value
    procedure :: compute_deriv
  end type htc_bndry_func

contains

  subroutine init(this, mesh)
    class(htc_bndry_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, htc, tamb, setids, stat, errmsg)
    class(htc_bndry_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: htc, tamb
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(htc)
    call this%glist%append(tamb)
  end subroutine add

  subroutine add_complete(this)
    class(htc_bndry_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%f)
    call scalar_func_list_to_box_array(this%glist, this%g)
  end subroutine add_complete

  subroutine compute_value(this, t, u, value)
    class(htc_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, j
    real(r8) :: args(0:size(this%mesh%x,dim=1)), c
    ASSERT(allocated(this%index))
    allocate(value(size(this%index)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
          end associate
          c = this%mesh%area(index(j)) * this%f(n)%eval(args)
          value(j) = c * (u(index(j)) - this%g(n)%eval(args))
        end do
      end associate
    end do
  end subroutine compute_value

  subroutine compute_deriv(this, t, u, deriv)
    class(htc_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: deriv(:)
    integer :: n, j
    real(r8) :: args(0:size(this%mesh%x,dim=1))
    ASSERT(allocated(this%index))
    allocate(deriv(size(this%index)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => deriv(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
          end associate
          deriv(j) = this%mesh%area(index(j)) * this%f(n)%eval(args)
        end do
      end associate
    end do
  end subroutine compute_deriv

end module htc_bndry_func_type
