!!
!! HTC_BNDRY_FUNC_TYPE
!!
!! This module defines a heat-transfer-coefficient boundary flux function on a
!! subset of mesh boundary faces. The function stores the sparse face set and
!! evaluates user-supplied HTC and ambient-temperature scalar functions for
!! thermal transport models.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The returned flux value is area*h*(T-Tamb). The HTC function h receives
!! anonymous arguments (T,time,x,y,z), while the ambient-temperature function
!! receives (time,x,y,z).
!!
!! COMPUTE_DERIV returns the derivative with respect to T. It includes a
!! centered finite-difference approximation of the temperature dependence of h
!! but treats Tamb as independent of T.

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
    real(r8) :: args(-1:size(this%mesh%x,dim=1))
    real(r8) :: args0(0:size(this%mesh%x,dim=1)), c
    ASSERT(allocated(this%index))
    allocate(value(size(this%index)))
    args(0) = t
    args0(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          args(-1) = u(index(j))
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
            args0(1:) = args(1:)
          end associate
          c = this%mesh%area(index(j)) * this%f(n)%eval(args)
          value(j) = c * (u(index(j)) - this%g(n)%eval(args0))
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
    real(r8) :: args(-1:size(this%mesh%x,dim=1))
    real(r8) :: args0(0:size(this%mesh%x,dim=1))
    real(r8) :: area, dh, h, tamb
    ASSERT(allocated(this%index))
    allocate(deriv(size(this%index)))
    args(0) = t
    args0(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => deriv(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          args(-1) = u(index(j))
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
            args0(1:) = args(1:)
          end associate
          area = this%mesh%area(index(j))
          h = this%f(n)%eval(args)
          tamb = this%g(n)%eval(args0)
          call compute_first_arg_deriv(this%f(n), args, dh)
          deriv(j) = area * (h + dh * (u(index(j)) - tamb))
        end do
      end associate
    end do
  end subroutine compute_deriv

  subroutine compute_first_arg_deriv(f, args, deriv)
    type(scalar_func_box), intent(in) :: f
    real(r8), intent(inout) :: args(:)
    real(r8), intent(out) :: deriv
    real(r8) :: fdinc, fm, fp, u

    ! Differentiate with respect to the first function argument.
    u = args(1)
    fdinc = max(1.0_r8, abs(u)) * sqrt(epsilon(1.0_r8))
    fdinc = scale(1.0_r8, exponent(fdinc))
    args(1) = u + fdinc
    fp = f%eval(args)
    args(1) = u - fdinc
    fm = f%eval(args)
    args(1) = u
    deriv = (fp - fm) / (2*fdinc)
  end subroutine compute_first_arg_deriv

end module htc_bndry_func_type
