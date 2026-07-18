!!
!! MTC_BNDRY_FUNC_TYPE
!!
!! This module defines a mass-transfer-coefficient boundary flux function on a
!! subset of mesh boundary faces. The function stores the sparse face set and
!! evaluates user-supplied MTC and ambient concentration scalar functions for
!! one-field or two-field species transport models.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The returned flux value is area*h*(C-Camb). The MTC function h receives
!! anonymous arguments (C,time,x,y,z) for a one-field evaluation and
!! (C,T,time,x,y,z) for a two-field evaluation. The ambient concentration
!! function receives (time,x,y,z) or (T,time,x,y,z); it does not receive C.
!!
!! COMPUTE_DERIV returns the derivative with respect to C in a one-field
!! evaluation, and COMPUTE_DERIV1 returns it in a two-field evaluation. These
!! derivatives include a centered finite-difference approximation of the
!! concentration dependence of h, but treat Camb as independent of C.
!!

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
    procedure :: compute_value_1
    procedure :: compute_value_2
    procedure :: compute_deriv
    procedure :: compute_deriv1
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
    call scalar_func_list_to_box_array(this%flist, this%f)
    call scalar_func_list_to_box_array(this%glist, this%g)
  end subroutine add_complete

  subroutine compute_value_1(this, t, u, value)
    class(mtc_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1)), args0(0:size(this%mesh%x,dim=1)), c
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
  end subroutine compute_value_1

  subroutine compute_value_2(this, t, u1, u2, value)
    class(mtc_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u1(:), u2(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, j
    real(r8) :: args(-2:size(this%mesh%x,dim=1)), args2(-1:size(this%mesh%x,dim=1)), c
    ASSERT(allocated(this%index))
    allocate(value(size(this%index)))
    args(0) = t
    args2(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          args(-2) = u1(index(j))
          args(-1) = u2(index(j))
          args2(-1) = u2(index(j))
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
            args2(1:) = args(1:)
          end associate
          c = this%mesh%area(index(j)) * this%f(n)%eval(args)
          value(j) = c * (u1(index(j)) - this%g(n)%eval(args2))
        end do
      end associate
    end do
  end subroutine compute_value_2

  subroutine compute_deriv(this, t, u, deriv)
    class(mtc_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: deriv(:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1)), args0(0:size(this%mesh%x,dim=1))
    real(r8) :: c, camb, df, fdinc, fm, fp, u0
    ASSERT(allocated(this%index))
    allocate(deriv(size(this%index)))
    args(0) = t
    args0(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => deriv(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          u0 = u(index(j))
          args(-1) = u0
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
            args0(1:) = args(1:)
          end associate
          c = this%mesh%area(index(j)) * this%f(n)%eval(args)
          camb = this%g(n)%eval(args0)
          fdinc = max(1.0_r8, abs(u0)) * sqrt(epsilon(1.0_r8))
          fdinc = scale(1.0_r8, exponent(fdinc))
          args(-1) = u0 + fdinc
          fp = this%f(n)%eval(args)
          args(-1) = u0 - fdinc
          fm = this%f(n)%eval(args)
          df = (fp - fm) / (2*fdinc)
          deriv(j) = c + this%mesh%area(index(j)) * df * (u0 - camb)
        end do
      end associate
    end do
  end subroutine compute_deriv

  subroutine compute_deriv1(this, t, u1, u2, deriv1)
    class(mtc_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u1(:), u2(:)
    real(r8), allocatable, intent(out) :: deriv1(:)
    integer :: n, j
    real(r8) :: args(-2:size(this%mesh%x,dim=1)), args2(-1:size(this%mesh%x,dim=1))
    real(r8) :: c, camb, df, fdinc, fm, fp, u10
    ASSERT(allocated(this%index))
    allocate(deriv1(size(this%index)))
    args(0) = t
    args2(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv1 => deriv1(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          u10 = u1(index(j))
          args(-2) = u10
          args(-1) = u2(index(j))
          args2(-1) = u2(index(j))
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
            args2(1:) = args(1:)
          end associate
          c = this%mesh%area(index(j)) * this%f(n)%eval(args)
          camb = this%g(n)%eval(args2)
          fdinc = max(1.0_r8, abs(u10)) * sqrt(epsilon(1.0_r8))
          fdinc = scale(1.0_r8, exponent(fdinc))
          args(-2) = u10 + fdinc
          fp = this%f(n)%eval(args)
          args(-2) = u10 - fdinc
          fm = this%f(n)%eval(args)
          df = (fp - fm) / (2*fdinc)
          deriv1(j) = c + this%mesh%area(index(j)) * df * (u10 - camb)
        end do
      end associate
    end do
  end subroutine compute_deriv1

end module mtc_bndry_func_type
