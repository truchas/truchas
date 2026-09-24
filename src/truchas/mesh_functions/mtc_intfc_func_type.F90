!!
!! MTC_INTFC_FUNC_TYPE
!!
!! This module defines a mass-transfer-coefficient interface flux function on
!! a subset of mesh interface face pairs. The function stores the sparse
!! interface-link set and evaluates user-supplied MTC scalar functions for
!! one-field or two-field species transport models.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Each value is the flux leaving INDEX(1,J) and entering INDEX(2,J). The MTC
!! coefficient is evaluated using the maximum concentration across the two
!! faces. The one-field user-function arguments are (max(C(face1),C(face2)),
!! time,x,y,z). The two-field arguments are (max(C(face1),C(face2)),
!! max(T(face1),T(face2)),time,x,y,z).
!!
!! COMPUTE_DERIV returns derivatives with respect to C in a one-field
!! evaluation, and COMPUTE_DERIV1 returns them in a two-field evaluation. These
!! derivatives include a centered finite-difference approximation of the MTC
!! coefficient's dependence on C.
!!

#include "f90_assert.fpp"

module mtc_intfc_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use intfc_multifield_func_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use intfc_link_group_builder_type
  implicit none
  private

  type, extends(intfc_multifield_func), public :: mtc_intfc_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: f(:)
    ! temporaries used during construction
    type(intfc_link_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute_value_1
    procedure :: compute_value_2
    procedure :: compute_deriv
    procedure :: compute_deriv1
  end type mtc_intfc_func

contains

  subroutine init(this, mesh)
    class(mtc_intfc_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, mtc, setids, stat, errmsg)
    class(mtc_intfc_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: mtc
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_link_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(mtc)
  end subroutine add

  subroutine add_complete(this)
    class(mtc_intfc_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_link_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%f)
  end subroutine add_complete

  subroutine compute_value_1(this, t, u, value)
    class(mtc_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1)), c1
    ASSERT(allocated(this%index))
    allocate(value(size(this%index,2)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                value => value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => u(index(1,j)), v2 => u(index(2,j)))
            args(-1) = max(v1, v2)
            c1 = this%f(n)%eval(args) * this%mesh%area(index(1,j))
            value(j) = c1*(v1 - v2)
          end associate
        end do
      end associate
    end do
  end subroutine compute_value_1

  subroutine compute_value_2(this, t, u1, u2, value)
    class(mtc_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, u1(:), u2(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, j
    real(r8) :: args(-2:size(this%mesh%x,dim=1)), c1
    ASSERT(allocated(this%index))
    allocate(value(size(this%index,2)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                value => value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => u1(index(1,j)), v2 => u1(index(2,j)))
            args(-2) = max(v1, v2)
            args(-1) = max(u2(index(1,j)), u2(index(2,j)))
            c1 = this%f(n)%eval(args) * this%mesh%area(index(1,j))
            value(j) = c1*(v1 - v2)
          end associate
        end do
      end associate
    end do
  end subroutine compute_value_2

  subroutine compute_deriv(this, t, u, deriv)
    class(mtc_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), allocatable, intent(out) :: deriv(:,:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1))
    real(r8) :: c1, c2, df, umax
    ASSERT(allocated(this%index))
    allocate(deriv(2,size(this%index,2)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => deriv(:,this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => u(index(1,j)), v2 => u(index(2,j)))
            umax = max(v1, v2)
            args(-1) = umax
            c1 = this%f(n)%eval(args) * this%mesh%area(index(1,j))
            call compute_first_arg_deriv(this%f(n), args, df)
            c2 = df * (v1 - v2) * this%mesh%area(index(1,j))
            deriv(1,j) =  c1 + merge(c2, 0.0_r8, v1 > v2)
            deriv(2,j) = -c1 + merge(c2, 0.0_r8, v2 > v1)
          end associate
        end do
      end associate
    end do
  end subroutine compute_deriv

  subroutine compute_deriv1(this, t, u1, u2, deriv1)
    class(mtc_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, u1(:), u2(:)
    real(r8), allocatable, intent(out) :: deriv1(:,:)
    integer :: n, j
    real(r8) :: args(-2:size(this%mesh%x,dim=1))
    real(r8) :: c1, c2, df, umax
    ASSERT(allocated(this%index))
    allocate(deriv1(2,size(this%index,2)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                deriv1 => deriv1(:,this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => u1(index(1,j)), v2 => u1(index(2,j)))
            umax = max(v1, v2)
            args(-2) = umax
            args(-1) = max(u2(index(1,j)), u2(index(2,j)))
            c1 = this%f(n)%eval(args) * this%mesh%area(index(1,j))
            call compute_first_arg_deriv(this%f(n), args, df)
            c2 = df * (v1 - v2) * this%mesh%area(index(1,j))
            deriv1(1,j) =  c1 + merge(c2, 0.0_r8, v1 > v2)
            deriv1(2,j) = -c1 + merge(c2, 0.0_r8, v2 > v1)
          end associate
        end do
      end associate
    end do
  end subroutine compute_deriv1

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

end module mtc_intfc_func_type
