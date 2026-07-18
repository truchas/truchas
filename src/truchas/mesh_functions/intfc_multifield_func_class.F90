!!
!! INTFC_MULTIFIELD_FUNC_CLASS
!!
!! This module defines an abstract sparse interface function for scalar
!! interface data evaluated with one or two mesh-wide fields.
!! Concrete implementations retain the interface entity-pair index set and
!! return transient values and derivatives through evaluation arguments. They
!! typically employ user-specified scalar functions to compute the interface
!! data; the details of that evaluation are left to the concrete implementation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The index array is persistent topology owned by the function object. For
!! current users, each column INDEX(:,J) identifies one pair of matching
!! interface faces managed by this function object. The interpretation of those
!! pairs remains an implicit contract with the physics code that holds the
!! polymorphic object.
!!
!! Evaluation is local to each managed face pair: the result for pair J may
!! depend on the supplied field values at INDEX(:,J), but not on field values at
!! entities in any other pair. Thus there is no field coupling between distinct
!! pairs.
!!
!! The COMPUTE_VALUE generic supports evaluation with either one or two
!! mesh-wide fields. Results are returned through allocatable arrays ordered to
!! match INDEX. The concrete implementation allocates value arrays with one
!! entry per managed face pair and derivative arrays with shape
!! (2,size(INDEX,2)). COMPUTE_DERIV returns derivatives with respect to the
!! field in a one-field evaluation, while COMPUTE_DERIV1 returns derivatives
!! with respect to U1 in a two-field evaluation. The first dimension of a
!! derivative array corresponds to the two entities in each interface pair.
!!

module intfc_multifield_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: intfc_multifield_func
    integer, allocatable, public :: index(:,:)
  contains
    procedure(compute_value_1), deferred :: compute_value_1
    procedure(compute_value_2), deferred :: compute_value_2
    procedure(compute_deriv),   deferred :: compute_deriv
    procedure(compute_deriv1),  deferred :: compute_deriv1
    generic :: compute_value => compute_value_1, compute_value_2
  end type

  abstract interface
    subroutine compute_value_1(this, t, u, value)
      import r8, intfc_multifield_func
      class(intfc_multifield_func), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), allocatable, intent(out) :: value(:)
    end subroutine

    subroutine compute_value_2(this, t, u1, u2, value)
      import r8, intfc_multifield_func
      class(intfc_multifield_func), intent(in) :: this
      real(r8), intent(in) :: t, u1(:), u2(:)
      real(r8), allocatable, intent(out) :: value(:)
    end subroutine

    subroutine compute_deriv(this, t, u, deriv)
      import r8, intfc_multifield_func
      class(intfc_multifield_func), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), allocatable, intent(out) :: deriv(:,:)
    end subroutine

    subroutine compute_deriv1(this, t, u1, u2, deriv1)
      import r8, intfc_multifield_func
      class(intfc_multifield_func), intent(in) :: this
      real(r8), intent(in) :: t, u1(:), u2(:)
      real(r8), allocatable, intent(out) :: deriv1(:,:)
    end subroutine
  end interface

end module intfc_multifield_func_class
