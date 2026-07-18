!!
!! BNDRY_FUNC3_CLASS
!!
!! This module defines an abstract sparse boundary function interface for
!! scalar boundary data evaluated with one or two mesh-wide fields.
!! Concrete implementations retain the boundary entity index set and return
!! transient values and derivatives through evaluation arguments.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The index array is persistent topology owned by the function object. The
!! interpretation of the indices is an implicit contract with the physics code
!! that holds the polymorphic object, typically boundary faces or nodes.
!!
!! The COMPUTE_VALUE generic supports evaluation with either one or two
!! mesh-wide fields. Results are returned through caller-provided arrays ordered
!! to match INDEX. COMPUTE_DERIV returns derivatives with respect to the field
!! in a one-field evaluation, while COMPUTE_DERIV1 returns derivatives with
!! respect to U1 in a two-field evaluation.
!!

module bndry_func3_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: bndry_func3
    integer,  allocatable, public :: index(:)
  contains
    procedure(compute_value_1_iface), deferred :: compute_value_1
    procedure(compute_value_2_iface), deferred :: compute_value_2
    procedure(compute_deriv_iface),   deferred :: compute_deriv
    procedure(compute_deriv1_iface),  deferred :: compute_deriv1
    generic :: compute_value => compute_value_1, compute_value_2
  end type bndry_func3

  abstract interface
    subroutine compute_value_1_iface(this, t, u, value)
      import r8, bndry_func3
      class(bndry_func3), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), intent(out) :: value(:)
    end subroutine

    subroutine compute_value_2_iface(this, t, u1, u2, value)
      import r8, bndry_func3
      class(bndry_func3), intent(in) :: this
      real(r8), intent(in) :: t, u1(:), u2(:)
      real(r8), intent(out) :: value(:)
    end subroutine

    subroutine compute_deriv_iface(this, t, u, deriv)
      import r8, bndry_func3
      class(bndry_func3), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), intent(out) :: deriv(:)
    end subroutine

    subroutine compute_deriv1_iface(this, t, u1, u2, deriv1)
      import r8, bndry_func3
      class(bndry_func3), intent(in) :: this
      real(r8), intent(in) :: t, u1(:), u2(:)
      real(r8), intent(out) :: deriv1(:)
    end subroutine
  end interface

end module bndry_func3_class
