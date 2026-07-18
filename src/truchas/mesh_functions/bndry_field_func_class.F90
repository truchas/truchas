!!
!! BNDRY_FIELD_FUNC_CLASS
!!
!! This module defines an abstract sparse boundary function for scalar boundary
!! data that depends on one mesh-wide field. Concrete implementations retain
!! the boundary entity index set and return transient values and derivatives
!! through evaluation arguments. They typically employ user-specified scalar
!! functions to compute the boundary data; the details of that evaluation are
!! left to the concrete implementation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The index array is persistent topology owned by the function object. For
!! current users, INDEX(J) identifies a boundary face managed by this function
!! object. The interpretation of the indices remains an implicit contract with
!! the physics code that holds the polymorphic object.
!!
!! Evaluation is local to each managed face: the result for entry J may depend
!! on the field value at INDEX(J), but not on field values at other entities.
!! Thus there is no field coupling between distinct entries.
!!
!! Evaluation results are returned through allocatable arrays ordered to match
!! INDEX. The concrete implementation allocates VALUE and DERIV with one entry
!! per managed face. DERIV gives the derivative with respect to the field value
!! on that face.
!!

module bndry_field_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: bndry_field_func
    integer, allocatable, public :: index(:)
  contains
    procedure(compute_value), deferred :: compute_value
    procedure(compute_deriv), deferred :: compute_deriv
  end type

  abstract interface
    subroutine compute_value(this, t, u, value)
      import r8, bndry_field_func
      class(bndry_field_func), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), allocatable, intent(out) :: value(:)
    end subroutine

    subroutine compute_deriv(this, t, u, deriv)
      import r8, bndry_field_func
      class(bndry_field_func), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), allocatable, intent(out) :: deriv(:)
    end subroutine
  end interface

end module bndry_field_func_class
