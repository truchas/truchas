!!
!! INTFC_FIELD_FUNC_CLASS
!!
!! This module defines an abstract sparse interface function for scalar
!! interface data that depends on one mesh-wide field. Concrete implementations
!! retain the interface entity-pair index set and return transient values and
!! derivatives through evaluation arguments. They typically employ
!! user-specified scalar functions to compute the interface data; the details
!! of that evaluation are left to the concrete implementation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The index array is persistent topology owned by the function object. Its
!! columns are unique. For current users, each column INDEX(:,J) identifies one
!! pair of matching interface faces managed by this function object. The
!! interpretation of those pairs remains an implicit contract with the physics
!! code that holds the polymorphic object.
!!
!! Evaluation is local to each managed face pair: the result for pair J may
!! depend on the field values at INDEX(:,J), but not on field values at entities
!! in any other pair. Thus there is no field coupling between distinct pairs.
!!
!! Evaluation results are returned through allocatable arrays ordered to match
!! INDEX. The concrete implementation allocates VALUE with one entry per
!! managed face pair and DERIV with shape (2,size(INDEX,2)). The first dimension
!! of DERIV corresponds to the two faces in each interface pair and gives the
!! derivative with respect to the field value on that face.
!!

module intfc_field_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: intfc_field_func
    integer, allocatable, public :: index(:,:)
  contains
    procedure(compute_value), deferred :: compute_value
    procedure(compute_deriv), deferred :: compute_deriv
  end type

  abstract interface
    subroutine compute_value(this, t, u, value)
      import r8, intfc_field_func
      class(intfc_field_func), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), allocatable, intent(out) :: value(:)
    end subroutine

    subroutine compute_deriv(this, t, u, deriv)
      import r8, intfc_field_func
      class(intfc_field_func), intent(in) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), allocatable, intent(out) :: deriv(:,:)
    end subroutine
  end interface

end module intfc_field_func_class
