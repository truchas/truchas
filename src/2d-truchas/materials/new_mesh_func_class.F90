!!
!! NEW_MESH_FUNC_CLASS
!!
!! This module defines the NEW_MESH_FUNC abstract type for evaluating scalar
!! mesh functions. State is indexed as (entity, state variable).
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module new_mesh_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: new_mesh_func
  contains
    procedure(compute_value), deferred :: compute_value
    procedure(compute_value_cell), deferred :: compute_value_cell
    procedure(compute_deriv), deferred :: compute_deriv
  end type

  abstract interface
    subroutine compute_value(this, state, value)
      import :: new_mesh_func, r8
      class(new_mesh_func), intent(in) :: this
      real(r8), intent(in) :: state(:,:)
      real(r8), intent(out) :: value(:)
    end subroutine

    subroutine compute_value_cell(this, cell, state, value)
      import :: new_mesh_func, r8
      class(new_mesh_func), intent(in) :: this
      integer, intent(in) :: cell
      real(r8), intent(in) :: state(:)
      real(r8), intent(out) :: value
    end subroutine

    subroutine compute_deriv(this, state, index, value)
      import :: new_mesh_func, r8
      class(new_mesh_func), intent(in) :: this
      real(r8), intent(in) :: state(:,:)
      integer, intent(in) :: index
      real(r8), intent(out) :: value(:)
    end subroutine
  end interface

end module new_mesh_func_class
