!!
!! VF_MATRIX_CLASS
!!
!! A common interface to distributed face-based view factor matrices used by the
!! enclosure radiosity solver.
!!
!! David Neill-Asanza <dhna@lanl.gov>, July 2019
!! Neil N. Carlson <nnc@lanl.gov>, refactoring June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!   This module defines the abstract base class VF_MATRIX which defines a
!!   generic face-base distributed view factor matrix and provides a common
!!   interface for operating on concrete instances of this class.  Extensions
!!   to this class must accept and return face-based arguments, regardless of
!!   how the enclosure patches were defined for the view factor calculation.
!!   In particular, if the VF_MATRIX includes ambient view factor data, then
!!   THIS%AMB_VF must be a face-based array.
!!
!!   Application code is expected to use polymorphic variables of this type and
!!   not work directly with extensions of the type.  The base type defines the
!!   following type bound procedures:
!!
!!   MATVEC(X, Y) computes the matrix-vector product Y = PHI*X, where PHI is
!!   the view factor matrix, and X and Y are distributed face-based vectors.
!!

module vf_matrix_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: vf_matrix
    integer :: nface, nface_tot
    real, allocatable :: amb_vf(:)
    logical :: has_ambient
    integer, allocatable :: face_map(:)
  contains
    procedure(matvec), deferred :: matvec
  end type

  abstract interface
    subroutine matvec(this, x, y)
      import :: vf_matrix, r8
      class(vf_matrix), intent(in) :: this
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: y(:)
    end subroutine
  end interface

end module vf_matrix_class
