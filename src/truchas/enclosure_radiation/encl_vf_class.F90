!!
!! ENCL_VF_CLASS
!!
!! A common interface to distributed face-based enclosure view factors used
!! by the radiative heat transfer solver.
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
!! The radiation view factors for an enclosure consist of the matrix of view
!! factors between faces of the enclosure surface mesh, and the vector of
!! ambient view factors between the faces and the ambient environment in the
!! case of a partial enclosure. The abstract base class ENCL_VF defines an
!! interface to this data that is used by the radiative heat transfer solver.
!! It consists of matrix-vector product method and ambient view factor vector
!! component, together with additional metadata describing the mapping and
!! distribution of enclosure faces. The specific internal details of how the
!! view factor matrix is stored is left to concrete implementations of this
!! base class.
!!
!!  MATVEC(X, Y) computes the matrix-vector product Y = PHI*X, where PHI is
!!  the view factor matrix, and X and Y are distributed face-based vectors.
!!

module encl_vf_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: encl_vf
    integer :: nface, nface_tot
    real, allocatable :: amb_vf(:)
    logical :: has_ambient
    integer, allocatable :: face_map(:)
  contains
    procedure(matvec), deferred :: matvec
  end type

  abstract interface
    subroutine matvec(this, x, y)
      import :: encl_vf, r8
      class(encl_vf), intent(in) :: this
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: y(:)
    end subroutine
  end interface

end module encl_vf_class
