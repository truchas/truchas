!!
!! INTFC_FUNC3_CLASS
!!
!! An abstract base class that defines an interface used by physics kernels
!! to access transient mesh-based data associated with interface conditions
!! that depend on the value of a state variable. Concrete implementations of
!! interface conditions will extend this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The base class defines three public array components, INDEX(:,:), VALUE(:),
!!  and DERIV(:,:). INDEX is a list of integer index tuples and VALUE and DERIV
!!  are corresponding lists of real values.  While this is a fairly generic way
!!  of describing a function and its derivative on a finite integer tuple set,
!!  the expectation is that the index tuples are matching pairs of mesh faces
!!  or nodes on internal mesh interfaces.  Application code is expected to use
!!  polymorphic variables of this type, and not work directly with extensions.
!!  The type bound subroutine COMPUTE_VALUE(T, VAR) is expected to fill the
!!  VALUE array with the function values at time T and state variable value
!!  VAR. The type bound subroutine COMPUTE_DERIV(T, VAR) is expected to fill
!!  the DERIV array with the function derivative with respect to the state
!!  variable at time T and state variable value VAR. The COMPUTE(T, VAR) fills
!!  both arrays. It is acceptable for implementations of COMPUTE_VALUE and
!!  COMPUTE_DERIV to simply be renames of COMPUTE; they exist for cases where
!!  computation of the function value and its derivative are sufficiently
!!  different and costly to warrant having separate calls when only one array
!!  or the other is needed.
!!

module intfc_func3_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: intfc_func3
    integer, allocatable, public :: index(:,:)
    real(r8), allocatable, public :: value(:)
    real(r8), allocatable, public :: deriv2(:,:)
  contains
    procedure(compute), deferred :: compute
    procedure(compute), deferred :: compute_value
    procedure(compute), deferred :: compute_deriv2
  end type

  abstract interface
    subroutine compute(this, t, var1, var2)
      import r8, intfc_func3
      class(intfc_func3), intent(inout) :: this
      real(r8), intent(in) :: t, var1(:), var2(:)
    end subroutine
  end interface

end module intfc_func3_class
