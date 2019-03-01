!!
!! BNDRY_VFUNC_CLASS
!!
!! An abstract base class that defines an interface used by physics kernels
!! to access transient mesh-based data associated with boundary conditions.
!! Concrete implementations of boundary conditions will extend this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The base class defines two public array components, INDEX(:) and VALUE(:,:).
!!  INDEX is a list of integer indices and VALUE is a corresponding list of
!!  real values.  While this is a fairly generic way of describing a function
!!  on a finite integer set, the expectation is that the indices refer to mesh
!!  faces or nodes on the boundary of the mesh.  Application code is expected
!!  to use polymorphic variables of this type, and not work directly with its
!!  extensions.  The type bound subroutine COMPUTE(T) is expected to fill the
!!  VALUE array with the data values corresponding to the time T.
!!

module bndry_vfunc_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: bndry_vfunc
    integer,  allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:,:)
  contains
    procedure(compute), deferred :: compute
  end type bndry_vfunc

  abstract interface
    subroutine compute(this, t)
      import r8, bndry_vfunc
      class(bndry_vfunc), intent(inout) :: this
      real(r8), intent(in) :: t
    end subroutine
  end interface

end module bndry_vfunc_class
