!!
!! BNDRY_FUNC_CLASS
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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The base class defines two public array components, INDEX(:) and VALUE(:).
!!  INDEX is a list of integer indices and VALUE is a corresponding list of
!!  real values.  While this is a fairly generic way of describing a function
!!  on a finite integer set, the expectation is that the indices refer to mesh
!!  faces or nodes on the boundary of the mesh.  Application code is expected
!!  to use polymorphic variables of this type, and not work directly with its
!!  extensions.  The type bound subroutine COMPUTE(T) is expected to fill the
!!  VALUE array with the data values corresponding to the time T.
!!

module tdep_bndry_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: tdep_bndry_func
    integer,  allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:)
    real(r8), allocatable, public :: deriv(:)
  contains
    procedure(compute), deferred :: compute
    procedure(compute), deferred :: compute_value
    procedure(compute), deferred :: compute_deriv
   end type tdep_bndry_func

  abstract interface
    subroutine compute(this, t, var)
      import r8, tdep_bndry_func
      class(tdep_bndry_func), intent(inout) :: this
      real(r8), intent(in) :: t,var(:)
    end subroutine
  end interface

end module tdep_bndry_func_class
