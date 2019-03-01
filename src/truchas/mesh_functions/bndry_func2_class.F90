!!
!! BNDRY_FUNC2_CLASS
!!
!! An abstract base class that defines an interface used by physics kernels
!! to access transient mesh-based data associated with boundary conditions
!! that depend on the value of a state variable. Concrete implementations of
!! boundary conditions will extend this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bndry_func2_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: bndry_func2
    integer,  allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:)
    real(r8), allocatable, public :: deriv(:)
  contains
    procedure(compute), deferred :: compute
    procedure(compute), deferred :: compute_value
    procedure(compute), deferred :: compute_deriv
  end type bndry_func2

  abstract interface
    subroutine compute(this, t, var)
      import r8, bndry_func2
      class(bndry_func2), intent(inout) :: this
      real(r8), intent(in) :: t, var(:)
    end subroutine
  end interface

end module bndry_func2_class
