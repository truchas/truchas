!!
!! BNDRY_CFUNC1_CLASS
!!
!! An abstract base class that defines an interface used by physics kernels
!! to access transient mesh-based data associated with boundary conditions.
!! Concrete implementations of boundary conditions will extend this type.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bndry_cfunc1_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: bndry_cfunc1
    integer, allocatable :: index(:)
    complex(r8), allocatable :: value(:)
  contains
    procedure(compute), deferred :: compute
  end type

  abstract interface
    subroutine compute(this, t)
      import r8, bndry_cfunc1
      class(bndry_cfunc1), intent(inout) :: this
      real(r8), intent(in) :: t
    end subroutine
  end interface

end module bndry_cfunc1_class
