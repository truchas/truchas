!!
!! FPTR_COMPLEX_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class COMPLEX_SCALAR_FUNC.
!! This implementation calls a user-provided parametrized function.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!! December 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fptr_complex_scalar_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_scalar_func_class
  implicit none
  private

  type, extends(complex_scalar_func), public :: fptr_complex_scalar_func
    procedure(fptr_func), nopass, pointer :: fptr => null()
    real(r8), allocatable :: p(:) ! function parameters
  contains
    procedure :: eval
  end type

  !public :: fptr_func, alloc_fptr_complex_scalar_func
  public :: alloc_fptr_complex_scalar_func

  abstract interface
    function fptr_func(x, p) result(fx)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      complex(r8) :: fx
    end function
  end interface

  !! Defined constructor
  interface fptr_complex_scalar_func
    procedure fptr_complex_scalar_func_value
  end interface

contains

  subroutine alloc_fptr_complex_scalar_func(f, fptr, p)
    class(complex_scalar_func), allocatable, intent(out) :: f
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    allocate(f, source=fptr_complex_scalar_func(fptr, p))
  end subroutine

  !! Constructor for FPTR_COMPLEX_SCALAR_FUNC objects.
  function fptr_complex_scalar_func_value(fptr, p) result(f)
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    type(fptr_complex_scalar_func) :: f
    f%fptr => fptr
    if (present(p)) then
      f%p = p
    else
      allocate(f%p(0))
    end if
  end function

  function eval(this, x) result(fx)
    class(fptr_complex_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx
    fx = this%fptr(x, this%p)
  end function

end module fptr_complex_scalar_func_type
