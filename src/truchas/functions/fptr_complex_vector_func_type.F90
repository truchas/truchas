!!
!! FPTR_COMPLEX_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class COMPLEX_VECTOR_FUNC.
!! This implementation calls a user-provided parametrized function.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!! September 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fptr_complex_vector_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_vector_func_class
  implicit none
  private

  type, extends(complex_vector_func), public :: fptr_complex_vector_func
    private
    procedure(fptr_func), nopass, pointer :: fptr => null()
    real(r8), allocatable :: p(:) ! function parameters
  contains
    procedure :: eval
  end type

  public :: fptr_func, alloc_fptr_complex_vector_func

  abstract interface
    function fptr_func(x, p, dim) result(fx)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      integer, value :: dim
      complex(r8) :: fx(dim)
    end function
  end interface

  !! Defined constructor
  interface fptr_complex_vector_func
    procedure fptr_complex_vector_func_value
  end interface

contains

  subroutine alloc_fptr_complex_vector_func(f, dim, fptr, p)
    class(complex_vector_func), allocatable, intent(out) :: f
    integer, intent(in) :: dim
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    allocate(f, source=fptr_complex_vector_func(dim, fptr, p))
  end subroutine

  !! Constructor for FPTR_vector_FUNC objects.
  function fptr_complex_vector_func_value(dim, fptr, p) result(f)
    integer, intent(in) :: dim
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    type(fptr_complex_vector_func) :: f
    f%dim = dim
    f%fptr => fptr
    if (present(p)) then
      f%p = p
    else
      allocate(f%p(0))
    end if
  end function

  function eval(this, x) result(fx)
    class(fptr_complex_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx(this%dim)
    fx = this%fptr(x, this%p, this%dim)
  end function

end module fptr_complex_vector_func_type
