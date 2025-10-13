!!
!! CONST_COMPLEX_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class COMPLEX_SCALAR_FUNC.
!! This implementation defines a constant-valued function.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module const_complex_scalar_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_scalar_func_class
  implicit none
  private

  public const_complex_scalar_func, alloc_const_complex_scalar_func

  type, extends(complex_scalar_func) :: const_complex_scalar_func
    private
    complex(r8) :: const = 0.0_r8
  contains
    procedure :: eval
  end type

  !! Defined constructor
  interface const_complex_scalar_func
    procedure const_complex_scalar_func_value
  end interface

contains

  subroutine alloc_const_complex_scalar_func(f, const)
    class(complex_scalar_func), allocatable, intent(out) :: f
    complex(r8), intent(in) :: const
    allocate(f, source=const_complex_scalar_func(const))
  end subroutine

  !! Constructor for CONST_SCALAR_FUNC objects.
  function const_complex_scalar_func_value(const) result(f)
    complex(r8), intent(in) :: const
    type(const_complex_scalar_func) :: f
    f%const = const
  end function 

  function eval(this, x) result(fx)
    class(const_complex_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx
    fx = this%const
  end function 

end module const_complex_scalar_func_type
