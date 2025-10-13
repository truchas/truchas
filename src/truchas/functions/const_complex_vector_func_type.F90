!!
!! CONST_COMPLEX_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class COMPLEX_VECTOR_FUNC.
!! This implementation defines a constant-valued function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module const_complex_vector_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_vector_func_class
  implicit none
  private

  public :: alloc_const_complex_vector_func

  type, extends(complex_vector_func), public :: const_complex_vector_func
    private
    complex(r8), allocatable :: const(:)
  contains
    procedure :: eval
  end type

  interface const_complex_vector_func
    procedure const_complex_vector_func_value
  end interface

contains

  subroutine alloc_const_complex_vector_func(f, const)
    class(complex_vector_func), allocatable, intent(out) :: f
    complex(r8), intent(in) :: const(:)
    allocate(f, source=const_complex_vector_func(const))
  end subroutine

  !! Constructor for CONST_COMPLEX_VECTOR_FUNC objects
  function const_complex_vector_func_value(const) result(f)
    complex(r8), intent(in) :: const(:)
    type(const_complex_vector_func) :: f
    f%dim = size(const)
    f%const = const
  end function

  function eval(this, x) result(fx)
    class(const_complex_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx(this%dim)
    fx = this%const
  end function

end module const_complex_vector_func_type
