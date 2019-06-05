!!
!! CONST_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.
!! This implementation defines a constant-valued function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module const_scalar_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  implicit none
  private

  type, extends(scalar_func), public :: const_scalar_func
    !private  ! scalar_func_tools needs access
    real(r8) :: const = 0.0_r8
  contains
    procedure :: eval
  end type const_scalar_func

  !! Defined constructor
  interface const_scalar_func
    procedure const_scalar_func_value
  end interface

contains

  !! Constructor for CONST_SCALAR_FUNC objects.
  function const_scalar_func_value(const) result(f)
    real(r8), intent(in) :: const
    type(const_scalar_func) :: f
    f%const = const
  end function const_scalar_func_value

  function eval(this, x) result(fx)
    class(const_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = this%const
  end function eval

end module const_scalar_func_type
