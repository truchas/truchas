!!
!! FPTR_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class VECTOR_FUNC.
!! This implementation calls a user-provided parametrized function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! October 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fptr_vector_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_func_class
  implicit none
  private

  type, extends(vector_func), public :: fptr_vector_func
    procedure(fptr_func), nopass, pointer :: fptr => null()
    real(r8), allocatable :: p(:) ! function parameters
  contains
    procedure :: eval
    procedure :: eval_comp
  end type fptr_vector_func

  public :: fptr_func

  abstract interface
    function fptr_func(x, p, dim) result(fx)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      integer, value :: dim
      real(r8) :: fx(dim)
    end function
  end interface

  !! Defined constructor
  interface fptr_vector_func
    procedure fptr_vector_func_value
  end interface

contains

  !! Constructor for FPTR_vector_FUNC objects.
  function fptr_vector_func_value(dim, fptr, p) result(f)
    integer, intent(in) :: dim
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    type(fptr_vector_func) :: f
    f%dim = dim
    f%fptr => fptr
    if (present(p)) then
      f%p = p
    else
      allocate(f%p(0))
    end if
  end function fptr_vector_func_value

  function eval(this, x) result(fx)
    class(fptr_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%dim)
    fx = this%fptr(x, this%p, this%dim)
  end function eval

  function eval_comp(this, i, x) result(fxi)
    class(fptr_vector_func), intent(in) :: this
    integer,  intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fxi, fx(this%dim)
    fx = this%fptr(x, this%p, this%dim)
    fxi = fx(i)
  end function eval_comp

end module fptr_vector_func_type
