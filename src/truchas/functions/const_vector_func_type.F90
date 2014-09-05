!!
!! CONST_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class VECTOR_FUNC.
!! This implementation defines a constant-valued function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2014
!!

#include "f90_assert.fpp"

module const_vector_func_type

  use kinds, only: r8
  use vector_func_class
  implicit none
  private

  type, extends(vector_func), public :: const_vector_func
    real(r8), allocatable :: const(:)
  contains
    procedure :: eval
    procedure :: eval_comp
  end type const_vector_func

  interface const_vector_func
    procedure const_vector_funcx
  end interface

contains

  !! Constructor
  function const_vector_funcx (const) result (f)
    real(r8), intent(in) :: const(:)
    type(const_vector_func) :: f
    f%dim = size(const)
    f%const = const
  end function const_vector_funcx

  function eval (this, x) result (fx)
    class(const_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%dim)
    fx = this%const
  end function eval

  function eval_comp (this, i, x) result (fx)
    class(const_vector_func), intent(in) :: this
    integer,  intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = this%const(i)
  end function eval_comp

end module const_vector_func_type
