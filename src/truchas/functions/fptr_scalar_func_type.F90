!!
!! FPTR_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.
!! This implementation calls a user-provided parametrized function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!

#include "f90_assert.fpp"

module fptr_scalar_func_type

  use kinds, only: r8
  use scalar_func_class
  implicit none
  private

  type, extends(scalar_func), public :: fptr_scalar_func
    procedure(fptr_func), nopass, pointer :: fptr => null()
    real(r8), allocatable :: p(:) ! function parameters
  contains
    procedure :: eval
  end type fptr_scalar_func

  public :: fptr_func

  abstract interface
    function fptr_func (x, p) result (fx)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      real(r8) :: fx
    end function
  end interface

  !! Defined constructor
  interface fptr_scalar_func
    procedure fptr_scalar_func_value
  end interface

contains

  !! Constructor for FPTR_SCALAR_FUNC objects.
  function fptr_scalar_func_value (fptr, p) result (f)
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    type(fptr_scalar_func) :: f
    f%fptr => fptr
    if (present(p)) then
      f%p = p
    else
      allocate(f%p(0))
    end if
  end function fptr_scalar_func_value

  function eval (this, x) result (fx)
    class(fptr_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = this%fptr(x, this%p)
  end function eval

end module fptr_scalar_func_type
