!!
!! COMPLEX_SCALAR_FUNC_CLASS
!!
!! This module provides the abstract base class COMPLEX_SCALAR_FUNC that defines
!! an interface to a general complex scalar-valued function of a vector argument.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module complex_scalar_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: complex_scalar_func
  contains
    procedure(eval), deferred :: eval
  end type

  abstract interface
    function eval(this, x) result(fx)
      import :: complex_scalar_func, r8
      class(complex_scalar_func), intent(in) :: this
      real(r8), intent(in) :: x(:)
      complex(r8) :: fx
    end function
  end interface

end module complex_scalar_func_class
