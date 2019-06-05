!!
!! SCALAR_FUNC_CLASS
!!
!! This module defines the abstract base class SCALAR_FUNC that provides
!! an interface to a general scalar-valued function of a vector argument.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scalar_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: scalar_func
  contains
    procedure(sf_eval), deferred :: eval
  end type scalar_func

  abstract interface
    function sf_eval(this, x) result(fx)
      import :: scalar_func, r8
      class(scalar_func), intent(in) :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: fx
    end function
  end interface

end module scalar_func_class
