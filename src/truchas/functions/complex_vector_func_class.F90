!!
!! COMPLEX_VECTOR_FUNC_CLASS
!!
!! This module provides the abstract base class COMPLEX_VECTOR_FUNC that defines
!! an interface to a general complex vector-valued function of a vector argument.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module complex_vector_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: complex_vector_func
    integer :: dim  ! the size of the vector value
  contains
    procedure(eval), deferred :: eval
  end type

  abstract interface
    function eval(this, x) result(fx)
      import :: complex_vector_func, r8
      class(complex_vector_func), intent(in) :: this
      real(r8), intent(in) :: x(:)
      complex(r8) :: fx(this%dim)
    end function
  end interface

end module complex_vector_func_class
