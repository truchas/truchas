!!
!! VECTOR_FUNC_CLASS
!!
!! This module defines the abstract base class VECTOR_FUNC that provides
!! an interface to a general vector-valued function of a vector argument.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vector_func_class

  use kinds, only: r8
  implicit none
  private

  type, abstract, public :: vector_func
    !private
    integer :: dim  ! the size of the vector value
  contains
    procedure(vf_eval), deferred :: eval
    procedure(vf_eval_comp), deferred :: eval_comp
  end type vector_func

  abstract interface
    function vf_eval (this, x) result (fx)
      import :: vector_func, r8
      class(vector_func), intent(in) :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: fx(this%dim)
    end function
  end interface

  abstract interface
    function vf_eval_comp (this, i, x) result (fx)
      import :: vector_func, r8
      class(vector_func), intent(in) :: this
      integer, intent(in) :: i
      real(r8), intent(in) :: x(:)
      real(r8) :: fx
    end function
  end interface

end module vector_func_class
