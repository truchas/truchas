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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scalar_func_class

  use kinds, only: r8
  implicit none
  private

  type, abstract, public :: scalar_func
  contains
    procedure(sf_eval), deferred :: eval
  end type scalar_func

  abstract interface
    function sf_eval (this, x) result (fx)
      import :: scalar_func, r8
      class(scalar_func), intent(in) :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: fx
    end function
  end interface

end module scalar_func_class
