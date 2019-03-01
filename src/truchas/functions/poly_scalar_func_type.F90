!!
!! POLY_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.
!! This implementation defines a general polynomial with user-specified
!! coefficients and integral exponents.
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

module poly_scalar_func_type

  use kinds, only: r8
  use scalar_func_class
  implicit none
  private

  type, extends(scalar_func), public :: poly_scalar_func
    !private  ! scalar_func_tools needs access
    integer :: emin = 0 ! minimum exponent
    integer :: emax = 0 ! maximum exponent
    real(r8) :: x0 = 0.0_r8 ! reference point
    real(r8), allocatable :: c(:) ! array of coefficients
  contains
    procedure :: eval
  end type poly_scalar_func

  !! Defined constructor
  interface poly_scalar_func
    procedure poly_scalar_func_value
  end interface

contains

  !! Constructor for POLY_SCALAR_FUNC objects.
  function poly_scalar_func_value (c, e, x0) result (f)

    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:)
    real(r8), intent(in), optional :: x0
    type(poly_scalar_func) :: f

    integer :: j

    INSIST(size(c) > 0 .and. size(c) == size(e))

    f%emin = min(0, minval(e))
    f%emax = max(0, maxval(e))

    allocate(f%c(f%emin:f%emax))

    f%c = 0.0_r8
    do j = 1, size(e)
      f%c(e(j)) = f%c(e(j)) + c(j)
    end do

    if (present(x0)) f%x0 = x0

  end function poly_scalar_func_value

  function eval (this, x) result (fx)

    class(poly_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx

    integer :: j
    real(r8) :: z, w

    !! Polynomial terms with non-negative exponents.
    fx = this%c(this%emax)
    if (this%emax > 0) then
      z = x(1) - this%x0
      do j = this%emax, 1, -1
        fx = this%c(j-1) + z*fx
      end do
    end if

    !! Polynomial terms with negative exponents.
    if (this%emin < 0) then
      w = this%c(this%emin)
      z = 1.0_r8 / (x(1) - this%x0)
      do j = this%emin, -2
        w = this%c(j+1) + z*w
      end do
      fx = fx + z*w
    end if

  end function eval

end module poly_scalar_func_type
