!!
!! NEAR_ZERO_FUNCTION
!!
!! This module provides a function for checking whether reals
!! are within some threshhold of zero.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module near_zero_function

  use kinds, only: r8
  implicit none
  private

  public :: near_zero

  interface near_zero
    module procedure near_zero_r8, near_zero_r8a, near_zero_r8aa
  end interface near_zero

contains

  pure logical function near_zero_r8 (x,tol)
    real(r8),           intent(in) :: x
    real(r8), optional, intent(in) :: tol

    real(r8) :: tolh

    tolh = 1e4_r8*epsilon(1.0_r8)
    if (present(tol)) tolh = tol

    near_zero_r8 = abs(x) < tolh
  end function near_zero_r8

  pure function near_zero_r8a (x,tol)
    real(r8),           intent(in) :: x(:)
    real(r8), optional, intent(in) :: tol
    logical                        :: near_zero_r8a(size(x))

    real(r8) :: tolh

    tolh = 1e4_r8*epsilon(1.0_r8)
    if (present(tol)) tolh = tol

    near_zero_r8a = abs(x) < tolh
  end function near_zero_r8a

  pure function near_zero_r8aa (x,tol)
    real(r8),           intent(in) :: x(:,:)
    real(r8), optional, intent(in) :: tol
    logical                        :: near_zero_r8aa(size(x,dim=1),size(x,dim=2))

    real(r8) :: tolh

    tolh = 1e4_r8*epsilon(1.0_r8)
    if (present(tol)) tolh = tol

    near_zero_r8aa = abs(x) < tolh
  end function near_zero_r8aa

end module near_zero_function
