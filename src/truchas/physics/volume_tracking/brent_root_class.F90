!!
!! BRENT_ROOT_CLASS
!!
!! This module provides an interface to a Brent's root-finding algorithm.
!! Note this could be replaced with some canned package,
!! like the GNU Scientific Library, netlib, or something else.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! July 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module brent_root_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: brent_root
    real(r8) :: eps = 0, feps = 0
    integer  :: maxitr, numitr = 0
  contains
    procedure, non_overridable :: find_root
    procedure(func), deferred :: f
  end type brent_root

  abstract interface
    function func (this, x) result (fx)
      import brent_root, r8
      class(brent_root), intent(inout) :: this
      real(r8), intent(in) :: x
      real(r8) :: fx
    end function func
  end interface

contains

  subroutine find_root (this, xmin, xmax, root, stat, guess)

    class(brent_root), intent(inout) :: this
    real(r8), intent(in) :: xmin, xmax
    real(r8), intent(inout) :: root
    integer, intent(out) :: stat
    logical, intent(in), optional :: guess

    real(r8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,tol2,xm
    integer :: i
    logical :: guess_

    guess_ = .false.
    if (present(guess)) guess_ = guess

    a = xmin
    b = xmax
    c = xmax
    fa = this%f(a)
    fb = this%f(b)
    fc = fb

    ! If guess is specified as true by the calling code, the argument 'root' is expected
    ! to carry the initial guess when this function is called. At the end of the function,
    ! it gets modified to the final solution obtained by Brent's iteration

    if (guess_) then
      if (root < xmax .and. root > xmin) then
        b = root
        fb = this%f(b)
        e = 0.5_r8 * (c-b)
        d = e
      end if
    end if

    stat = 0

    root=huge(1.0_r8)
    do i = 1,this%maxitr
      if ((fb > 0.0_r8 .and. fc > 0.0_r8) .or. (fb < 0.0_r8 .and. fc < 0.0_r8)) then
        c = a
        fc = fa
        e = b-a
        d = e
      end if

      if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
      end if

      tol1 = 2.0_r8*epsilon(1.0_r8)*abs(b) + 0.5_r8*this%eps
      tol2 = 2.0_r8*epsilon(1.0_r8)*abs(fb) + 0.5_r8*this%feps
      xm = 0.5_r8*(c-b)

      if (abs(xm) <= tol1 .or. abs(fb) <= tol2) then
        root=b
        return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (a==c) then
          p = 2.0_r8 * xm * s
          q = 1.0_r8 - s
        else
          q = fa / fc
          r = fb / fc
          p = s*(2.0_r8*xm*q*(q-r)-(b-a)*(r-1.0_r8))
          q = (q-1.0_r8)*(r-1.0_r8)*(s-1.0_r8)
        end if
        if (p > 0.0_r8) q = -q
        p = abs(p)
        if (2.0_r8*p < min(3.0_r8*xm*q-abs(tol1*q), abs(e*q))) then
          e=d
          d=p/q
        else
          d=xm
          e=d
        end if
      else
        d=xm
        e=d
      end if
      a=b
      fa=fb
      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1,xm)
      end if
      fb = this%f(b)
    end do

    root = b
    stat = 1

  end subroutine find_root

end module brent_root_class
