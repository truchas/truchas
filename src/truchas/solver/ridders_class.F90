!!
!! RIDDERS_CLASS
!!
!! An abstract base class that implements Ridders method for solving f(x) = 0.
!! Concrete extensions implement the deferred type bound function f.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  Ridders' method is a robust root finding algorithm based on the false
!!  position method with performance similar to Brent's method.  The root is
!!  bracketed, and the length of the bracketing interval at least halves each
!!  iteration.  To apply the method to f(x) = 0, define an extension of the
!!  RIDDERS_CLASS type that defines the deferred type bound function F.  That
!!  function implements f(x) using the interface specified below.
!!
!!  To initialize an instance of this extended type, define the components
!!
!!    %EPS -- error tolerance for the root
!!    %MAXITR -- maximum number of iterations allowed
!!
!!  Once initialized the following method is called to compute the root.
!!
!!  FIND_ROOT(XMIN, XMAX, ROOT, STAT) computes an approximate value ROOT for
!!  a root of f(x).  The specified interval [XMIN, XMAX] must bracket a root
!!  and satisfy f(XMIN)*f(XMAX) <= 0.  The error in ROOT is at most the
!!  specified EPS.  STAT returns 0 if no errors were encountered.  If XMIN
!!  and XMAX do not bracket a root, STAT returns -1.  Otherwise, if the
!!  iteration fails to converge to the specified error tolerance within
!!  MAXITR iterations, STAT returns 1.
!!
!!  The following (read-only) components return info about the FIND_ROOT call:
!!
!!    %NUMITR -- the number of iterations taken
!!    %ERROR  -- an upper bound on the error in ROOT; literally the length
!!               of the final bracketing interval
!!
!!  If the function is strictly monotone, the following method can be used to
!!  determine a bracketing interval required by FIND_ROOT.
!!
!!  BRACKET_ROOT(XMIN, XMAX, MAXADJ, STAT) attempts to discover an interval
!!  that brackets the root. It assumes the function is strictly monotone, and
!!  the interval is successively shifted to the left or right, its length
!!  doubled with each attempt until the interval brackets a root.  MAXADJ
!!  specifies the maximum number of adjustments allowed. The input values
!!  of XMIN and XMAX are an initial guess for the interval, and their output
!!  values are the discovered bracketing interval if successful (STAT==0).
!!
!!  The read-only %NUMADJ component returns the number of interval adjustments
!!  made by the last BRACKET_ROOT call.
!!
!! IMPLEMENTATION NOTE
!!
!!  Why not make find_root a simple procedure that takes a function f(x) as
!!  one of its arguments?  The problem with that approach is that typical
!!  functions aren't static, hardcoded things, but depend on parameters, often
!!  in complicated ways.  The problem is how to supply this additional context
!!  data to the function.  The chosen method of using an abstract base class
!!  is an elegant solution to that problem; the context data is added as data
!!  components of the extended type.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ridders_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: ridders
    real(r8) :: eps     ! convergence tolerance
    integer  :: maxitr  ! maximum number of iterations allowed
    real(r8) :: error = 0.0_r8  ! estimate of the error in the root
    integer  :: numitr = 0      ! number of iterations taken
    integer  :: numadj = 0      ! number of bracketing interval adjustments
  contains
    procedure, non_overridable :: find_root
    procedure :: bracket_root
    procedure(func), deferred  :: f
  end type ridders

  abstract interface
    function func (this, x) result (fx)
      import ridders, r8
      class(ridders), intent(in) :: this
      real(r8), intent(in) :: x
      real(r8) :: fx
    end function
  end interface

contains

  subroutine find_root (this, xmin, xmax, root, stat)

    class(ridders), intent(inout) :: this
    real(r8), intent(in) :: xmin, xmax
    real(r8), intent(out) :: root
    integer, intent(out) :: stat

    real(r8) :: a, b, c, m, fa, fb, fc, fm

    ASSERT(xmin <= xmax)
    a = xmin
    b = xmax

    this%error = 0.0_r8
    this%numitr = 0
    stat = 0

    fa = this%f(a)
    if (fa == 0.0_r8) then
      root = a
      return
    end if

    fb = this%f(b)
    if (fb == 0.0_r8) then
      root = b
      return
    end if

    if (fa == sign(fa,fb)) then ! interval doesn't bracket a root
      stat = -1
      return
    end if

    do  ! until converged
      if (this%numitr == this%maxitr) then
        stat = 1
        return
      end if
      this%numitr = this%numitr + 1
      !! Next approximate root C.
      m = 0.5d0*(a + b)
      fm = this%f(m)
      c = m + (m-a) * fm * sign(1.0d0/sqrt(fm*fm - fa*fb), fa)
      !! First convergence check.
      if (this%numitr > 1) then
        this%error = abs(c - root)
        if (this%error <= this%eps) then
          root = c
          return
        end if
      end if
      root = c
      !! Update the interval bracketing the root.
      fc = this%f(c)
      if (fc == 0.0_r8) then
        this%error = 0.0_r8
        return
      else if (fm /= sign(fm,fa)) then  ! [a,m]
        if (fm /= sign(fm,fc)) then ! [c,m]
          a = c; fa = fc
          b = m; fb = fm
        else ! [a,c]
          b = c; fb = fc
        end if
      else  ! [m,b]
        if (fm /= sign(fm,fc)) then ! [m,c]
          a = m; fa = fm
          b = c; fb = fc
        else  ! [c,b]
          a = c; fa = fc
        end if
      endif
      ASSERT(a < b)
      ASSERT(fa*fb < 0.0_r8)
      !! Second convergence check: [a,b] brackets the root.
      this%error = b - a
      if (this%error <= this%eps) return
    end do

  end subroutine find_root

  subroutine bracket_root(this, xmin, xmax, maxadj, stat)

    class(ridders), intent(inout) :: this
    real(r8), intent(inout) :: xmin, xmax
    integer,  intent(in) :: maxadj
    integer,  intent(out) :: stat

    integer :: n
    real(r8) :: a, b, d, fa, fb

    ASSERT(xmin < xmax)

    stat = 0
    this%numadj = 0
    a = xmin; b = xmax
    fa = this%f(a); fb = this%f(b)
    if (fa /= sign(fa,fb)) return ! initial interval brackets root

    d = 2*(b - a)
    if (fa == sign(fa,fb-fa)) then  ! shift interval left
      do n = 1, maxadj
        b = a; fb = fa
        a = b - d; fa = this%f(a)
        if (fa /= sign(fa,fb)) exit ! [a,b] brackets root
        d = 2*d
      end do
    else  ! shift interval right
      do n = 1, maxadj
        a = b; fa = fb
        b = a + d; fb = this%f(b)
        if (fb /= sign(fb,fa)) exit ! [a,b] brackets root
        d = 2*d
      end do
    end if

    xmin = a; xmax = b
    if (n <= maxadj) then
      this%numadj = n
      stat = 0
    else
      stat = -1
    end if

  end subroutine bracket_root

end module ridders_class
