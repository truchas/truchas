!!
!! TIME_STEP_SYNC_TYPE
!!
!! This module defines a derived type that provides a method for gradually
!! adjusting the time steps so that one coincides with a given target time.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines the TIME_STEP_SYNC derived type.  It has a single type
!! bound function:
!!
!!  NEXT_TIME(TSTAR, T0, H0, H1) returns the time T1 for the next step.  TSTAR
!!    is the target time; T0 is the time of the current step, T0 < TSTAR; H0 is
!!    the last time step size; and, H1 is the next proposed time step size.
!!    All arguments are intent-in.  Nominally the returned value for T1 would
!!    be T0 + H1, however the function returns an adjusted time whenever the
!!    current time is within range of the target time.  When used consecutively
!!    to define the time for the next step, the function ensures that the
!!    target time is hit precisely while keeping the ratio of consecutive step
!!    sizes close to one.  The resulting step size will never exceed H1 in any
!!    event.  How gradually the step sizes are adjusted depends on how far in
!!    advance of the target time the times are adjusted. This is the look-ahead
!!    range, which is measured in terms of the number of steps given a current
!!    step size.  Note that the function does not modify H1; if the returned
!!    T1 is used for the next step, then H1 should be overwritten with T1 - T0.
!!
!! A variable of this type is initialized by assignment from a structure
!! constructor.  The constructor TIME_STEP_SYNC(N) returns an object with
!! look-ahead range set to N steps.
!!
!! NB: NEXT_TIME is carefully implemented so that it will return precisely
!! the target value TSTAR when intended and not merely something very close
!! (within a few ulps, for example).  As a result it is entirely proper to
!! use '==' when comparing the result to TSTAR.
!!
!! To give a feel for the behavior as a function of the range N, consider an
!! extreme example with a step size H and TSTAR - T0 = (N-1)*H + epsilon.
!! In other words step N-1 would fall just shy of the target time by a small
!! amount and so require N steps to reach TSTAR.  The function will reduce
!! the effective step size by a factor R<1 each step so that TSTAR is reached
!! in N steps.  Here are the values of R for various ranges.  Note that if
!! step N-1 misses TSTAR by a wider margin the associated R is closer to 1.
!!
!!     N |  R
!!   ----+------
!!     1 | 0.00
!!     2 | 0.50
!!     3 | 0.62
!!     4 | 0.81
!!     5 | 0.89
!!     6 | 0.93
!!     7 | 0.95
!!     8 | 0.96
!!     9 | 0.97
!!    10 | 0.98
!!
!! The result for N=1, effectively no look-ahead, is indicative of the bad
!! thing that can happen in the naive approach of simply cutting the step size
!! at the last moment so that TSTAR is not passed -- the time step can be
!! arbitrarily small.
!!

#include "f90_assert.fpp"

module time_step_sync_type

  use kinds, only: r8
  implicit none
  private

  type, public :: time_step_sync
    private
    integer :: n = 2
  contains
    procedure :: next_time
  end type time_step_sync

  !! Defined constructor
  interface time_step_sync
    procedure time_step_sync_init
  end interface

contains

  !! Defined constructor
  pure function time_step_sync_init (n) result (obj)
    integer, intent(in) :: n
    type(time_step_sync) :: obj
    obj%n = n
  end function time_step_sync_init

  function next_time (this, tstar, t0, h0, h1) result (t1)

    class(time_step_sync), intent(in) :: this
    real(r8), intent(in) :: tstar ! target time to hit exactly
    real(r8), intent(in) :: t0    ! current time
    real(r8), intent(in) :: h0    ! last time step
    real(r8), intent(in) :: h1    ! next time step (proposed)
    real(r8)             :: t1    ! next time

    integer :: n0, n1
    real(r8) :: a0, a1, r

    ASSERT(tstar > t0)

    a0 = min((tstar - t0)/h0, this%n + 1.0_r8)
    n0 = ceiling(a0)

    a1 = min((tstar - t0)/h1, this%n + 1.0_r8)
    n1 = ceiling(a1)

    if (n0 <= this%n) then ! capture for soft landing
      r = braking_ratio(n0,a0)
      if (h0*r <= h1) then
        t1 = merge(tstar, t0 + h0*r, n0==1)
      else if (n1 <= this%n) then
        r = braking_ratio(n1,a1)
        t1 = merge(tstar, t0 + h1*r, n1==1)
      else
        t1 = t0 + h1
      end if
    else if (n1 <= this%n) then
      r = braking_ratio(n1,a1)
      t1 = merge(tstar, t0 + h1*r, n1==1)
    else
      t1 = t0 + h1
    end if

  end function next_time

  !! Auxiliary function returns the r solving r + r^2 + ... + r^n = a,
  !! where n-1 < a <= n.

  pure function braking_ratio (n, a) result (r)
    integer, intent(in) :: n
    real(r8), intent(in) :: a
    real(r8) :: r, f, s
    integer :: k
    select case (n)
    case (1)
      r = a
    case (2)
      r = (sqrt(1 + 4*a) - 1) / 2
    case (3:) ! use Newton iteration
      r = 1 + 2*(a-n)/(n*(n+1))
      do  ! until converged to machine precision
        f = r
        do k = 2, n
          f = r*(1 + f)
        end do
        if (abs(f-a) <= 2*spacing(f)) exit
        s = n
        do k = n-1, 1, -1
          s = k + r*s
        end do
        r = r + (a-f)/s
      end do
    end select
  end function braking_ratio

end module time_step_sync_type
