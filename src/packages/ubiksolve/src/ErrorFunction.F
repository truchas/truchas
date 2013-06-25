! -*- Mode: f90 -*-

!*******************************************************************
! $Id: ErrorFunction.F,v 1.3 2002/12/03 05:33:21 turner Exp $
!
! <para>
! Use rational approximation to compute error function.
!
! erf x =
!    1 - (a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5) e^(-x^2) + eps
!
! where:
!               1
!        t = ---------  ,  ABS(eps) <= 1.5e-7 ,
!             1 + p*x
!
! and p, a1, a2, a3, a4, and a5 are constants.
!
! Reference:
!   Handbook of Mathematical Functions
!   M. Abramowitz and I. E. Stegun, ed.
!   Tenth Printing, 12/72, p. 299
! </para>
!
! <PARAMETER LIST>
!
!  Input:
!   x - argument
!
!  Output:
!   ErrorFunction - error function evaluated at x
!
!***************************************************************
module Ubik_ErrorFunction_module
  implicit none

  private
  public :: ErrorFunction

contains

  function ErrorFunction(x)
    use ubik_types_module

    ! <arguments>
    real(real_type), intent(in) :: x
    real(real_type) :: ErrorFunction
    ! </arguments>

    ! Local:
    real(real_type) :: t, t2, t3, t4
    real(real_type), parameter ::  &
         one = 1.0,  &
         p = 0.3275911,  &
         a1 =  0.254829592,  &
         a2 = -0.284496736,  &
         a3 =  1.421413741,  &
         a4 = -1.453152027,  &
         a5 =  1.061405429

    t = one / (one + p*x)
    t2 = t*t
    t3 = t*t2
    t4 = t*t3

    ErrorFunction = one - ( a1*t + a2*t2 + a3*t3  + a4*t4 + a5*t4*t )*EXP(-x*x)

    return
  end function ErrorFunction

end module Ubik_ErrorFunction_module
