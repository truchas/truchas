!!
!!  ELLIPTIC_INTEGRALS
!!
!!    Neil N. Carlson <nnc@newmexico.com> 29 Aug 2003
!!
!!  This module provides functions to compute the complete elliptic integrals
!!  of the first and second kind, and Carlson's (no relation!) elliptic
!!  integrals of the first and second kind.  The code is based on Carlson's
!!  algorithms (1979), code by Carlson and Notis (1981), and the GSL (Gnu
!!  Scientific Library) implementation.  The module adopts the interface used
!!  by the IMSL special functions library.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  NB: The routines ELRF and ELRD have not been thoroughly verified.  The
!!  complete elliptic integrals (which use ELRF and ELRD) have been generally
!!  verified, but not thoroughly so in the asymptotic regimes.
!!
!!  ROUTINES:
!!  The following functions return a real result of the same kind as the
!!  argument.  At present only double precision versions are defined.
!!  The argument X to the complete elliptic integral functions is the elliptic
!!  'parameter m'; other standard forms use the elliptic 'modulus k', m = k^2,
!!  or the 'modular angle alpha', m = sin(alpha)^2.
!!
!!    ELK(X) returns the complete elliptic integral of the first kind K(x).
!!      The argument X must satisfy 0 <= X < 1.
!!
!!    ELE(X) returns the complete elliptic integral of the second kind E(x).
!!      The argument X must satisfy 0 <= X < 1.
!!
!!    ELRF(X,Y,Z) returns Carlson's elliptic function of the first kind
!!      R_F(x,y,z).  The arguments must be nonnegative.
!!
!!    ELRD(X,Y,Z) returns Carlson's elliptic function of the second kind
!!      R_D(x,y,z).  The arguments must be nonnegative.
!!
!!  REFERENCES:
!!  B. C. Carlson, Computing elliptic integrals by duplication,
!!    Numer. Math. 33 (1979).
!!  B. C. Carlson and E. M. Notis, Algorithms for incomplete elliptic
!!    integrals, ACM Transactions on Mathematical Software, September (1981).
!!

module elliptic_integrals

  implicit none
  private
  
  public :: elk, ele, elrf, elrd
  
  interface elk
    module procedure d_elk
  end interface

  interface ele
    module procedure d_ele
  end interface

  interface elrf
    module procedure d_elrf
  end interface
  
  interface elrd
    module procedure d_elrd
  end interface
  
contains

  !!
  !!  COMPLETE ELLIPTIC INTEGRALS OF THE FIRST AND SECOND KINDS.
  !!
  !!  These are computed via Carlson's symmetric elliptic integrals.  The GSL implementation
  !!  checks for domain error, and uses a separate asymptotic expansion (Abramowitz-Stegun)
  !!  in the case the elliptic parameter x (conventionally written 'm') is close to 1.
  !!  I won't bother with either:  I'll let the symmetric integral functions handle domain
  !!  errors, and their relative accuracy should be equivalent in the asymptotic regime,
  !!  though perhaps at greater expense.
  !!
  
  double precision function d_elk (x)
    double precision, intent(in) :: x
    d_elk = elrf(0.0d0, 1.0d0-x, 1.0d0)
  end function d_elk
  
  double precision function d_ele (x)
    double precision, intent(in) :: x
    d_ele = elrf(0.0d0, 1.d0-x, 1.0d0) - (x/3.0d0) * elrd(0.0d0, 1.d0-x, 1.0d0)
  end function d_ele
  
  !!
  !!  CARLSON'S INCOMPLETE ELLIPTIC INTEGRAL OF THE FIRST KIND R_F(X,Y,Z)
  !!
  !!  According to Carlson, the relative error due to truncation is bounded by
  !!  ERRTOL**6/(4(1-ERRTOL)), and he gives the following table:
  !!
  !!                  ERRTOL   Rel Error
  !!                  1.D-3    3.D-19
  !!                  3.D-3    2.D-16
  !!                  1.D-2    3.D-13
  !!                  3.D-2    2.D-10
  !!                  1.D-1    3.D-7
  !!
  !! We set ERRTOL = 1.0D-3.
  !!
  
  double precision function d_elrf (x, y, z)
  
    double precision, intent(in) :: x, y, z
  
    double precision :: xn, yn, zn, xndev, yndev, zndev, xnroot, ynroot, znroot, mu, eps, lambda, e2, e3, s
    double precision, parameter :: LOLIM = 5.0d0 * tiny(1.0d0)
    double precision, parameter :: UPLIM = 0.2d0 * huge(1.0d0)
    double precision, parameter :: C1 = 1.d0/24.0d0, C2 = 3.0d0/44.0d0, C3 = 1.0d0/14.0d0
    double precision, parameter :: ERRTOL = 1.0d-3
  
    if (min(x,y,z)<0.0d0 .or. max(x,y,z)>UPLIM .or. min(x+y,x+z,y+z)<LOLIM) then
    
      d_elrf = d_elrf_raise_invalid()
      d_elrf = huge(1.0d0)  ! IMSL function behavior
      
    else
    
      xn = x
      yn = y
      zn = z
      do
        mu = (xn + yn + zn) / 3.0d0
        xndev = 2.0d0 - (mu + xn) / mu
        yndev = 2.0d0 - (mu + yn) / mu
        zndev = 2.0d0 - (mu + zn) / mu
        eps = max(abs(xndev), abs(yndev), abs(zndev))
        if (eps < ERRTOL) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lambda = xnroot * (ynroot + znroot) + ynroot * znroot
        xn = (xn + lambda) * 0.25d0
        yn = (yn + lambda) * 0.25d0
        zn = (zn + lambda) * 0.25d0
      end do
      e2 = xndev * yndev - zndev * zndev
      e3 = xndev * yndev * zndev
      s = 1.0d0 + (C1 * e2 - 0.1d0 - C2 * e3) * e2 + C3 * e3
      d_elrf = s / sqrt(mu)
      
    end if
    
  contains
  
    !!  We want to raise an ieee invalid exception, such as would be raised by sqrt(-1.0);
    !!  however we must be a bit obtuse in doing it, else a compiler may catch it at
    !!  compile time.
    
    function d_elrf_raise_invalid ()
      double precision :: d_elrf_raise_invalid
      double precision, save :: x = 1.0d0
      d_elrf_raise_invalid = sqrt(1.0d0 - 2.0d0 * x)
    end function d_elrf_raise_invalid
      
  end function d_elrf    

  !!
  !!  CARLSON'S INCOMPLETE ELLIPTIC INTEGRAL OF THE SECOND KIND R_D(X,Y,Z).
  !!
  !!  According to Carlson, the relative error due to truncation is bounded by
  !!  3*ERRTOL**6/(1-ERRTOL)**(3/2), and he gives the following table:
  !!
  !!                      ERRTOL   Rel Error
  !!                      1.D-3    4.D-18
  !!                      3.D-3    3.D-15
  !!                      1.D-2    4.D-12
  !!                      3.D-2    3.D-9
  !!                      1.D-1    4.D-6
  !!
  !!  We choose ERRTOL = 1.0D-3.  
  !!
  
  double precision function d_elrd (x, y, z)
  
    double precision, intent(in) :: x, y, z
  
    double precision :: xn, yn, zn, xndev, yndev, zndev, xnroot, ynroot, znroot
    double precision :: mu, lambda, eps, sigma, power4, ea, eb, ec, ed, ef, s1, s2
    double precision :: LOLIM, UPLIM
    
    double precision, parameter :: C1 = 3.0d0 / 14.0d0, C2 = 1.0d0 /  6.0d0
    double precision, parameter :: C3 = 9.0d0 / 22.0d0, C4 = 3.0d0 / 26.0d0
    double precision, parameter :: ERRTOL = 1.0d-3
  
    !! Because of the powers, we can't make these parameters.
    LOLIM = 2.0d0 / huge(1.0d0)**(2.0d0/3.0d0)
    UPLIM = (0.1d0 * ERRTOL / tiny(1.0d0))**(2.0d0/3.0d0)
    
    if (min(x,y)<0.0d0 .or. min(x+y,z)<LOLIM.or. max(x,y,z)>UPLIM) then
    
      d_elrd = d_elrd_raise_invalid()
      d_elrd = huge(1.0d0)  ! IMSL function behavior
      
    else
          
      xn = x
      yn = y
      zn = z
      sigma  = 0.0d0
      power4 = 1.0d0
      do
        mu = (xn + yn + 3.0d0 * zn) * 0.2d0
        xndev = (mu - xn) / mu
        yndev = (mu - yn) / mu
        zndev = (mu - zn) / mu
        eps = max(abs(xndev), abs(yndev), abs(zndev))
        if (eps < ERRTOL) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lambda = xnroot * (ynroot + znroot) + ynroot * znroot
        sigma  = sigma + power4 / (znroot * (zn + lambda))
        power4 = power4 * 0.25d0
        xn = (xn + lambda) * 0.25d0
        yn = (yn + lambda) * 0.25d0
        zn = (zn + lambda) * 0.25d0
      end do
      ea = xndev * yndev
      eb = zndev * zndev
      ec = ea - eb
      ed = ea - 6.0d0 * eb
      ef = ed + ec + ec
      s1 = ed * (- C1 + 0.25d0 * C3 * ed - 1.5d0 * C4 * zndev * ef)
      s2 = zndev * (C2 * ef + zndev * (-C3 * ec + zndev * C4 * ea))
      d_elrd = 3.0d0 * sigma + power4 * (1.0d0 + s1 + s2) / (mu * sqrt(mu))
      
    end if
    
  contains
  
    !! We want to raise an ieee invalid exception, such as would be raised by sqrt(-1.0);
    !! however we must be a bit obtuse in doing it, else a compiler may catch it at
    !! compile time.
    
    function d_elrd_raise_invalid ()
      double precision :: d_elrd_raise_invalid
      double precision, save :: x = 1.0d0
      d_elrd_raise_invalid = sqrt(1.0d0 - 2.0d0 * x)
    end function d_elrd_raise_invalid
      
  end function d_elrd    
  
end module elliptic_integrals
