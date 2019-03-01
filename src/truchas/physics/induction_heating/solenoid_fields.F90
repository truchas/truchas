!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solenoid_fields

  use kinds, only: r8
  implicit none
  private
  
  public :: H_coil
  
  !! Don't use; not ready for prime time!
  public :: H_sheet
  
contains

  !!
  !! This routine computes the H-field at a point due to current flow in a
  !! cylindrical, z-axial, n-turn coil centered at the origin.  The coil is
  !! approximated by an array of coaxial current loops centered at equi-spaced
  !! points along the z-axis from z=-a to z=a.  The result assumes a unit
  !! current in each loop that flows ccw with respect to coil's axis; the
  !! H-field is linear in the current.  If n is negative, the result corresponds
  !! to a |n|-turn coil with a cw flowing current.  If n=0, the result is a
  !! zero vector.
  !!
  
  function H_coil (x, radius, a, n) result (H)
  
    real(kind=r8), intent(in) :: x(3)
    real(kind=r8), intent(in) :: radius
    real(kind=r8), intent(in) :: a
    integer,       intent(in) :: n
    
    integer :: j
    real(kind=r8) :: H(3), dz, y(3)
    
    select case (abs(n))
    case (0)  ! No loops.
      H = 0.0_r8
    case (1)  ! Just a single loop centered at the origin.
      H = H_loop (x, radius)
      if (n < 0) H = -H
    case (2:) ! Multiple loops distributed over [-a,a] on the z-axis.
      dz = 2.0_r8 * a / real(abs(n)-1, kind=r8)
      y = x
      y(3) = y(3) + a
      H = 0.0_r8
      do j = 1, abs(n)
        H = H + H_loop(y, radius)
        y(3) = y(3) - dz
      end do
      if (n < 0) H = -H
    end select
    
  end function H_coil
  
  !!
  !! This routine computes the H-field at a point due to a circular current
  !! loop in the x/y-plane centered at the origin.  The result assumes a
  !! unit current flowing ccw; to obtain the result for a current J, simply
  !! multiply the result by J.
  !!
  !! NB: The general expression is indeterminate on the loop axis (divide by
  !! zero).  The routine makes a half-assed effort to remove this numerical
  !! singularity by using the analytic limiting expression for points exactly
  !! on the axis.  We really ought to be using a series expansion for points
  !! near the axis, with a smooth transition between the two expressions.
  !!
  
  function H_loop (x, radius) result (H)
  
    use elliptic_integrals, only: elk, ele
    real(kind=r8), intent(in) :: x(3), radius
    
    real(kind=r8) :: H(3), a, b, c, r, q, E, K, hz, hr, m
    real(kind=r8), parameter :: TWOPI = 6.2831853071795864769_r8
    
    r = sqrt(x(1)**2 + x(2)**2)
    if (r > 0.0_r8) then
      a = r / radius
      b = x(3) / radius
      c = x(3) / r
      q = (1.0_r8 + a)**2 + b**2
      m = (4.0_r8*a/q) ! elliptic parameter m
      E = ele(m)
      K = elk(m)
      hz =     (E * (1.0_r8 - a**2 - b**2)/(q - 4.0_r8*a) + K) / (TWOPI * radius * sqrt(q))
      hr = c * (E * (1.0_r8 + a**2 + b**2)/(q - 4.0_r8*a) - K) / (TWOPI * radius * sqrt(q))
      H(1) = hr * (x(1)/r)
      H(2) = hr * (x(2)/r)
      H(3) = hz
    else
      H(1) = 0.0_r8
      H(2) = 0.0_r8
      H(3) = 0.5_r8 * radius**2 / (radius**2 + x(3)**2)**1.5_r8
    end if
    
  end function H_loop
  
  !!
  !! This routine computes the H-field at a point due to a current in a
  !! cylindrical sheet solenoid.  This should be viewed as the limiting
  !! case of an n-turn coil as n -> infinity.  The routine uses the
  !! adaptive Romberg integration procedure.
  !!
  !! NB: This is still in the evaluation phase.
  !!
  
  function H_sheet (x, radius, a) result (H)
  
    real(kind=r8), intent(in) :: x(:), radius, a
    
    integer, parameter :: KMAX=10
    real(kind=r8) :: H(3), y1(3), y2(3), RTable(3,0:KMAX,0:KMAX)
    integer :: j, k
    
    y1 = x
    y2 = x
    y1(3) = y1(3) + a
    y2(3) = y2(3) - a
    RTable(:,0,0) = 0.5_r8 * (H_loop(y1, radius) + H_loop(y2, radius))
    print *, 0, RTable(:,0,0)
    
    do j = 1, KMAX
      RTable(:,j,0) = 0.5_r8 * (RTable(:,j-1,0) + H_coil(x, radius, a, 2**(j-1)))
      do k = 1, j
        RTable(:,j,k) = (4.0_r8**k * RTable(:,j,k-1) - RTable(:,j-1,k-1)) / (4.0_r8**k - 1.0_r8)
      end do
      print *, j, RTable(:,j,j)
      H = RTable(:,j,j)
      if (sqrt(sum((H-RTable(:,j-1,j-1))**2)) < 1.0e-4_r8 * sqrt(sum(H**2))) exit
    end do

  end function H_sheet
  
end module solenoid_fields
