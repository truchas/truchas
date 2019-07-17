!!
!! plane_2d_type
!!
!! This module defines a 2d plane type (a line), along with routines for
!! calculating intersection points and distance.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module plane_2d_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  ! dot(n,x) - rho = 0
  type, public :: plane
    real(r8) :: rho, normal(2) ! plane constant and normal
  contains
    procedure :: signed_distance
    procedure :: intersects
    procedure :: intersection_point
    procedure :: print_data
  end type plane

  ! This cut-off for local precision comes from Hopcroft, J. E., & Kahn, P. J. (1992).
  ! A paradigm for robust geometric algorithms. Algorithmica, 7(1-6), 339-380.
  real(r8), parameter, public :: alpha = 1e-9_r8 ! local precision cutoff

contains

  ! calculates the signed distance from a plane
  real(r8) function signed_distance (this,x)

    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(:)

    ASSERT(size(x)==size(this%normal))

    signed_distance = dot_product(x,this%normal) - this%rho

    ! set distance to zero if it lies within alpha of the plane
    signed_distance = merge(signed_distance, 0.0_r8, abs(signed_distance) > alpha)
  end function signed_distance

  ! check if the plane lies between two points in space
  logical function intersects (this,x)

    use near_zero_function

    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(:,:) ! tuple of x positions

    real(r8) :: d1,d2

    ASSERT(size(x, dim=1)==size(this%normal))
    ASSERT(size(x, dim=2)==2)

    d1 = this%signed_distance(x(:,1))
    d2 = this%signed_distance(x(:,2))

    ! if the signed distances have opposite signs, the two points are on opposite sides of the plane
    intersects = sign(1.0_r8,d1)/=sign(1.0_r8,d2) .or. near_zero (d1,alpha) .or. near_zero (d2,alpha)
  end function intersects

  ! return the point where the line between x1 & x2 intersects with the given plane
  subroutine intersection_point (this, intx, on_point, x)

    use near_zero_function

    class(plane), intent(in) :: this
    real(r8), intent(out) :: intx(:)
    integer, intent(out) :: on_point
    real(r8), intent(in) :: x(:,:)

    real(r8) :: dx(2),d1,d2

    ASSERT(all(shape(x)==[2,2]))

    if (.not.this%intersects (x)) call TLS_panic('edge does not intersect plane')

    d1 = this%signed_distance(x(:,1))
    d2 = this%signed_distance(x(:,2))

    if (near_zero(d1,alpha)) then
      intx = x(:,1)
      on_point = 1
    else if (near_zero(d2,alpha)) then
      intx = x(:,2)
      on_point = 2
    else
      dx = x(:,2)-x(:,1)
      intx = x(:,1) - (d1/sum(dx*this%normal)) * dx
      on_point = 0
    end if

  end subroutine intersection_point

  subroutine print_data (this)
    class(plane), intent(in) :: this
    print '(a,4es30.20)', 'plane n, rho: ',this%normal, this%rho
  end subroutine print_data

end module plane_2d_type
