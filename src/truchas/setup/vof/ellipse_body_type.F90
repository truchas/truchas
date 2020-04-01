!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Defines an infinitely long elliptic cylinder aligned with the z-axis.

#include "f90_assert.fpp"

module ellipse_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: ellipse_body
    private
    real(r8), allocatable :: center(:), coeffs(:)
    logical :: fill_outside
  contains
    procedure :: eval
    procedure :: signed_distance
  end type ellipse_body

  interface ellipse_body
    procedure ellipse_body_value
  end interface ellipse_body

contains

  !! constructor for ELLIPSE_BODY objects
  function ellipse_body_value(xc, coeffs, fill_inside) result(r)
    real(r8), intent(in) :: xc(:), coeffs(:)
    logical, intent(in) :: fill_inside
    type(ellipse_body) :: r
    ASSERT(size(xc) >= 2 .and. size(coeffs) == 2)
    r%center = xc(:2)
    r%coeffs = coeffs
    r%fill_outside = .not.fill_inside
  end function ellipse_body_value

  logical function eval(this, x, cellid)
    class(ellipse_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = norm2((x(:2) - this%center) / this%coeffs) <= 1
    if (this%fill_outside) eval = .not.eval
  end function eval

  real(r8) function signed_distance(this, x)
    class(ellipse_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = sqrt(product(this%coeffs)) * (norm2((x(:2) - this%center) / this%coeffs) - 1)
    if (this%fill_outside) signed_distance = -signed_distance
  end function signed_distance

end module ellipse_body_type
