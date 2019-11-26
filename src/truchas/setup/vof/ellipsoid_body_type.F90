!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ellipsoid_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: ellipsoid_body
    private
    real(r8), allocatable :: center(:), axes(:)
  contains
    procedure :: eval
    procedure :: signed_distance
  end type ellipsoid_body

  interface ellipsoid_body
    procedure ellipsoid_body_value
  end interface ellipsoid_body

contains

  !! constructor for ELLIPSOID_BODY objects
  function ellipsoid_body_value(xc, axes) result(r)
    real(r8), intent(in) :: xc(:), axes(:)
    type(ellipsoid_body) :: r
    ASSERT(size(xc)==size(axes))
    r%center = xc
    r%axes = axes
  end function ellipsoid_body_value

  logical function eval(this, x, cellid)
    class(ellipsoid_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = norm2((x - this%center) / this%axes) < 1
  end function eval

  real(r8) function signed_distance(this, x)
    class(ellipsoid_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = -sqrt(product(this%axes)) * (1 - norm2((x - this%center) / this%axes))
  end function signed_distance

end module ellipsoid_body_type
