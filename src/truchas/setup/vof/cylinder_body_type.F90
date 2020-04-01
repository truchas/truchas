!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module cylinder_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: cylinder_body
    private
    real(r8), allocatable :: center(:), axis(:)
    real(r8) :: radius, length
    logical :: fill_outside
  contains
    procedure :: eval
    procedure :: signed_distance
  end type cylinder_body

  interface cylinder_body
    procedure cylinder_body_value
  end interface cylinder_body

contains

  !! constructor for CYLINDER_BODY objects
  function cylinder_body_value(xc, axis, length, radius, fill_inside) result(r)
    real(r8), intent(in) :: xc(:), axis(:), length, radius
    logical, intent(in) :: fill_inside
    type(cylinder_body) :: r
    INSIST(norm2(axis) > 0)
    r%center = xc
    r%axis = axis / norm2(axis)
    r%length = length
    r%radius = radius
    r%fill_outside = .not.fill_inside
  end function cylinder_body_value

  logical function eval(this, x, cellid)

    class(cylinder_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid

    real(r8) :: xl, yl, x1(3), xo(3)

    ! get the local coordinates xl, yl
    xo = x - this%center
    x1 = dot_product(xo, this%axis) * this%axis
    yl = norm2(x1)
    xl = norm2(xo - x1)
    eval = yl <= this%length / 2 .and. xl <= this%radius
    if (this%fill_outside) eval = .not.eval

  end function eval

  real(r8) function signed_distance(this, x)
    class(cylinder_body), intent(in) :: this
    real(r8), intent(in) :: x(:)

    real(r8) :: xl, yl, x1(3), xo(3)

    ! get the local coordinates xl, yl
    xo = x - this%center
    x1 = dot_product(xo, this%axis) * this%axis
    yl = norm2(x1)
    xl = norm2(xo - x1)

    !signed_distance = 0 ! not implemented yet
    signed_distance = xl - this%radius
    if (this%fill_outside) signed_distance = -signed_distance

  end function signed_distance

end module cylinder_body_type
