!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module sphere_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: sphere_body
    private
    real(r8), allocatable :: center(:)
    real(r8) :: radius
    logical :: fill_outside
  contains
    procedure :: eval
    procedure :: signed_distance
  end type sphere_body

  interface sphere_body
    procedure sphere_body_value
  end interface sphere_body

contains

  !! constructor for SPHERE_BODY objects
  function sphere_body_value(xc, radius, fill_inside) result(r)
    real(r8), intent(in) :: xc(:), radius
    logical, intent(in) :: fill_inside
    type(sphere_body) :: r
    r%center = xc
    r%radius = radius
    r%fill_outside = .not.fill_inside
  end function sphere_body_value

  logical function eval(this, x, cellid)
    class(sphere_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = norm2(x-this%center) <= this%radius
    if (this%fill_outside) eval = .not.eval
  end function eval

  real(r8) function signed_distance(this, x)
    class(sphere_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = norm2(x-this%center) - this%radius
    if (this%fill_outside) signed_distance = -signed_distance
  end function signed_distance

end module sphere_body_type
