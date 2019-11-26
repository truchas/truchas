!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module plane_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: plane_body
    private
    real(r8), allocatable :: normal(:)
    real(r8) :: plane_const
  contains
    procedure :: eval
    procedure :: signed_distance
  end type plane_body

  interface plane_body
    procedure plane_body_value
  end interface plane_body

contains

  !! constructor for PLANE_BODY objects
  function plane_body_value(n, p) result(r)
    real(r8), intent(in) :: n(:), p
    type(plane_body) :: r
    r%normal = n
    r%plane_const = p
  end function plane_body_value

  logical function eval(this, x, cellid)
    class(plane_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = this%signed_distance(x) < 0
  end function eval

  real(r8) function signed_distance(this, x)
    class(plane_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = dot_product(x, this%normal) - this%plane_const
  end function signed_distance

end module plane_body_type
