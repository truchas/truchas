!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module background_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: background_body
  contains
    procedure :: eval
    procedure :: signed_distance
  end type background_body

contains

  logical function eval(this, x, cellid)
    class(background_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = .true.
  end function eval

  real(r8) function signed_distance(this, x)
    class(background_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = 0
  end function signed_distance

end module background_body_type
