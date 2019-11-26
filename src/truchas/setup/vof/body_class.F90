!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module body_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: body
  contains
    procedure(point_is_inside_body), deferred :: eval
    procedure(signed_distance_to_surface), deferred :: signed_distance
  end type body

  type, public :: body_box
    class(body), allocatable :: f
  contains
    procedure :: eval
    procedure :: signed_distance
  end type body_box

  abstract interface
    logical function point_is_inside_body(this, x, cellid)
      import :: body, r8
      class(body), intent(in) :: this
      real(r8), intent(in) :: x(:)
      integer, intent(in) :: cellid
    end function point_is_inside_body

    real(r8) function signed_distance_to_surface(this, x)
      import :: body, r8
      class(body), intent(in) :: this
      real(r8), intent(in) :: x(:)
    end function signed_distance_to_surface
  end interface

contains

  logical function eval(this, x, cellid)
    class(body_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = this%f%eval(x, cellid)
  end function eval

  real(r8) function signed_distance(this, x)
    class(body_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = this%f%signed_distance(x)
  end function signed_distance

end module body_class
