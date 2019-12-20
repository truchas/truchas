!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module box_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  implicit none
  private

  type, extends(body), public :: box_body
    private
    real(r8), allocatable :: upper(:), lower(:)
  contains
    procedure :: eval
    procedure :: signed_distance
  end type box_body

  interface box_body
    procedure box_body_value
  end interface box_body

contains

  !! constructor for BOX_BODY objects
  function box_body_value(upper, lower) result(r)
    real(r8), intent(in) :: upper(:), lower(:)
    type(box_body) :: r
    r%upper = upper
    r%lower = lower
  end function box_body_value

  logical function eval(this, x, cellid)
    class(box_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = all(x <= this%upper .and. x >= this%lower)
  end function eval

  ! Note this does not take into account corner distance
  real(r8) function signed_distance(this, x)
    class(box_body), intent(in) :: this
    real(r8), intent(in) :: x(:)

    integer :: i

    signed_distance = 0 ! not implemented yet
    ! signed_distance = huge(1.0_r8)
    ! do i = 1, size(x)
    !   signed_distance = minmag(signed_distance, x(i) - this%upper(i), this%lower(i) - x(i))
    ! end do

  end function signed_distance

  real(r8) pure function minmag(a, b, c)
    real(r8), intent(in) :: a, b, c
    real(r8) :: a2, b2, c2
    a2 = abs(a)
    b2 = abs(b)
    c2 = abs(c)
    minmag = merge(merge(a, c, a2 < c2), merge(b, c, b2 < c2), a2 < b2)
  end function minmag

end module box_body_type
