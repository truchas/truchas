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
    real(r8), allocatable :: upper(:), lower(:), xl(:), xc(:)
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
    r%xc = (upper + lower) / 2
    r%xl = (upper - lower) / 2
    ASSERT(all(r%xl > 0))
  end function box_body_value


  logical function eval(this, x, cellid)
    class(box_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    eval = all(x <= this%upper .and. x >= this%lower)
  end function eval


  real(r8) function signed_distance(this, x)

    class(box_body), intent(in) :: this
    real(r8), intent(in) :: x(:)

    integer :: i
    real(r8) :: q(size(x))

    q = abs(x - this%xc) - this%xl
    signed_distance = min(maxval(q), 0.0_r8)
    do i = 1, size(q)
      q(i) = max(q(i), 0.0_r8)
    end do
    signed_distance = signed_distance + norm2(q)

  end function signed_distance

end module box_body_type
