!#include "f90_assert.fpp"

module region2d_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: region2d
  contains
    procedure(encloses), deferred :: encloses
    procedure :: ifunc
    !generic :: encloses => encloses_mesh, encloses_geom
    !procedure :: encloses_mesh, encloses_geom
  end type

  abstract interface
    pure logical function encloses(this, x, bitmask)
      import region2d, r8
      class(region2d), intent(in) :: this
      real(r8), intent(in) :: x(:)
      integer, intent(in) :: bitmask
    end function
    pure integer function ifunc(this, x, bitmask)
      import region2d, r8
      class(region2d), intent(in) :: this
      real(r8), intent(in) :: x(:)
      integer, intent(in) :: bitmask
    end function
  end interface

contains

  pure integer function ifunc(this, x, bitmask)
    class(region2d), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask
    ifunc = merge(1, 0, this%encloses(x, bitmask))
  end function

!contains
!
!  logical function encloses_mesh(this, bitmask)
!    class(region2d), intent(in) :: this
!    integer, intent(in) :: bitmask
!    INSIST(.false.)
!    encloses_mesh = .false.
!  end function
!
!  logical function encloses_geom(this, x)
!    class(region2d), intent(in) :: this
!    real(r8), intent(in) :: x(:)
!    INSIST(.false.)
!    encloses_geom = .false.
!  end function
!
end module region2d_class
