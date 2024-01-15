!FIXME: This mimics the existing "body" implementation but specialized for 2D.
! The term "body" is replaced by "region" to prevent module name clashes
! (Though I prefer region to body for this. I think associating a material with
! a region forms a body in the modeling sense.)  However there is very little
! different between 2D and 3D, and the two should try to be merged somehow.

module region_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: region
  contains
    procedure(encloses), deferred :: encloses
    procedure :: ifunc
  end type

  abstract interface
    pure logical function encloses(this, x, bitmask)
      import region, r8
      class(region), intent(in) :: this
      real(r8), intent(in) :: x(:)
      integer, intent(in) :: bitmask
    end function
  end interface

contains

  pure integer function ifunc(this, x, bitmask)
    class(region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask
    ifunc = merge(1, 0, this%encloses(x, bitmask))
  end function

end module region_class
