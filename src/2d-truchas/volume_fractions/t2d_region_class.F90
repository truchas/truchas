!FIXME: This mimics the existing "body" implementation but specialized for 2D.
! The term "body" is replaced by "region" to prevent module name clashes
! (Though I prefer region to body for this. I think associating a material with
! a region forms a body in the modeling sense.)  However there is very little
! different between 2D and 3D, and the two should try to be merged somehow.

module t2d_region_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: t2d_region
  contains
    procedure(encloses), deferred :: encloses
    procedure :: ifunc
  end type

  abstract interface
    pure logical function encloses(this, x, bitmask)
      import t2d_region, r8
      class(t2d_region), intent(in) :: this
      real(r8), intent(in) :: x(:)
      integer, intent(in) :: bitmask
    end function
  end interface

  type, public :: t2d_region_box
    class(t2d_region), allocatable :: reg
  contains
    procedure :: encloses => region_box_encloses
    procedure :: ifunc => region_box_ifunc
  end type

  !! A special region that contains all points/cells
  type, extends(t2d_region) :: t2d_background_region
  contains
    procedure :: encloses => background_encloses
  end type

  public :: alloc_background_region

contains

  pure integer function ifunc(this, x, bitmask)
    class(t2d_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask
    ifunc = merge(1, 0, this%encloses(x, bitmask))
  end function

  pure logical function region_box_encloses(this, x, bitmask)
    class(t2d_region_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask
    region_box_encloses = this%reg%encloses(x, bitmask)
  end function

  pure integer function region_box_ifunc(this, x, bitmask)
    class(t2d_region_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask
    region_box_ifunc = this%reg%ifunc(x, bitmask)
  end function

  subroutine alloc_background_region(this)
    class(t2d_region), allocatable, intent(out) :: this
    allocate(t2d_background_region :: this)
  end subroutine

  pure logical function background_encloses(this, x, bitmask)
    class(t2d_background_region), intent(in) :: this
    real(r8), intent(in) :: x(:)    ! unused for this type
    integer, intent(in) :: bitmask  ! unused for this type
    background_encloses = .true.
  end function
  
end module t2d_region_class
