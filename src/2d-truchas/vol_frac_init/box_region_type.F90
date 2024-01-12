module box_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region2d_class
  implicit none

  public :: alloc_box_region

  type, extends(region2d) :: box_region
    private
    real(r8) :: lower(2), upper(2)
    logical  :: outside = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_box_region(this, lower, upper, outside)
    class(region2d), allocatable, intent(out) :: this
    real(r8), intent(in) :: lower(:), upper(:)
    logical, intent(in), optional :: outside
    allocate(this, source=box_region(lower, upper, outside))
  end subroutine

  pure logical function encloses(this, x, bitmask)
    class(box_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask  ! irrelevant to this type
    encloses = all(x >= this%lower .and. x <= this%upper) .xor. this%outside
  end function

end module box_region_type
