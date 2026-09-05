module t2d_box_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use t2d_region_class
  implicit none

  public :: alloc_box_region

  type, extends(t2d_region) :: t2d_box_region
    private
    real(r8) :: lower(2), upper(2)
    logical  :: complement = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_box_region(this, lower, upper, complement)
    class(t2d_region), allocatable, intent(out) :: this
    real(r8), intent(in) :: lower(:), upper(:)
    logical, intent(in), optional :: complement
    allocate(this, source=t2d_box_region(lower, upper, complement))
  end subroutine

  pure logical function encloses(this, x, bitmask)
    class(t2d_box_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask  ! unused for this type
    encloses = all(x >= this%lower .and. x <= this%upper) .neqv. this%complement
  end function

end module t2d_box_region_type
