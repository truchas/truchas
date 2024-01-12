module cell_set_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region2d_class
  implicit none

  public :: alloc_cell_set_region

  type, extends(region2d) :: cell_set_region
    private
    integer :: bitmask
    logical :: outside = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_cell_set_region(this, bitmask, outside)
    class(region2d), allocatable, intent(out) :: this
    integer, intent(in) :: bitmask
    integer, intent(in), optional :: outside
    allocate(this, source=cell_set_region(bitmask, outside))
  end subroutine

  pure logical function encloses(this, x, bitmask)
    class(cell_set_region), intent(in) :: this
    real(r8), intent(in) :: x(:)  ! irrelevant to this type
    integer,  intent(in) :: bitmask
    encloses = (iand(bitmask, this%bitmask) /= 0) .xor. this%outside
  end function

end module cell_set_region_type
