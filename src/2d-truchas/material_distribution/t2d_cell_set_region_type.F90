module t2d_cell_set_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use t2d_region_class
  implicit none

  public :: alloc_cell_set_region

  type, extends(t2d_region) :: t2d_cell_set_region
    private
    integer :: bitmask
    logical :: complement = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_cell_set_region(this, bitmask, complement)
    class(t2d_region), allocatable, intent(out) :: this
    integer, intent(in) :: bitmask
    logical, intent(in), optional :: complement
    allocate(this, source=t2d_cell_set_region(bitmask, complement))
  end subroutine

  pure logical function encloses(this, x, bitmask)
    class(t2d_cell_set_region), intent(in) :: this
    real(r8), intent(in) :: x(:)  ! unused for this type
    integer,  intent(in) :: bitmask
    encloses = (iand(bitmask, this%bitmask) /= 0) .neqv. this%complement
  end function

end module t2d_cell_set_region_type
