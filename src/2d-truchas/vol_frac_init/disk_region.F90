module disk_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_class
  implicit none

  public :: alloc_disk_region

  type, extends(region) :: disk_region
    private
    real(r8) :: center(2), radius
    logical  :: complement = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_disk_region(this, center, radius, complement)
    class(region), allocatable, intent(out) :: this
    real(r8), intent(in) :: center(:), radius
    logical, intent(in), optional :: complement
    allocate(this, source=disk_region(center, radius, complement))
  end subroutine

  pure logical function encloses(this, x, bitmask)
    class(disk_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask  ! irrelevant to this type
    encloses = (norm2(x-this%center) <= this%radius) .xor. this%complement
  end function

end module disk_region_type
