module half_plane_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_class
  implicit none

  public :: alloc_half_plane_region

  type, extends(region) :: half_plane_region
    private
    real(r8) :: point(2), normal(2)
    logical  :: complement = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_half_plane_region(this, point, normal, complement)
    class(region), allocatable, intent(out) :: this
    real(r8), intent(in) :: point(:), normal(:)
    logical, intent(in), optional :: complement
    allocate(this, source=half_plane_region(point, normal, complement))
  end subroutine

  pure logical function encloses(this, x, bitmask)
    class(half_plane_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask  ! unused for this type
    encloses = (dot_product(this%normal, x-this%point) <= 0) .neqv. this%complement
  end function

end module half_plane_region_type
