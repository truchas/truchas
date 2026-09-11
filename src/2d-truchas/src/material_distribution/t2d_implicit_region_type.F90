module t2d_implicit_region_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  use t2d_region_class
  implicit none
  private

  public :: alloc_implicit_region

  type, extends(t2d_region) :: t2d_implicit_region
    private
    class(scalar_func), allocatable :: f
    logical :: complement = .false.
  contains
    procedure :: encloses
  end type

contains

  subroutine alloc_implicit_region(this, f, complement)
    class(t2d_region), allocatable, intent(out) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    logical, intent(in), optional :: complement

    allocate(t2d_implicit_region :: this)
    select type (this)
    type is (t2d_implicit_region)
      call move_alloc(f, this%f)
      if (present(complement)) this%complement = complement
    end select
  end subroutine

  logical function encloses(this, x, bitmask)
    class(t2d_implicit_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask  ! unused for this type

    encloses = (this%f%eval(x) <= 0.0_r8) .neqv. this%complement
  end function

end module t2d_implicit_region_type
