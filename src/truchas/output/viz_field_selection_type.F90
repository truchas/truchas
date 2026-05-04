!! Shared representation of a provider-local visualization field selection.
!!
!! VTKHDF providers receive this value after manager-level name resolution.  A
!! selection may mean "write all provider defaults" or "write this selected
!! name list". Package-local VTKHDF writers use it to build
!! small output plans so temporal datasets are registered only for fields that
!! will actually be written.
module viz_field_selection_type

  implicit none
  private

  integer, parameter, public :: FIELD_MODE_ALL      = 1
  integer, parameter, public :: FIELD_MODE_SELECTED = 2

  type, public :: viz_field_selection
    integer :: mode = FIELD_MODE_ALL
    character(:), allocatable :: names(:)
  contains
    procedure :: init
    procedure :: write_all
    procedure :: has_selected_fields
  end type

contains

  subroutine init(this, field_names)

    class(viz_field_selection), intent(out) :: this
    character(*), intent(in), optional :: field_names(:)

    if (present(field_names)) then
      this%names = field_names
      this%mode = FIELD_MODE_SELECTED
    end if

  end subroutine init

  logical function write_all(this)
    class(viz_field_selection), intent(in) :: this
    write_all = (this%mode == FIELD_MODE_ALL)
  end function

  logical function has_selected_fields(this)
    class(viz_field_selection), intent(in) :: this
    has_selected_fields = (this%mode == FIELD_MODE_SELECTED) .and. allocated(this%names)
    if (has_selected_fields) has_selected_fields = has_selected_fields .and. size(this%names) > 0
  end function

end module viz_field_selection_type
