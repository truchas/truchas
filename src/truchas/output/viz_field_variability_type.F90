!! Variability and output-mode helpers for VTKHDF field registration.
!!
!! Providers classify fields by whether their values are constant over the run or
!! may change over time. The stream layer then combines that variability with
!! block output requirements to choose the concrete VTKHDF encoding.
module viz_field_variability_type

  implicit none
  private

  integer, parameter, public :: VIZ_FIELD_CONSTANT = 1
  integer, parameter, public :: VIZ_FIELD_DYNAMIC = 2

  integer, parameter, public :: VIZ_FIELD_OUTPUT_STATIC = 1
  integer, parameter, public :: VIZ_FIELD_OUTPUT_TEMPORAL_CONSTANT = 2
  integer, parameter, public :: VIZ_FIELD_OUTPUT_TEMPORAL_DYNAMIC = 3

  public :: get_viz_field_output_mode

contains

  pure integer function get_viz_field_output_mode(variability, requires_temporal_output)

    integer, intent(in) :: variability
    logical, intent(in) :: requires_temporal_output

    select case (variability)
    case (VIZ_FIELD_DYNAMIC)
      get_viz_field_output_mode = VIZ_FIELD_OUTPUT_TEMPORAL_DYNAMIC
    case default
      if (requires_temporal_output) then
        get_viz_field_output_mode = VIZ_FIELD_OUTPUT_TEMPORAL_CONSTANT
      else
        get_viz_field_output_mode = VIZ_FIELD_OUTPUT_STATIC
      end if
    end select

  end function get_viz_field_output_mode

end module viz_field_variability_type
