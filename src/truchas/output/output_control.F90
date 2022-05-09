!!
!! OUTPUT_CONTROL
!!
!! This module was created to hold the variables from OUTPUT_DATA_MODULE and
!! OUTPUT_MODULE that control when the solution was output.  Everything else
!! in those modules is tied to TBrook and they are slated for deletion when
!! TBrook is ultimately deleted from the code.  I expect this module to exist
!! only temporarily until the control of output is refactored/reimplemented.
!!
!! Neil Carlson <nnc@lanl.gov>
!! October 2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module output_control

  use kinds, only: r8
  use parameter_module, only: mops
  use toolpath_type
  implicit none
  private

  public :: output_init, add_output_events

  real(r8), save, public :: output_dt(mops)
  real(r8), save, public :: output_t(mops+1)
  integer,  save, public :: output_dt_multiplier(mops) = 1
  integer,  save, public :: nops

  integer, allocatable, public :: part(:)
  type(toolpath), allocatable, target, public :: part_path
  logical, public :: write_mesh_partition

contains

  subroutine output_init(t)
    real(r8), intent(in) :: t
    if (allocated(part_path)) call part_path%set_segment(t)
  end subroutine

  subroutine add_output_events(eventq)
    use sim_event_queue_type
    use toolpath_event_type
    type(sim_event_queue), intent(inout) :: eventq
    integer :: j
    if (allocated(part_path)) then
      block
        real(r8), allocatable :: times(:)
        call part_path%get_segment_starts(times)
        do j = 1, size(times)
          call eventq%add_event(times(j), toolpath_event(part_path))
        end do
      end block
    end if
  end subroutine

end module output_control
