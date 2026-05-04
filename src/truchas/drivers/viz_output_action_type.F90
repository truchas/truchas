!! Event-queue action used to trigger VTKHDF output.
!!
!! viz_stream_state schedules these actions at the stream's requested output
!! times.  The existing driver dispatch path recognizes this action type and
!! calls the TVO output facade, which writes the selected stream timestep.
module viz_output_action_type

  use sim_event_queue_type, only: event_action
  implicit none
  private

  type, extends(event_action), public :: viz_output_action
    integer :: stream_id = 0
  end type

end module viz_output_action_type
