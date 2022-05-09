!!
!! TOOLPATH_EVENT_TYPE
!!
!! This module defines the derived type TOOLPATH_EVENT which extends the
!! abstract base class EVENT_ACTION. Events of this type provide an action
!! for advancing a toolpath from one segment to the next.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module toolpath_event_type

  use toolpath_type
  use sim_event_queue_type, only: event_action
  implicit none
  private

  type, extends(event_action), public :: toolpath_event
    private
    type(toolpath), pointer :: tp => null()
  contains
    procedure :: next_toolpath_segment
  end type

  !! User-defined constructor
  interface toolpath_event
    procedure tp_event
  end interface

contains

  !! Defined constructor for PATH_EVENT objects.
  function tp_event(tp) result(event)
    type(toolpath), intent(in), target :: tp
    type(toolpath_event) :: event
    event%tp => tp
  end function

  subroutine next_toolpath_segment(this)
    class(toolpath_event), intent(in) :: this
    call this%tp%next_segment
  end subroutine

end module toolpath_event_type
