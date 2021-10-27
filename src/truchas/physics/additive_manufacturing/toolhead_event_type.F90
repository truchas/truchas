!! TOOLHEAD_EVENT_TYPE
!!
!! This module defines the derived type TOOLHEAD_EVENT which extends the
!! abstract base class EVENT_ACTION. Events of this type provide an action
!! for advancing the toolpath associated with toolhead to the next path
!! segment, and making any necessary updates to the toolhead sources that
!! result from switching path segments.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module toolhead_event_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sim_event_queue_type
  use toolhead_type
  implicit none
  private

  public :: add_toolhead_events

  integer, parameter :: DT_POLICY_NONE   = 0
  integer, parameter :: DT_POLICY_NEXT   = 1
  integer, parameter :: DT_POLICY_FACTOR = 2
  integer, parameter :: DT_POLICY_VALUE  = 3

  type, extends(event_action), public :: toolhead_event
    private
    type(toolhead), pointer :: th => null()
    integer  :: dt_policy = DT_POLICY_NONE
    real(r8) :: c = 0.0_r8
  contains
    procedure :: update_toolhead
    procedure :: init_dt
  end type

  !! User-defined constructor
  interface toolhead_event
    procedure th_event
  end interface

contains

  !! Constructor for TOOLHEAD_EVENT objects.
  function th_event(th, dt_policy, c) result(event)
    type(toolhead), intent(in), target :: th
    integer, intent(in), optional :: dt_policy
    real(r8), intent(in), optional :: c
    type(toolhead_event) :: event
    event%th => th
    if (present(dt_policy)) event%dt_policy = dt_policy
    if (present(c)) event%c = c
  end function

  subroutine update_toolhead(this, t)
    class(toolhead_event), intent(in) :: this
    real(r8), intent(in) :: t
    call this%th%next_tp_segment(t)
  end subroutine

  pure function init_dt(this, dt_last, dt_next) result(dt)
    class(toolhead_event), intent(in) :: this
    real(r8), intent(in) :: dt_last, dt_next
    real(r8) :: dt
    select case (this%dt_policy)
    case (DT_POLICY_NONE, DT_POLICY_NEXT)
      dt = dt_next
    case (DT_POLICY_FACTOR)
      dt = this%c*dt_last
    case (DT_POLICY_VALUE)
      dt = this%c
    case default
      dt = 0.0_r8 ! should never be here
    end select
  end function

  subroutine add_toolhead_events(th, eventq)

    type(toolhead), intent(in), target :: th
    type(sim_event_queue), intent(inout) :: eventq

    integer :: j
    real(r8), allocatable :: times(:)
    logical,  allocatable :: discont(:)

    call th%tp%get_segment_starts(times, discont)
    do j = 1, size(times)
      if (discont(j)) then
        call eventq%add_event(times(j), toolhead_event(th, DT_POLICY_FACTOR, 1.0e-1_r8), rank=80)
      else
        call eventq%add_event(times(j), toolhead_event(th, DT_POLICY_NEXT), rank=80)
      end if
    end do

  end subroutine

end module
