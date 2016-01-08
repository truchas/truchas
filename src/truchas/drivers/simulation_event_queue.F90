!!
!! SIMULATION_EVENT_QUEUE
!!
!! This module implements a rudimentary "event queue" that holds a time-ordered
!! queue of simulation events.  This is currently limited to the starts of
!! simulation phases, but anticipates a comprehensive system handling all types
!! of simulation events currently handled elsewhere in an ad hoc manner.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module holds a time-based "event queue" as private data. In the present
!!  implementation the events are limited to the start of simulation phases.
!!  This data is initialized through namelist input via a call to the following
!!  subroutine:
!!
!!  READ_SIMULATION_CONTROL_NAMELIST (LUN) reads the SIMULATION_CONTROL
!!    namelist from the file open on unit LUN.  The namelist comprises the
!!    real array PHASE_START_TIMES and real scalars PHASE_INIT_DT and
!!    PHASE_INIT_DT_FACTOR.  PHASE_START_TIMES is an array of times marking
!!    the starting time of one or more simulation phases.  The remaining
!!    variables specify the initial time step size for each simulation phase:
!!    PHASE_INIT_DT gives the absolute step size, and PHASE_INIT_DT_FACTOR the
!!    ratio of the step size to the final step size of the preceding phase.
!!
!!  After the module is initialized, the function NEXT_EVENT(T) returns the next
!!  event, in order, with time > T.  The returned value is of type SIM_EVENT
!!  which has the following type bound functions:
!!
!!  TIME() returns the time of the event
!!
!!  INIT_DT(DT_LAST, DT_NEXT) returns the initial step size.  DT_LAST is the
!!    preceding step size (concluding the preceding phase), and DT_NEXT is
!!    the next proposed time step (determined by the client).  The returned
!!    value may, or may not, depend on these values according to the internal
!!    configuration of the event.
!!
!! IMPLEMENTATION NOTES
!!
!!  While this module was written to meet immediate needs (a capability for
!!  multi-phase simulations), it is the kernel of a larger design concept to
!!  address a more general problem.  The simulation process is punctuated
!!  by numerous events besides the start of a simulation phase. These include
!!  output of field data and checkpoint data, in situ result analysis, etc.
!!  The concept of an "event queue" aims to unify these diverse actions as
!!  specific types of events, ordered in time, with the simulation proceding
!!  from one event to the next.
!!

module simulation_event_queue

  use kinds, only: r8
  use parameter_list_type
  implicit none
  private

  public :: read_simulation_control_namelist, next_event

  type, public :: sim_event
    private
    real(r8) :: t, c
    integer ::  dt_policy
  contains
    procedure :: time => sim_event_time
    procedure :: init_dt => sim_event_init_dt
  end type sim_event

  integer, parameter :: DT_POLICY_NEXT   = 1
  integer, parameter :: DT_POLICY_FACTOR = 2
  integer, parameter :: DT_POLICY_VALUE  = 3

  type(sim_event), allocatable, target :: event_queue(:)
  integer :: top

contains

  pure function sim_event_time (this) result (t)
    class(sim_event), intent(in) :: this
    real(r8) :: t
    t = this%t
  end function sim_event_time

  pure function sim_event_init_dt (this, dt_last, dt_next) result (dt)
    class(sim_event), intent(in) :: this
    real(r8), intent(in) :: dt_last, dt_next
    real(r8) :: dt
    select case (this%dt_policy)
    case (DT_POLICY_NEXT)
      dt = dt_next
    case (DT_POLICY_FACTOR)
      dt = this%c*dt_last
    case (DT_POLICY_VALUE)
      dt = this%c
    case default
      dt = 0.0_r8 ! should never be here
    end select
  end function sim_event_init_dt

  function next_event (t) result (event)
    real(r8), intent(in), optional :: t
    type(sim_event), pointer :: event
    integer :: j
    event => null()
    if (.not.allocated(event_queue)) return
    if (present(t)) then
      do j = top, size(event_queue)
        if (t < event_queue(j)%time()) exit
      end do
      top = j
    end if
    if (top > size(event_queue)) return
    event => event_queue(top)
    top = top + 1
  end function next_event

  subroutine read_simulation_control_namelist (lun)

    use kinds, only: r8
    use input_utilities, only: seek_to_namelist, NULL_R
    use sort_module, only: sort
    use string_utilities, only: i_to_c, raise_case
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios
    logical :: found
    real(r8), allocatable :: array(:)

    real(r8) :: phase_init_dt, phase_init_dt_factor, phase_start_times(500)
    namelist /simulation_control/ phase_start_times, phase_init_dt, phase_init_dt_factor

    !! Locate the optional SIMULATION_CONTROL namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'SIMULATION_CONTROL', found, iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast (found)
    if (.not.found) return

    call TLS_info ('')
    call TLS_info ('Reading SIMULATION_CONTROL namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      phase_init_dt = NULL_R
      phase_init_dt_factor = NULL_R
      phase_start_times = NULL_R
      read(lun,nml=simulation_control,iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading SIMULATION_CONTROL namelist')

    !! Broadcast the namelist variables.
    call broadcast (phase_init_dt)
    call broadcast (phase_init_dt_factor)
    call broadcast (phase_start_times)

    !! Check the variables.

    allocate(array(count(phase_start_times /= NULL_R)))
    array = pack(phase_start_times, mask=(phase_start_times /= NULL_R))
    if (size(array) == 0) call TLS_fatal ('no values assigned to PHASE_START_TIMES')
    call sort (array) !NB: should check for, or cull, duplicates
    !call params%set ('phase-start-times', array)
    allocate(event_queue(size(array)))
    event_queue%t = array
    top = 1

    if (phase_init_dt == NULL_R .eqv. phase_init_dt_factor == NULL_R) then
      call TLS_fatal ('Either PHASE_INIT_DT or PHASE_INIT_DT_FACTOR must be defined but not both')
    else if (phase_init_dt /= NULL_R) then
      if (phase_init_dt <= 0.0_r8) call TLS_fatal ('PHASE_INIT_DT must be > 0')
      !call params%set ('phase-init-dt', phase_init_dt)
      event_queue%dt_policy = DT_POLICY_VALUE
      event_queue%c = phase_init_dt
    else if (phase_init_dt_factor /= NULL_R) then
      if (phase_init_dt_factor <= 0.0_r8) call TLS_fatal ('PHASE_INIT_DT_FACTOR must be > 0')
      !call params%set ('phase-init-dt-factor', phase_init_dt_factor)
      event_queue%dt_policy = DT_POLICY_FACTOR
      event_queue%c = phase_init_dt_factor
    end if

  end subroutine read_simulation_control_namelist

end module simulation_event_queue
