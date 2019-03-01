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
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
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

#include "f90_assert.fpp"

module simulation_event_queue

  use kinds, only: r8
  use parameter_list_type
  implicit none
  private

  public :: read_simulation_control_namelist, add_event, next_event

  integer, parameter, public :: DT_POLICY_NONE   = 0
  integer, parameter, public :: DT_POLICY_NEXT   = 1
  integer, parameter, public :: DT_POLICY_FACTOR = 2
  integer, parameter, public :: DT_POLICY_VALUE  = 3

  type, public :: sim_event
    private
    real(r8) :: t
    integer  :: dt_policy = DT_POLICY_NONE
    real(r8) :: c
  contains
    procedure :: time => sim_event_time
    procedure :: init_dt => sim_event_init_dt
  end type sim_event

  !! Defined constructor
  interface sim_event
    module procedure sim_event_values
  end interface

  type :: sim_event_queue
    type(queue_item), pointer :: top => null()
  contains
    procedure :: add
    procedure :: pop
    procedure :: is_empty
  end type sim_event_queue

  type :: queue_item
    type(sim_event) :: event
    type(queue_item), pointer :: next => null()
  end type queue_item

  type(sim_event_queue) :: event_queue

contains

  function sim_event_values(t, dt_policy, c) result(this)
    real(r8), intent(in) :: t
    integer, intent(in) :: dt_policy
    real(r8), intent(in), optional :: c
    type(sim_event) :: this
    this%t = t
    this%dt_policy = dt_policy
    if (present(c)) this%c = c
  end function sim_event_values

  subroutine add(this, event)
    use truchas_logging_services
    class(sim_event_queue), intent(inout) :: this
    type(sim_event), intent(in) :: event
    type(queue_item), pointer :: item, new_item
    item => find_queue_item(this, event%t)
    if (associated(item)) then
      if (item%event%t == event%t) then ! arbitrate on policy
        if (item%event%dt_policy == DT_POLICY_NONE) then
          item%event = event
        else
          call TLS_warn('SIMULATION_EVENT_QUEUE: conflicting DT policy')
        end if
      else ! insert after item
        allocate(new_item)
        new_item%event = event
        new_item%next => item%next
        item%next => new_item
      end if
    else ! insert at the front of the list
      allocate(new_item)
      new_item%event = event
      new_item%next => this%top
      this%top => new_item
    end if
  end subroutine add

  function pop(this) result(event)
    class(sim_event_queue), intent(inout) :: this
    type(sim_event) :: event
    type(queue_item), pointer :: item
    INSIST(associated(this%top))
    event = this%top%event
    item => this%top
    this%top => item%next
    deallocate(item)
  end function pop

  logical function is_empty(this)
    class(sim_event_queue), intent(in) :: this
    is_empty = .not.associated(this%top)
  end function is_empty

  !! Returns a pointer to the QUEUE_ITEM in the ordered list positioned at the
  !! given time.  This is the insertion point of an event with this time, but
  !! note that the event time of the item may equal the given time.  A null
  !! pointer is returned if the insertion point is at the head of the list.
  function find_queue_item(this, t) result(item)
    class(sim_event_queue), intent(in) :: this
    real(r8), intent(in) :: t
    type(queue_item), pointer :: item
    item => this%top
    if (.not.associated(item)) return
    if (t < item%event%t) then
      item => null()
      return
    end if
    do while (associated(item%next))
      if (t < item%next%event%t) return
      item => item%next
    end do
  end function find_queue_item

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
    case (DT_POLICY_NONE, DT_POLICY_NEXT)
      dt = dt_next
    case (DT_POLICY_FACTOR)
      dt = this%c*dt_last
    case (DT_POLICY_VALUE)
      dt = this%c
    case default
      dt = 0.0_r8 ! should never be here
    end select
  end function sim_event_init_dt

  subroutine add_event (event)
    type(sim_event), intent(in) :: event
    call event_queue%add(event)
  end subroutine add_event

  function next_event (t) result (event)
    real(r8), intent(in), optional :: t
    type(sim_event), pointer :: event
    event => null()
    if (event_queue%is_empty()) return
    allocate(event)
    event = event_queue%pop()
    if (present(t)) then
      do while (event%time() <= t)
        if (event_queue%is_empty()) then
          deallocate(event)
          return
        end if
        event = event_queue%pop()
      end do
    end if
  end function next_event

  subroutine read_simulation_control_namelist (lun)

    use kinds, only: r8
    use input_utilities, only: seek_to_namelist, NULL_R
    use sort_module, only: sort
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios, j, dt_policy
    logical :: found
    real(r8) :: c
    real(r8), allocatable :: array(:)
#ifdef INTEL_COMPILER_WORKAROUND
    type(sim_event) :: event
#endif

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

    if (phase_init_dt == NULL_R .eqv. phase_init_dt_factor == NULL_R) then
      call TLS_fatal ('Either PHASE_INIT_DT or PHASE_INIT_DT_FACTOR must be defined but not both')
    else if (phase_init_dt /= NULL_R) then
      if (phase_init_dt <= 0.0_r8) call TLS_fatal ('PHASE_INIT_DT must be > 0')
      !call params%set ('phase-init-dt', phase_init_dt)
      dt_policy = DT_POLICY_VALUE
      c = phase_init_dt
    else if (phase_init_dt_factor /= NULL_R) then
      if (phase_init_dt_factor <= 0.0_r8) call TLS_fatal ('PHASE_INIT_DT_FACTOR must be > 0')
      !call params%set ('phase-init-dt-factor', phase_init_dt_factor)
      dt_policy = DT_POLICY_FACTOR
      c = phase_init_dt_factor
    end if

    do j = size(array), 1, -1
#ifdef INTEL_COMPILER_WORKAROUND
      event = sim_event(array(j), dt_policy, c)
      call event_queue%add(event)
#else
      call event_queue%add(sim_event(array(j), dt_policy, c))
#endif
    end do

  end subroutine read_simulation_control_namelist

end module simulation_event_queue
