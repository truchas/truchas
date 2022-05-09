!!
!! SIMULATION_EVENT_QUEUE
!!
!! This module implements an "event queue" that holds a time-ordered queue of
!! simulation events.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2015; expanded April 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module simulation_event_queue

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sim_event_queue_type
  use parameter_list_type
  use toolpath_event_type
  use toolhead_event_type, only: toolhead_event
  use diffusion_solver, only: vf_event
  implicit none
  private

  public :: read_simulation_control_namelist
  public :: init_sim_event_queue

  integer, parameter, public :: DT_POLICY_NONE   = 0
  integer, parameter, public :: DT_POLICY_NEXT   = 1
  integer, parameter, public :: DT_POLICY_FACTOR = 2
  integer, parameter, public :: DT_POLICY_VALUE  = 3

  type, extends(event_action), public :: phase_event
    private
    integer  :: dt_policy = DT_POLICY_NONE
    real(r8) :: c
  contains
    procedure :: init_dt => phase_event_init_dt
  end type phase_event

  type, extends(event_action), public :: output_event
#ifdef INTEL_BUG20180222
    integer :: dummy = 1
#endif
  end type

  type, extends(event_action), public :: short_edit_event
#ifdef INTEL_BUG20180222
    integer :: dummy = 1
#endif
  end type

  type, extends(event_action), public :: stop_event
#ifdef INTEL_BUG20180222
    integer :: dummy = 1
#endif
  end type

  type(sim_event_queue), public :: event_queue
  public :: toolhead_event, vf_event, toolpath_event

  type(parameter_list), public :: params

contains

  subroutine init_sim_event_queue(dt_min)

    use toolhead_driver, only: add_toolhead_events
    use diffusion_solver_data, only: ds_enabled
    use diffusion_solver, only: add_moving_vf_events
    use edit_module, only: short_edit, short_output_dt_multiplier
    use output_control

    real(r8), intent(in) :: dt_min

    integer :: i, j, n, dt_policy
    real(r8) :: c, t
    real(r8), allocatable :: array(:)

    call event_queue%set_time_resolution(dt_min)

    !! The order in which events are added can be significant when events
    !! occur at distinct but nearly the same time (within resolution).

    !! Events involving toolpaths are added first to avoid their times
    !! being adjusted to avoid small spacings.
    if (ds_enabled) call add_moving_vf_events(event_queue)
    call add_toolhead_events(event_queue)
    call add_output_events(event_queue) ! only for part_path right now

    !! Add user-specified phase start times
    if (params%is_parameter('phase-start-times')) then
      call params%get('phase-start-times', array)
      if (params%is_parameter('phase-init-dt')) then
        call params%get('phase-init-dt', c)
        dt_policy = DT_POLICY_VALUE
      else
        call params%get('phase-init-dt-factor', c)
        dt_policy = DT_POLICY_FACTOR
      end if
      do j = 1, size(array)
        call event_queue%add_event(array(j), phase_event(dt_policy, c))
      end do
    end if

    !! Add output times
    do j = 1, nops
      n = (output_t(j+1) - output_t(j) + 0.9*output_dt(j)) / output_dt(j)
      do i = 0, max(0, n-1)
        t = output_t(j) + i*output_dt(j)
        if (modulo(i,output_dt_multiplier(j)) == 0) then
          call event_queue%add_event(t, output_event())
        end if
        if (short_output_dt_multiplier(j) > 0) then
          if (modulo(i,short_output_dt_multiplier(j)) == 0) &
              call event_queue%add_event(t, short_edit_event())
        end if
      end do
    end do
    t = output_t(nops+1)
    call event_queue%add_event(t, output_event())
    if (short_edit) call event_queue%add_event(t, short_edit_event())
    call event_queue%add_event(t, stop_event(), rank=99)

  end subroutine init_sim_event_queue

  pure function phase_event_init_dt(this, dt_last, dt_next) result(dt)
    class(phase_event), intent(in) :: this
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
  end function phase_event_init_dt


  subroutine read_simulation_control_namelist (lun)

    use kinds, only: r8
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I
    use sort_module, only: sort
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios
    logical :: found
    character(128) :: iom
    real(r8), allocatable :: array(:)

    integer :: event_lookahead
    real(r8) :: phase_init_dt, phase_init_dt_factor, phase_start_times(500)
    namelist /simulation_control/ phase_start_times, phase_init_dt, phase_init_dt_factor, &
        event_lookahead

    !! Locate the optional SIMULATION_CONTROL namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'SIMULATION_CONTROL', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) return

    call TLS_info('Reading SIMULATION_CONTROL namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      phase_init_dt = NULL_R
      phase_init_dt_factor = NULL_R
      phase_start_times = NULL_R
      event_lookahead = NULL_I
      read(lun,nml=simulation_control,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading SIMULATION_CONTROL namelist: ' // trim(iom))

    !! Broadcast the namelist variables.
    call broadcast(phase_init_dt)
    call broadcast(phase_init_dt_factor)
    call broadcast(phase_start_times)
    call broadcast(event_lookahead)

    !! Check the variables.
    array = pack(phase_start_times, mask=(phase_start_times /= NULL_R))
    if (size(array) > 0) then
      call sort(array) !NB: should check for, or cull, duplicates
      call params%set('phase-start-times', array)
      if (phase_init_dt == NULL_R .eqv. phase_init_dt_factor == NULL_R) then
        call TLS_fatal('Either PHASE_INIT_DT or PHASE_INIT_DT_FACTOR must be defined but not both')
      else if (phase_init_dt /= NULL_R) then
        if (phase_init_dt <= 0.0_r8) call TLS_fatal('PHASE_INIT_DT must be > 0')
        call params%set('phase-init-dt', phase_init_dt)
      else if (phase_init_dt_factor /= NULL_R) then
        if (phase_init_dt_factor <= 0.0_r8) call TLS_fatal('PHASE_INIT_DT_FACTOR must be > 0')
        call params%set('phase-init-dt-factor', phase_init_dt_factor)
      end if
    end if

    if (event_lookahead /= NULL_I) then
      if (event_lookahead < 2) then
        call TLS_fatal('EVENT_LOOKAHEAD must be >= 2')
      end if
      call params%set('event-lookahead', event_lookahead)
    end if

  end subroutine read_simulation_control_namelist

end module simulation_event_queue
