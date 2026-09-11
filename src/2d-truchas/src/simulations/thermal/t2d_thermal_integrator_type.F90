!!
!! T2D_THERMAL_INTEGRATOR_TYPE
!!
!! This module defines the target-time integration layer for 2D thermal
!! transport.  It owns the adaptive time-step policy, output-time
!! synchronization, retry handling, and integration counters, and delegates
!! individual thermal steps to T2D_THERMAL_SOLVER.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_thermal_integrator_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use parameter_list_type
  use simulation_environment_type
  use signal_handler, only: read_signal, SIGUSR1
  use time_step_sync_type
  use t2d_thermal_model_type
  use t2d_thermal_solver_type
  implicit none
  private

  type, public :: t2d_thermal_integrator
    private
    type(t2d_thermal_solver), pointer :: solver => null()
    integer(int64) :: nstep = 0_int64
    real(r8) :: tlast, hlast, hnext
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    integer :: max_try
    character(8) :: hnext_cause = 'init'
    type(time_step_sync) :: ts_sync
    logical :: time_stepper_initialized = .false.
  contains
    final :: delete
    procedure :: init
    procedure :: init_time_stepper
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: last_time
    procedure :: num_steps
    procedure :: initial_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
  end type

contains

  subroutine init(this, env, model, params, stat, errmsg)

    class(t2d_thermal_integrator), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_thermal_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    allocate(this%solver)
    call this%solver%init(env, model, params, stat, errmsg)
  end subroutine init


  subroutine delete(this)
    type(t2d_thermal_integrator), intent(inout) :: this
    if (associated(this%solver)) deallocate(this%solver)
  end subroutine delete


  subroutine init_time_stepper(this, params, stat, errmsg)

    class(t2d_thermal_integrator), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: lookahead

    call params%get('initial-time-step', this%dt_init, stat, errmsg)
    if (stat /= 0) return
    call params%get('min-time-step', this%dt_min, stat, errmsg)
    if (stat /= 0) return
    call params%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-growth', this%dt_grow, default=5.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-lookahead', lookahead, default=3, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. lookahead < 1 .or. this%max_try <= 0) then
      stat = 1
      errmsg = 'require 0 < min-time-step <= initial-time-step <= max-time-step, ' // &
          'time-step-growth >= 1, time-step-lookahead >= 1, and max-try-at-step > 0'
      return
    end if
    this%ts_sync = time_step_sync(lookahead)
    this%time_stepper_initialized = .true.
    stat = 0
  end subroutine init_time_stepper


  subroutine set_initial_state(this, env, t, temp, stat, errmsg)

    class(t2d_thermal_integrator), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t, temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    ASSERT(this%time_stepper_initialized)
    call this%solver%set_initial_state(env, t, this%dt_init, temp, stat, errmsg)
    if (stat /= 0) return
    this%nstep = 0_int64
    this%tlast = t
    this%hnext = this%dt_init
    this%hnext_cause = 'init'
    this%hlast = this%hnext
  end subroutine set_initial_state


  subroutine integrate(this, env, tout, stat, errmsg)

    class(t2d_thermal_integrator), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, hproposed, hthermal
    logical :: sig_rcvd
    integer :: n
    character(256) :: line
    character(8) :: cause

    stat = 0
    ASSERT(this%time_stepper_initialized)
    t_n = this%tlast
    ASSERT(tout >= t_n)
    do while (t_n < tout)
      hproposed = this%hnext
      cause = this%hnext_cause
      t_np1 = this%ts_sync%next_time(tout, t_n, this%hlast, this%hnext)
      if (t_np1 < t_n + hproposed) cause = 'output'
      if (t_np1 - t_n < this%dt_min) then
        stat = -1
        errmsg = 'next time step is too small'
        return
      end if
      do n = 1, this%max_try
        write(line,'(a,i0,a,i0,a,es0.5,a,es0.5,a,a)') 'step=', this%nstep + 1_int64, &
            ' attempt=', n, ' t0=', t_n, ' dt=', t_np1 - t_n, ' cause=', trim(cause)
        call env%simlog%begin_section(trim(line))
        call this%solver%step(env, t_n, t_np1, stat, errmsg, hthermal)
        if (stat == 0) then
          call env%simlog%end_section('step-end status=accepted')
          exit
        end if
        t_np1 = t_n + hthermal
        if (t_np1 - t_n < this%dt_min) then
          stat = -1
          errmsg = 'next time step is too small'
          call env%simlog%end_section('step-end status=failed')
          return
        end if
        if (n == this%max_try) then
          stat = -2
          errmsg = 'unable to take a thermal time step'
          call env%simlog%end_section('step-end status=failed')
          return
        end if
        call env%simlog%end_section('step-end status=rejected')
      end do
      call this%solver%commit_step()
      t_n = t_np1
      this%nstep = this%nstep + 1_int64
      this%hlast = t_n - this%tlast
      this%tlast = t_n
      this%hnext = min(hthermal, this%dt_grow*this%hlast, this%dt_max)
      this%hnext_cause = 'thermal'
      if (this%hnext == this%dt_grow*this%hlast) this%hnext_cause = 'growth'
      if (this%hnext == this%dt_max) this%hnext_cause = 'max'

      call read_signal(SIGUSR1, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGUSR1 signal'
        return
      end if
    end do
  end subroutine integrate


  real(r8) function last_time(this)
    class(t2d_thermal_integrator), intent(in) :: this
    last_time = this%solver%last_time()
  end function last_time


  integer(int64) function num_steps(this)
    class(t2d_thermal_integrator), intent(in) :: this
    num_steps = this%nstep
  end function num_steps


  real(r8) function initial_time_step(this)
    class(t2d_thermal_integrator), intent(in) :: this
    ASSERT(this%time_stepper_initialized)
    initial_time_step = this%dt_init
  end function initial_time_step


  subroutine init_temporal_output(this, data)
    class(t2d_thermal_integrator), intent(in) :: this
    type(parameter_list), intent(inout) :: data
    call data%set('NStep', this%nstep)
  end subroutine init_temporal_output


  subroutine set_temporal_output(this, data)
    class(t2d_thermal_integrator), intent(in) :: this
    type(parameter_list), intent(inout) :: data
    call data%set('NStep', this%nstep)
  end subroutine set_temporal_output


  subroutine get_cell_heat_soln(this, enth)
    class(t2d_thermal_integrator), intent(in) :: this
    real(r8), intent(inout) :: enth(:)
    call this%solver%get_cell_heat_soln(enth)
  end subroutine get_cell_heat_soln


  subroutine get_cell_temp_soln(this, temp)
    class(t2d_thermal_integrator), intent(in) :: this
    real(r8), intent(inout) :: temp(:)
    call this%solver%get_cell_temp_soln(temp)
  end subroutine get_cell_temp_soln


end module t2d_thermal_integrator_type
