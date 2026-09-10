!!
!! T2D_FLOW_THERMAL_INTEGRATOR_TYPE
!!
!! This module defines T2D_FLOW_THERMAL_INTEGRATOR, the time-integration
!! layer for coupled two-dimensional incompressible flow and thermal
!! transport. It selects endpoint times, manages recoverable thermal-step
!! failures, and delegates each coupled step to T2D_FLOW_THERMAL_SOLVER.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flow_thermal_integrator_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use t2d_flow_model_type
  use t2d_thermal_model_type
  use t2d_flow_thermal_solver_type
  use time_step_sync_type
  implicit none
  private

  type, public :: t2d_flow_thermal_integrator
    private
    type(t2d_flow_thermal_solver) :: solver
    real(r8) :: tlast, hlast, hnext
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    integer :: max_try
    character(8) :: hnext_cause = 'init'
    logical :: time_stepper_initialized = .false.
    type(time_step_sync) :: ts_sync
  contains
    procedure :: init
    procedure :: init_time_stepper
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: last_time
    procedure :: num_steps
    procedure :: initial_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: get_cell_flow_active
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
  end type

contains

  subroutine init(this, env, flow_model, ht_model, matl_model, matl_dist, params, stat, errmsg)

    class(t2d_flow_thermal_integrator), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_flow_model), target, intent(inout) :: flow_model
    type(t2d_thermal_model), target, intent(in) :: ht_model
    type(material_model), intent(in) :: matl_model
    type(material_distribution), target, intent(in) :: matl_dist
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%solver%init(env, flow_model, ht_model, matl_model, matl_dist, params, stat, errmsg)
    if (stat /= 0) return
  end subroutine


  !! Initialize the coupled time-step policy from TIME-STEPPING parameters.
  !! The output schedule itself remains owned by the simulation driver.
  subroutine init_time_stepper(this, params, stat, errmsg)

    class(t2d_flow_thermal_integrator), intent(inout) :: this
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
    call params%get('time-step-growth', this%dt_grow, default=1.05_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-lookahead', lookahead, default=3, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. lookahead < 1) then
      stat = 1
      errmsg = 'invalid coupled time-step controls, including time-step-lookahead >= 1'
      return
    end if
    if (this%max_try <= 0) then
      stat = 1
      errmsg = 'maximum coupled step attempts must be > 0'
      return
    end if
    this%ts_sync = time_step_sync(lookahead)
    this%time_stepper_initialized = .true.
  end subroutine


  subroutine set_initial_state(this, env, matl_model, time, velocity, temp, stat, errmsg)

    class(t2d_flow_thermal_integrator), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: time, velocity(:,:), temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8) :: hlimit

    ASSERT(this%time_stepper_initialized)
    call this%solver%set_initial_state(env, matl_model, time, this%dt_init, velocity, temp, stat, errmsg)
    if (stat /= 0) return
    this%tlast = time
    this%hnext = this%dt_init
    this%hnext_cause = 'init'
    hlimit = this%solver%courant_time_step()
    if (hlimit < this%hnext) then
      this%hnext = hlimit
      this%hnext_cause = 'cfl'
    end if
    this%hlast = this%hnext
  end subroutine


  subroutine integrate(this, env, matl_model, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(t2d_flow_thermal_integrator), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, thermal_hnext, hproposed
    logical :: sig_rcvd
    integer :: n
    character(256) :: line
    character(8) :: cause, attempt_cause

    stat = 0
    ASSERT(this%time_stepper_initialized)
    t_n = this%tlast
    ASSERT(tout >= t_n)
    do while (t_n < tout)
      hproposed = this%hnext
      cause = this%hnext_cause
      t_np1 = this%ts_sync%next_time(tout, t_n, this%hlast, this%hnext)
      if (t_np1 < t_n + hproposed) cause = 'output'
      do n = 1, this%max_try
        attempt_cause = cause
        if (n > 1) attempt_cause = 'thermal'
        write(line,'(a,i0,a,i0,a,es0.5,a,es0.5,a,a)') 'step=', this%solver%num_steps() + 1_int64, &
            ' attempt=', n, ' t0=', t_n, ' dt=', t_np1 - t_n, ' cause=', trim(attempt_cause)
        call env%simlog%begin_section(trim(line))
        call this%solver%step(env, matl_model, t_n, t_np1, stat, errmsg, thermal_hnext)
        if (stat == 0) then
          call env%simlog%end_section('step-end status=accepted')
          exit
        else if (stat > 0) then
          t_np1 = t_n + thermal_hnext
          if (t_np1 - t_n < this%dt_min) then
            stat = -1
            errmsg = 'next coupled time step is too small'
            call env%simlog%end_section('step-end status=failed')
            return
          end if
          if (n == this%max_try) then
            stat = -2
            errmsg = 'unable to take a coupled time step'
            call env%simlog%end_section('step-end status=failed')
            return
          end if
          call env%simlog%end_section('step-end status=rejected')
        else
          call env%simlog%end_section('step-end status=failed')
          return
        end if
      end do
      this%hlast = t_np1 - t_n
      call select_step_cause(this, thermal_hnext, this%hnext, this%hnext_cause)
      t_n = t_np1
      this%tlast = t_n
      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if
    end do
  end subroutine


  subroutine select_step_cause(this, thermal_hnext, hnext, cause)
    class(t2d_flow_thermal_integrator), intent(in) :: this
    real(r8), intent(in) :: thermal_hnext
    real(r8), intent(out) :: hnext
    character(*), intent(out) :: cause

    real(r8) :: hlimit

    hnext = thermal_hnext
    cause = 'thermal'
    hlimit = this%dt_grow*this%hlast
    if (hlimit < hnext) then
      hnext = hlimit
      cause = 'growth'
    end if
    if (this%dt_max < hnext) then
      hnext = this%dt_max
      cause = 'max'
    end if
    hlimit = this%solver%courant_time_step()
    if (hlimit < hnext) then
      hnext = hlimit
      cause = 'cfl'
    end if
  end subroutine


  real(r8) function last_time(this)
    class(t2d_flow_thermal_integrator), intent(in) :: this

    last_time = this%tlast
  end function


  integer(int64) function num_steps(this)
    class(t2d_flow_thermal_integrator), intent(in) :: this

    num_steps = this%solver%num_steps()
  end function


  real(r8) function initial_time_step(this)
    class(t2d_flow_thermal_integrator), intent(in) :: this

    ASSERT(this%time_stepper_initialized)
    initial_time_step = this%dt_init
  end function


  subroutine init_temporal_output(this, data)
    class(t2d_flow_thermal_integrator), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call this%solver%init_temporal_output(data)
  end subroutine


  subroutine set_temporal_output(this, data)
    class(t2d_flow_thermal_integrator), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call this%solver%set_temporal_output(data)
  end subroutine


  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(t2d_flow_thermal_integrator), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%solver%get_cell_flow_soln(pressure, velocity)
  end subroutine


  subroutine get_face_velocity(this, velocity)
    class(t2d_flow_thermal_integrator), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%solver%get_face_velocity(velocity)
  end subroutine


  subroutine get_cell_flow_active(this, active)
    class(t2d_flow_thermal_integrator), target, intent(in) :: this
    logical, pointer, intent(out) :: active(:)

    call this%solver%get_cell_flow_active(active)
  end subroutine


  subroutine get_cell_heat_soln(this, enth)
    class(t2d_flow_thermal_integrator), intent(in) :: this
    real(r8), intent(inout) :: enth(:)

    call this%solver%get_cell_heat_soln(enth)
  end subroutine


  subroutine get_cell_temp_soln(this, temp)
    class(t2d_flow_thermal_integrator), intent(in) :: this
    real(r8), intent(inout) :: temp(:)

    call this%solver%get_cell_temp_soln(temp)
  end subroutine

end module t2d_flow_thermal_integrator_type
