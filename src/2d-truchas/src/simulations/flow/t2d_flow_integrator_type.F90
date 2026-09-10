!!
!! T2D_FLOW_INTEGRATOR_TYPE
!!
!! This module defines T2D_FLOW_INTEGRATOR, the isothermal incompressible
!! Navier--Stokes time-integration layer. It selects endpoint times and
!! delegates each complete flow step to T2D_FLOW_SOLVER. The flow solver owns
!! material transport and flow mechanics; this type owns only the standalone
!! time-step policy and target-time progression.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flow_integrator_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use t2d_flow_model_type
  use t2d_flow_solver_type
  use time_step_sync_type
  implicit none
  private

  type, public :: t2d_flow_integrator
    private
    type(t2d_flow_solver) :: solver
    real(r8) :: tlast, hlast, hnext
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    logical :: time_stepper_initialized = .false.
    type(time_step_sync) :: ts_sync
  contains
    procedure :: init
    procedure :: init_time_stepper
    procedure :: set_volume_fractions
    procedure :: set_initial_material_state
    procedure :: get_reduced_volume_fractions
    procedure :: update_material_distribution
    procedure :: set_buoyancy_temperature
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_cell_flow_active
    procedure :: get_face_velocity
    procedure :: integrate
    procedure :: last_time
    procedure :: num_steps
    procedure :: initial_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
  end type

contains

  subroutine init(this, env, model, matl_model, params, stat, errmsg)
    class(t2d_flow_integrator), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_flow_model), target, intent(inout) :: model
    type(material_model), intent(in) :: matl_model
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%solver%init(env, model, matl_model, params, stat, errmsg)
  end subroutine


  !! Initialize the standalone time-step policy from TIME-STEPPING parameters.
  !! The output schedule itself remains owned by the simulation driver.
  subroutine init_time_stepper(this, params, stat, errmsg)
    class(t2d_flow_integrator), intent(inout) :: this
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
    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. lookahead < 1) then
      stat = 1
      errmsg = 'require 0 < min-time-step <= initial-time-step <= max-time-step, ' // &
          'time-step-growth >= 1, and time-step-lookahead >= 1'
      return
    end if
    this%ts_sync = time_step_sync(lookahead)
    this%time_stepper_initialized = .true.
    stat = 0
    errmsg = ''
  end subroutine


  subroutine set_volume_fractions(this, vfrac)
    class(t2d_flow_integrator), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    call this%solver%set_volume_fractions(vfrac)
  end subroutine


  subroutine set_initial_material_state(this, vfrac, temperature)
    class(t2d_flow_integrator), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%solver%set_initial_material_state(vfrac, temperature)
  end subroutine


  subroutine get_reduced_volume_fractions(this, matl_dist, vfrac)
    class(t2d_flow_integrator), intent(in) :: this
    type(material_distribution), intent(in) :: matl_dist
    real(r8), allocatable, intent(out) :: vfrac(:,:)

    call this%solver%get_reduced_volume_fractions(matl_dist, vfrac)
  end subroutine


  subroutine update_material_distribution(this, matl_dist)
    class(t2d_flow_integrator), intent(in) :: this
    type(material_distribution), intent(inout) :: matl_dist

    call this%solver%update_material_distribution(matl_dist)
  end subroutine


  subroutine set_buoyancy_temperature(this, temperature)
    class(t2d_flow_integrator), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    call this%solver%set_buoyancy_temperature(temperature)
  end subroutine


  subroutine set_initial_state(this, env, time, dt, velocity, stat)
    class(t2d_flow_integrator), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    ASSERT(this%time_stepper_initialized)
    call this%solver%set_initial_state(env, time, dt, velocity, stat)
    if (stat /= 0) return
    this%tlast = time
    this%hnext = min(this%dt_init, this%solver%courant_time_step())
    this%hlast = this%hnext
  end subroutine


  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(t2d_flow_integrator), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%solver%get_cell_flow_soln(pressure, velocity)
  end subroutine


  subroutine get_cell_flow_active(this, active)
    class(t2d_flow_integrator), target, intent(in) :: this
    logical, pointer, intent(out) :: active(:)

    call this%solver%get_cell_flow_active(active)
  end subroutine


  subroutine get_face_velocity(this, velocity)
    class(t2d_flow_integrator), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%solver%get_face_velocity(velocity)
  end subroutine


  subroutine integrate(this, env, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(t2d_flow_integrator), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, hproposed
    logical :: sig_rcvd
    character(256) :: line
    character(8) :: cause

    stat = 0
    ASSERT(this%time_stepper_initialized)
    t_n = this%tlast
    ASSERT(tout >= t_n)
    do while (t_n < tout)
      hproposed = this%hnext
      call select_step_cause(this, cause)
      t_np1 = this%ts_sync%next_time(tout, t_n, this%hlast, this%hnext)
      if (t_np1 < t_n + hproposed) cause = 'output'
      if (t_np1 - t_n < this%dt_min) then
        stat = -1
        errmsg = 'next time step is too small'
        return
      end if
      write(line,'(a,i0,a,es0.5,a,es0.5,a,a)') 'step=', this%solver%num_steps() + 1_int64, &
          ' attempt=1 t0=', t_n, ' dt=', t_np1 - t_n, ' cause=', trim(cause)
      call env%simlog%begin_section(trim(line))
      call this%solver%step(env, t_n, t_np1, stat, errmsg)
      if (stat /= 0) then
        if (.not.allocated(errmsg)) errmsg = 'Navier--Stokes solver step failed'
        call env%simlog%end_section('step-end status=failed')
        return
      end if
      call env%simlog%end_section('step-end status=accepted')
      this%hlast = t_np1 - t_n
      this%hnext = min(this%dt_grow*this%hlast, this%dt_max, this%solver%courant_time_step())
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


  subroutine select_step_cause(this, cause)
    class(t2d_flow_integrator), intent(in) :: this
    character(*), intent(out) :: cause

    real(r8) :: h, hlimit

    if (this%solver%num_steps() == 0_int64) then
      h = this%dt_init
      cause = 'init'
      hlimit = this%solver%courant_time_step()
      if (hlimit < h) then
        h = hlimit
        cause = 'cfl'
      end if
    else
      h = this%dt_grow*this%hlast
      cause = 'growth'
      if (this%dt_max < h) then
        h = this%dt_max
        cause = 'max'
      end if
      hlimit = this%solver%courant_time_step()
      if (hlimit < h) cause = 'cfl'
    end if
  end subroutine


  function last_time(this) result(time)
    class(t2d_flow_integrator), intent(in) :: this
    real(r8) :: time

    time = this%tlast
  end function


  integer(int64) function num_steps(this)
    class(t2d_flow_integrator), intent(in) :: this

    num_steps = this%solver%num_steps()
  end function


  function initial_time_step(this) result(dt)
    class(t2d_flow_integrator), intent(in) :: this
    real(r8) :: dt

    ASSERT(this%time_stepper_initialized)
    dt = this%dt_init
  end function


  subroutine init_temporal_output(this, data)
    class(t2d_flow_integrator), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call this%solver%init_temporal_output(data)
  end subroutine


  subroutine set_temporal_output(this, data)
    class(t2d_flow_integrator), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call this%solver%set_temporal_output(data)
  end subroutine

end module t2d_flow_integrator_type
