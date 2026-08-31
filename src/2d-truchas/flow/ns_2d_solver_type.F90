!!
!! NS_2D_SOLVER_TYPE
!!
!! This module defines NS_2D_SOLVER, the isothermal incompressible
!! Navier--Stokes orchestration layer.  It owns material transport, the
!! standalone time-step policy, and the count of successful steps.  It
!! delegates flow mechanics and flow-state management to FLOW_2D_SOLVER.
!! Momentum transport is an explicit first-order donor-cell contribution
!! supplied to that common solver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ns_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use flow_2d_model_type
  use flow_2d_solver_type
  use flow_2d_material_layout_type
  use flow_2d_material_transport_type
  use time_step_sync_type
  implicit none
  private

  type, public :: ns_2d_solver
    private
    type(flow_2d_solver) :: flow
    type(flow_2d_material_transport) :: material_transport
    real(r8), allocatable :: vfrac(:,:)
    integer(int64) :: nstep = 0_int64
    real(r8) :: tlast, hlast, hnext
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    logical :: time_stepper_initialized = .false.
    type(time_step_sync) :: ts_sync
  contains
    procedure :: init
    procedure :: init_time_stepper
    procedure :: set_volume_fractions
    procedure :: set_initial_material_state
    procedure :: get_volume_fractions
    procedure :: set_buoyancy_temperature
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: step
    procedure :: integrate
    procedure :: last_time
    procedure :: num_steps
    procedure :: initial_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: advance_momentum
    procedure :: commit_step
    procedure :: reject_step
    procedure :: courant_time_step
    final :: delete
  end type

contains

  subroutine init(this, env, model, momentum_params, projection_params, material_layout, tracking_params, courant_number)
    class(ns_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(in) :: model
    type(parameter_list), target, intent(in), optional :: momentum_params
    type(parameter_list), target, intent(in) :: projection_params
    type(flow_2d_material_layout), intent(in), optional :: material_layout
    type(parameter_list), target, intent(inout), optional :: tracking_params
    real(r8), intent(in), optional :: courant_number

    integer :: nrealfluid, nfluid, nmat
    integer, allocatable :: priority(:)
    character(:), allocatable :: algorithm

    if (present(courant_number)) then
      call this%flow%init(env, model, momentum_params, projection_params, courant_number)
    else
      call this%flow%init(env, model, momentum_params, projection_params)
    end if
    if (present(material_layout)) then
      nrealfluid = material_layout%num_real_fluid()
      nfluid = material_layout%num_fluid()
      nmat = material_layout%num_material()
      allocate(priority(nmat))
      call material_layout%get_priority(priority)
    else
      nrealfluid = 1
      nfluid = 1
      nmat = 1
    end if
    if (present(tracking_params)) then
      call tracking_params%get('algorithm', algorithm, default='simple')
    else
      algorithm = 'simple'
    end if
    allocate(this%vfrac(nmat,model%mesh%ncell))
    this%vfrac = 0.0_r8
    this%vfrac(1,:) = 1.0_r8
    if (present(material_layout)) then
      call this%material_transport%init(env, model%mesh, nrealfluid, nfluid, nmat, algorithm, priority)
    else
      call this%material_transport%init(env, model%mesh, 1, 1, 1)
    end if
  end subroutine


  !! Initialize the standalone time-step policy from SIM-CONTROL parameters.
  !! The output schedule itself remains owned by the simulation driver.
  subroutine init_time_stepper(this, params, stat, errmsg)
    class(ns_2d_solver), intent(inout) :: this
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
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    call this%flow%set_volume_fractions(vfrac)
  end subroutine


  subroutine set_initial_material_state(this, vfrac, temperature)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%flow%set_initial_material_state(vfrac, temperature)
    ASSERT(size(vfrac,1) == size(this%vfrac,1))
    ASSERT(size(vfrac,2) == size(this%vfrac,2))
    this%vfrac = vfrac
  end subroutine


  subroutine get_volume_fractions(this, vfrac)
    class(ns_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: vfrac(:,:)

    vfrac => this%vfrac
  end subroutine


  subroutine set_buoyancy_temperature(this, temperature)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    call this%flow%set_buoyancy_temperature(temperature)
  end subroutine


  subroutine delete(this)
    type(ns_2d_solver), intent(inout) :: this
  end subroutine


  !! Set STATE from an input velocity. The common initial-condition solver
  !! projects the velocity and computes an initial pressure with its temporary
  !! Stokes step, as mainline does when it omits initial momentum transport.
  subroutine set_initial_state(this, env, time, dt, velocity, stat)
    class(ns_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    call this%flow%set_initial_state(env, time, dt, velocity, stat)
    if (stat /= 0) return
    this%nstep = 0_int64
    if (this%time_stepper_initialized) then
      this%tlast = time
      this%hnext = min(this%dt_init, this%courant_time_step())
      this%hlast = this%hnext
    end if
  end subroutine


  !! Return no-copy views of the current cell-centered pressure and velocity.
  !! The current state is pending after a successful ADVANCE_MOMENTUM call;
  !! callers must reacquire these views after COMMIT_STEP or REJECT_STEP.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(ns_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%flow%get_cell_flow_soln(pressure, velocity)
  end subroutine


  !! Return a no-copy view of the current face-normal velocity.  The current
  !! state is pending after a successful ADVANCE_MOMENTUM call; callers must
  !! reacquire this view after COMMIT_STEP or REJECT_STEP.
  subroutine get_face_velocity(this, velocity)
    class(ns_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%flow%get_face_velocity(velocity)
  end subroutine


  !! Advance STATE from T_N to T_NP1. The time step is derived from the two
  !! endpoint times so callers retain exact target times. This is the
  !! isothermal wrapper: it first obtains material transport from the old face
  !! velocity and then advances momentum and pressure.
  subroutine step(this, env, t_n, t_np1, stat, errmsg, step_cause)
    class(ns_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg
    character(*), optional, intent(in) :: step_cause

    real(r8), pointer :: vfrac_trial(:,:), face_velocity(:)
    character(256) :: line
    character(8) :: cause

    cause = 'explicit'
    if (present(step_cause)) cause = step_cause
    write(line,'(a,i0,a,es0.5,a,es0.5,a,a)') 'step=', this%nstep + 1_int64, &
        ' attempt=1 t0=', t_n, ' dt=', t_np1 - t_n, ' cause=', trim(cause)
    call env%simlog%begin_section(trim(line))
    call env%timer%start('flow/material-transport')
    call this%flow%get_face_velocity(face_velocity)
    call this%material_transport%advance(env, t_n, t_np1, face_velocity, this%vfrac)
    call this%material_transport%get_trial_volume_fractions(vfrac_trial)
    call env%timer%stop('flow/material-transport')
    call this%flow%set_volume_fractions(vfrac_trial)
    call this%flow%advance_momentum(env, t_n, t_np1, stat, errmsg, this%material_transport%flux_volumes)
    if (stat /= 0) then
      call this%flow%reject_step()
      call this%flow%set_volume_fractions(this%vfrac)
      call env%simlog%end_section('step-end status=failed')
      return
    end if
    this%vfrac = vfrac_trial
    call this%commit_step()
    call env%simlog%end_section('step-end status=accepted')
  end subroutine


  !! Integrate the flow state from its current time to TOUT.  The endpoint
  !! time is primary; each step derives its size from the two endpoint times.
  subroutine integrate(this, env, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(ns_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, hproposed
    logical :: sig_rcvd
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
      call this%step(env, t_n, t_np1, stat, errmsg, step_cause=cause)
      if (stat /= 0) then
        if (.not.allocated(errmsg)) errmsg = 'Navier--Stokes solver step failed'
        return
      end if
      this%hlast = t_np1 - t_n
      this%hnext = min(this%dt_grow*this%hlast, this%dt_max, &
          this%courant_time_step())
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
    class(ns_2d_solver), intent(in) :: this
    character(*), intent(out) :: cause

    real(r8) :: h, hlimit

    if (this%nstep == 0_int64) then
      h = this%dt_init
      cause = 'init'
      hlimit = this%courant_time_step()
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
      hlimit = this%courant_time_step()
      if (hlimit < h) cause = 'cfl'
    end if
  end subroutine


  function last_time(this) result(time)
    class(ns_2d_solver), intent(in) :: this
    real(r8) :: time

    time = this%tlast
  end function


  integer(int64) function num_steps(this)
    class(ns_2d_solver), intent(in) :: this

    num_steps = this%nstep
  end function


  function initial_time_step(this) result(dt)
    class(ns_2d_solver), intent(in) :: this
    real(r8) :: dt

    ASSERT(this%time_stepper_initialized)
    dt = this%dt_init
  end function


  !! Advance momentum and pressure from T_N to T_NP1 using material-resolved
  !! flux volumes already constructed for the pending step. This separate
  !! operation lets a coupled solver advect material and thermal enthalpy
  !! before it updates the flow state.
  subroutine advance_momentum(this, env, t_n, t_np1, flux_volumes, stat, errmsg)
    class(ns_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    real(r8), intent(in) :: flux_volumes(:,:)
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg

    call this%flow%advance_momentum(env, t_n, t_np1, stat, errmsg, flux_volumes)
  end subroutine

  !! Commit the pending flow state and its current material properties.
  subroutine commit_step(this)
    class(ns_2d_solver), intent(inout) :: this

    call this%flow%commit_step()
    this%nstep = this%nstep + 1_int64
  end subroutine


  !! Declare the temporal scalar fields published by this solver.
  !! These fields are updated at each requested solution output and written
  !! by the simulation's output writer.
  subroutine init_temporal_output(this, data)
    class(ns_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Set the current values of the temporal scalar fields published by this
  !! solver.
  subroutine set_temporal_output(this, data)
    class(ns_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Reject the pending flow state, restoring the last accepted state.
  !! Material volume fractions are restored by the owning caller.
  subroutine reject_step(this)
    class(ns_2d_solver), intent(inout) :: this

    call this%flow%reject_step()
  end subroutine

  !! Return the maximum step size requested by the flow mechanics for the old
  !! face-normal velocity.
  function courant_time_step(this) result(dt)
    class(ns_2d_solver), intent(in) :: this
    real(r8) :: dt

    dt = this%flow%courant_time_step()
  end function

end module ns_2d_solver_type
