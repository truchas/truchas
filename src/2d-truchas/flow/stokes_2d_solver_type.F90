!!
!! STOKES_2D_SOLVER_TYPE
!!
!! This module defines STOKES_2D_SOLVER, the two-dimensional incompressible
!! Stokes orchestration layer.  It owns the material mapping, standalone
!! time-step policy, integration loop, and count of successful steps.  It
!! delegates the momentum and pressure-correction work, flow state, and
!! initial-condition solve to FLOW_2D_SOLVER.  The latter is also used by the
!! Navier--Stokes solver, so Stokes and Navier--Stokes share the same flow
!! mechanics and material-property path.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module stokes_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use unstr_2d_mesh_type
  use flow_2d_model_type
  use flow_2d_solver_type
  use flow_material_mapping_type
  use time_step_sync_type
  implicit none
  private

  type, public :: stokes_2d_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(flow_material_mapping) :: matl_map
    type(flow_2d_solver) :: flow
    integer(int64) :: nstep = 0_int64
    real(r8) :: tlast, hlast, hnext
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    logical :: time_stepper_initialized = .false.
    type(time_step_sync) :: ts_sync
  contains
    procedure :: init
    procedure :: init_time_stepper
    procedure :: set_initial_material_state
    procedure :: get_reduced_volume_fractions
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_cell_flow_active
    procedure :: get_face_velocity
    procedure :: step
    procedure :: integrate
    procedure :: last_time
    procedure :: num_steps
    procedure :: initial_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: reject_step
    procedure :: courant_time_step
  end type

contains

  subroutine init(this, env, model, matl_model, params, stat, errmsg)
    class(stokes_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(inout) :: model
    type(material_model), intent(in) :: matl_model
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: nrealfluid, nfluid
    integer, allocatable :: phase_ids(:)
    type(parameter_list), pointer :: momentum_params, projection_params
    real(r8) :: courant_number

    stat = 0
    this%mesh => model%mesh
    if (matl_model%nphase_real /= matl_model%nmatl_real) then
      stat = 1
      errmsg = 'Stokes flow requires single-phase materials'
      return
    end if
    call this%matl_map%init(matl_model, stat, errmsg)
    if (stat /= 0) return

    call params%get('courant-number', courant_number, default=0.5_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // params%path() // ': ' // errmsg
      return
    end if
    if (courant_number <= 0.0_r8 .or. courant_number > 1.0_r8) then
      stat = 1
      errmsg = 'processing ' // params%path() // ': "courant-number" must be in (0,1]'
      return
    end if

    if (.not.params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'requires a "projection-solver" sublist'
      return
    end if
    projection_params => params%sublist('projection-solver')
    if (.not.model%inviscid) then
      if (.not.params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a "momentum-solver" sublist'
        return
      end if
      momentum_params => params%sublist('momentum-solver')
    end if

    nrealfluid = this%matl_map%num_real_fluid()
    nfluid = this%matl_map%num_fluid()
    allocate(phase_ids(nrealfluid))
    call this%matl_map%get_real_fluid_phase_ids(phase_ids)
    call model%init_material(matl_model, phase_ids, stat, errmsg, nfluid=nfluid)
    if (stat /= 0) return

    if (model%inviscid) then
      call this%flow%init(env, model, projection_params=projection_params, courant_number=courant_number, &
          stat=stat, errmsg=errmsg)
    else
      call this%flow%init(env, model, momentum_params, projection_params, courant_number, stat, errmsg)
    end if
  end subroutine


  !! Initialize the standalone time-step policy from SIM-CONTROL parameters.
  !! The output schedule itself remains owned by the simulation driver.
  subroutine init_time_stepper(this, params, stat, errmsg)
    class(stokes_2d_solver), intent(inout) :: this
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


  subroutine set_initial_material_state(this, vfrac, temperature)
    class(stokes_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%flow%set_initial_material_state(vfrac, temperature)
  end subroutine


  subroutine get_reduced_volume_fractions(this, matl_dist, vfrac)
    class(stokes_2d_solver), intent(in) :: this
    type(material_distribution), intent(in) :: matl_dist
    real(r8), allocatable, intent(out) :: vfrac(:,:)

    allocate(vfrac(this%matl_map%num_material(), this%mesh%ncell))
    call this%matl_map%get_reduced_volume_fractions(matl_dist, vfrac)
  end subroutine


  subroutine set_initial_state(this, env, time, dt, velocity, stat)
    class(stokes_2d_solver), intent(inout) :: this
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


  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(stokes_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%flow%get_cell_flow_soln(pressure, velocity)
  end subroutine


  subroutine get_cell_flow_active(this, active)
    class(stokes_2d_solver), target, intent(in) :: this
    logical, pointer, intent(out) :: active(:)

    call this%flow%get_cell_flow_active(active)
  end subroutine


  subroutine get_face_velocity(this, velocity)
    class(stokes_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%flow%get_face_velocity(velocity)
  end subroutine


  subroutine step(this, env, t_n, t_np1, stat, errmsg, step_cause)
    class(stokes_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg
    character(*), optional, intent(in) :: step_cause

    call this%flow%step(env, t_n, t_np1, stat, errmsg, step_cause)
    if (stat == 0) this%nstep = this%nstep + 1_int64
  end subroutine


  subroutine integrate(this, env, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(stokes_2d_solver), intent(inout) :: this
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
        if (.not.allocated(errmsg)) errmsg = 'Stokes solver step failed'
        return
      end if
      this%hlast = t_np1 - t_n
      this%hnext = min(this%dt_grow*this%hlast, this%dt_max, this%courant_time_step())
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
    class(stokes_2d_solver), intent(in) :: this
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
    class(stokes_2d_solver), intent(in) :: this
    real(r8) :: time

    time = this%tlast
  end function


  integer(int64) function num_steps(this)
    class(stokes_2d_solver), intent(in) :: this

    num_steps = this%nstep
  end function


  function initial_time_step(this) result(dt)
    class(stokes_2d_solver), intent(in) :: this
    real(r8) :: dt

    ASSERT(this%time_stepper_initialized)
    dt = this%dt_init
  end function


  subroutine init_temporal_output(this, data)
    class(stokes_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  subroutine set_temporal_output(this, data)
    class(stokes_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  subroutine reject_step(this)
    class(stokes_2d_solver), intent(inout) :: this

    call this%flow%reject_step()
  end subroutine


  function courant_time_step(this) result(dt)
    class(stokes_2d_solver), intent(in) :: this
    real(r8) :: dt

    dt = this%flow%courant_time_step()
  end function

end module stokes_2d_solver_type
