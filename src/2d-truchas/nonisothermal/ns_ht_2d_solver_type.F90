!!
!! NS_HT_2D_SOLVER_TYPE
!!
!! This module defines NS_HT_2D_SOLVER, which coordinates an attempted
!! two-dimensional incompressible Navier--Stokes/thermal-transport step. It
!! advances material transport, converts its fluxes to an enthalpy rate,
!! attempts thermal transport, and then advances flow momentum and pressure.
!! The mesh, material distribution, and models remain sim-owned. The flow
!! solver owns the flow state and provides accessors for data needed by the
!! coupled simulation driver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ns_ht_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use unstr_2d_mesh_type
  use flow_2d_model_type
  use flow_material_mapping_type
  use flow_2d_material_transport_type
  use flow_2d_solver_type
  use ht_2d_model_type
  use ht_2d_solver_type
  use ns_ht_2d_enthalpy_advector_type
  use time_step_sync_type
  implicit none
  private

  type, public :: ns_ht_2d_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    type(material_distribution), pointer :: matl_dist => null() ! unowned reference
    type(flow_material_mapping) :: matl_map
    type(flow_2d_solver) :: flow
    type(flow_2d_material_transport) :: material_transport
    type(ns_ht_2d_enthalpy_advector) :: enthalpy_advector
    type(ht_2d_solver), pointer :: thermal => null()
    real(r8), allocatable :: temp(:), enthalpy_increment(:), flow_vfrac(:,:), flow_vfrac_old(:,:), &
        matl_vfrac_old(:,:)
    integer :: ncell_onP
    logical :: inertial = .true.
    integer(int64) :: nstep = 0_int64
    real(r8) :: dt_init, dt_min, dt_max, dt_grow, hnext, hlast
    integer :: max_try
    character(8) :: hnext_cause = 'init'
    type(time_step_sync) :: ts_sync
  contains
    procedure :: init
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
    final :: delete
  end type

contains

  subroutine init(this, env, flow_model, ht_model, matl_model, matl_dist, &
      params, stat, errmsg)

    class(ns_ht_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(inout) :: flow_model
    type(ht_2d_model), target, intent(in) :: ht_model
    type(material_model), intent(in) :: matl_model
    type(material_distribution), target, intent(in) :: matl_dist
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: lookahead
    real(r8) :: courant_number
    type(parameter_list), pointer :: flow_params, momentum_params, projection_params, thermal_params
    type(parameter_list), pointer :: tracking_params => null()
    character(:), allocatable :: tracking_algorithm
    integer, allocatable :: flow_pids(:), priority(:)
    logical :: simple_default

    stat = 0
    ASSERT(size(matl_dist%vfrac,1) == matl_model%nmatl)
    if (.not.params%is_sublist('flow') .or. .not.params%is_sublist('thermal')) then
      stat = 1
      errmsg = 'solver requires flow and thermal sublists'
      return
    end if
    flow_params => params%sublist('flow')
    thermal_params => params%sublist('thermal')
    call flow_params%get('inertial', this%inertial, default=.true., stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (flow_model%inviscid .and. .not.this%inertial) then
      stat = 1
      errmsg = 'inviscid flow is incompatible with non-inertial flow'
      return
    end if
    if (this%inertial) then
      call env%simlog%info('Using inertial momentum.')
    else
      call env%simlog%info('Using non-inertial momentum (Stokes).')
    end if
    simple_default = .false.
    if (matl_model%nmatl_real == 1 .and. matl_model%nphase_real == 1 .and. .not.matl_model%have_void) &
      simple_default = matl_model%is_fluid(1)
    tracking_algorithm = 'geometric'
    if (simple_default) tracking_algorithm = 'simple'
    if (flow_params%is_sublist('volume-tracking')) then
      tracking_params => flow_params%sublist('volume-tracking')
      call tracking_params%get('algorithm', tracking_algorithm, default=tracking_algorithm, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
    end if
    if (tracking_algorithm /= 'simple' .and. tracking_algorithm /= 'geometric') then
      stat = 1
      errmsg = 'solver.flow.volume-tracking.algorithm must be "simple" or "geometric"'
      return
    end if
    call this%matl_map%init(matl_model, stat, errmsg)
    if (stat /= 0) return
    if (associated(tracking_params)) then
      call this%matl_map%set_priority(tracking_params, stat, errmsg)
      if (stat /= 0) return
    end if
    if (this%matl_map%num_real_fluid() == 0) then
      stat = 1
      errmsg = 'non-isothermal flow requires at least one fluid material'
      return
    end if
    allocate(flow_pids(this%matl_map%num_real_fluid()))
    call this%matl_map%get_real_fluid_phase_ids(flow_pids)
    call flow_model%init_material(matl_model, flow_pids, stat, errmsg, boussinesq=.true., &
        nfluid=this%matl_map%num_fluid())
    if (stat /= 0) return
    call params%get('initial-time-step', this%dt_init, stat, errmsg)
    if (stat /= 0) return
    call params%get('min-time-step', this%dt_min, stat, errmsg)
    if (stat /= 0) return
    call params%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-growth', this%dt_grow, default=1.05_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call flow_params%get('courant-number', courant_number, default=0.5_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-lookahead', lookahead, default=3, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. courant_number <= 0.0_r8 .or. &
        courant_number > 1.0_r8 .or. lookahead < 1) then
      stat = 1
      errmsg = 'invalid coupled time-step controls, including time-step-lookahead >= 1'
      return
    end if
    if (this%max_try <= 0) then
      stat = 1
      errmsg = 'maximum coupled step attempts must be > 0'
      return
    end if
    if (.not.flow_params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'solver.flow requires a projection-solver sublist'
      return
    end if
    projection_params => flow_params%sublist('projection-solver')

    this%ncell_onP = flow_model%mesh%ncell_onP
    this%mesh => flow_model%mesh
    this%matl_dist => matl_dist
    this%ts_sync = time_step_sync(lookahead)
    allocate(this%temp(flow_model%mesh%ncell_onP), this%enthalpy_increment(flow_model%mesh%ncell_onP), &
        this%matl_vfrac_old(matl_model%nmatl,flow_model%mesh%ncell_onP), &
        this%flow_vfrac(this%matl_map%num_material(),flow_model%mesh%ncell), &
        this%flow_vfrac_old(this%matl_map%num_material(),flow_model%mesh%ncell))
    call this%matl_map%get_reduced_volume_fractions(matl_dist, this%flow_vfrac)
    call flow_model%mesh%cell_imap%gather_offp(this%flow_vfrac)
    this%flow_vfrac_old = this%flow_vfrac
    call flow_model%set_volume_fractions(this%flow_vfrac)
    if (flow_model%inviscid) then
      call this%flow%init(env, flow_model, projection_params=projection_params, courant_number=courant_number, &
          stat=stat, errmsg=errmsg)
    else
      if (.not.flow_params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a momentum-solver sublist'
        return
      end if
      momentum_params => flow_params%sublist('momentum-solver')
      call this%flow%init(env, flow_model, momentum_params, projection_params, courant_number, stat, errmsg)
    end if
    if (stat /= 0) return
    allocate(priority(this%matl_map%num_material()))
    call this%matl_map%get_priority(priority)
    call this%material_transport%init(env, flow_model%mesh, this%matl_map%num_real_fluid(), &
        this%matl_map%num_fluid(), this%matl_map%num_material(), algorithm=tracking_algorithm, &
        priority=priority)
    if (allocated(ht_model%bc_inflow)) then
      call this%enthalpy_advector%init(flow_model%mesh, matl_model, flow_pids, stat, errmsg, &
          inflow_temperature=ht_model%bc_inflow)
    else
      call this%enthalpy_advector%init(flow_model%mesh, matl_model, flow_pids, stat, errmsg)
    end if
    if (stat /= 0) return
    allocate(this%thermal)
    call this%thermal%init(env, ht_model, thermal_params, stat, errmsg)
  end subroutine


  subroutine delete(this)
    type(ns_ht_2d_solver), intent(inout) :: this

    if (associated(this%thermal)) deallocate(this%thermal)
  end subroutine


  !! Set the initial flow and thermal states at TIME. The configured initial
  !! time step is supplied to their respective initial-condition procedures.
  subroutine set_initial_state(this, env, matl_model, time, velocity, temp, stat, errmsg)

    class(ns_ht_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: time, velocity(:,:), temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8) :: hlimit

    call this%thermal%set_initial_state(env, time, temp, stat, errmsg, dt=this%dt_init)
    if (stat /= 0) return
    call this%thermal%get_cell_temp_soln(this%temp)
    call this%matl_map%get_phase_volume_fractions(matl_model, this%matl_dist, this%temp, this%flow_vfrac)
    call this%mesh%cell_imap%gather_offp(this%flow_vfrac)
    call this%flow%set_initial_material_state(this%flow_vfrac, this%temp)
    call this%flow%set_buoyancy_temperature(this%temp)
    call this%flow%set_initial_state(env, time, this%dt_init, velocity, stat)
    if (stat /= 0) then
      errmsg = 'initializing flow state failed'
      return
    end if
    this%nstep = 0_int64
    this%hnext = this%dt_init
    this%hnext_cause = 'init'
    hlimit = this%flow%courant_time_step()
    if (hlimit < this%hnext) then
      this%hnext = hlimit
      this%hnext_cause = 'cfl'
    end if
    this%hlast = this%hnext
  end subroutine


  !! Integrate the coupled state from its current time to TOUT. This owns
  !! output-time synchronization, thermal retry, time-step growth, and the
  !! convective Courant restriction.
  subroutine integrate(this, env, matl_model, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(ns_ht_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, t, thermal_hnext, hproposed
    logical :: sig_rcvd
    character(8) :: cause

    stat = 0
    t_n = this%thermal%last_time()
    ASSERT(tout >= t_n)
    do while (t_n < tout)
      hproposed = this%hnext
      cause = this%hnext_cause
      t_np1 = this%ts_sync%next_time(tout, t_n, this%hlast, this%hnext)
      if (t_np1 < t_n + hproposed) cause = 'output'
      call attempt_step(this, env, matl_model, t_n, t_np1, cause, t, thermal_hnext, stat, errmsg)
      if (stat /= 0) return
      this%hlast = t - t_n
      call select_step_cause(this, thermal_hnext, this%hnext, this%hnext_cause)
      t_n = t
      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if
    end do
  end subroutine


  !! Attempt a coupled step from T_N toward T_NP1. T is the accepted endpoint:
  !! it normally equals T_NP1, but is smaller after thermal recovery. A flow
  !! failure is non-recoverable.
  subroutine attempt_step(this, env, matl_model, t_n, t_np1, step_cause, t, hnext, stat, errmsg)

    class(ns_ht_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: t_n, t_np1
    character(*), intent(in) :: step_cause
    real(r8), intent(out) :: t, hnext
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: t_try
    real(r8), pointer :: face_velocity(:)
    character(256) :: line
    character(8) :: attempt_cause

    ASSERT(t_np1 > t_n)
    ASSERT(this%thermal%last_time() == t_n)
    t = t_n
    t_try = t_np1
    do n = 1, this%max_try
      attempt_cause = step_cause
      if (n > 1) attempt_cause = 'thermal'
      write(line,'(a,i0,a,i0,a,es0.5,a,es0.5,a,a)') 'step=', this%nstep + 1_int64, &
          ' attempt=', n, ' t0=', t_n, ' dt=', t_try - t_n, ' cause=', trim(attempt_cause)
      call env%simlog%begin_section(trim(line))
      call env%timer%start('flow/material-transport')
      this%flow_vfrac_old = this%flow_vfrac
      call this%flow%get_face_velocity(face_velocity)
      call this%material_transport%advance(env, t_n, t_try, face_velocity, this%flow_vfrac)
      call env%timer%stop('flow/material-transport')
      call this%thermal%get_cell_temp_soln(this%temp)
      this%matl_vfrac_old = this%matl_dist%vfrac
      call this%matl_map%apply_phase_fluxes(this%mesh, this%material_transport%flux_volumes, this%matl_dist)
      call this%matl_map%get_phase_volume_fractions(matl_model, this%matl_dist, this%temp, this%flow_vfrac)
      call this%mesh%cell_imap%gather_offp(this%flow_vfrac)
      call this%flow%set_volume_fractions(this%flow_vfrac)
      call this%flow%set_pre_solidification_state()
      call env%timer%start('thermal/advection')
      call this%enthalpy_advector%get_advected_enthalpy(t_n, this%temp, &
          this%material_transport%flux_volumes, this%enthalpy_increment)
      call env%timer%stop('thermal/advection')
      !! TODO: A future adaptive BDF1 thermal error estimate should measure
      !! the conduction update relative to this advected enthalpy state, so
      !! material-front motion does not masquerade as thermal truncation error.
      call this%thermal%set_ext_enthalpy_rate(this%enthalpy_increment / (t_try - t_n))
      call this%thermal%step(env, t_n, t_try, stat, errmsg, hnext=hnext)
      if (stat /= 0) then
        this%matl_dist%vfrac = this%matl_vfrac_old
        this%flow_vfrac = this%flow_vfrac_old
        call this%flow%set_volume_fractions(this%flow_vfrac_old)
        call this%flow%set_pre_solidification_state()
        t_try = t_n + hnext
        if (t_try - t_n < this%dt_min) then
          stat = -1
          errmsg = 'next coupled time step is too small'
          call env%simlog%end_section('step-end status=failed')
          return
        end if
        call env%simlog%end_section('step-end status=rejected')
        cycle
      end if
      call this%thermal%get_cell_temp_soln(this%temp)
      call this%matl_map%get_phase_volume_fractions(matl_model, this%matl_dist, this%temp, this%flow_vfrac)
      call this%mesh%cell_imap%gather_offp(this%flow_vfrac)
      call this%flow%set_volume_fractions(this%flow_vfrac)
      call this%flow%set_buoyancy_temperature(this%temp)
      if (this%inertial) then
        call this%flow%advance_momentum(env, t_n, t_try, stat, errmsg, this%material_transport%flux_volumes)
      else
        call this%flow%advance_momentum(env, t_n, t_try, stat, errmsg)
      end if
      if (stat /= 0) then
        call this%flow%reject_step()
        call this%thermal%reject_step()
        call this%thermal%get_cell_temp_soln(this%temp)
        this%matl_dist%vfrac = this%matl_vfrac_old
        this%flow_vfrac = this%flow_vfrac_old
        call this%flow%set_volume_fractions(this%flow_vfrac_old)
        call this%flow%set_pre_solidification_state()
        call this%flow%set_buoyancy_temperature(this%temp)
        stat = -3
        if (.not.allocated(errmsg)) errmsg = 'flow momentum update failed'
        call env%simlog%end_section('step-end status=failed')
        return
      end if
      call this%flow%commit_step()
      call this%thermal%commit_step
      this%nstep = this%nstep + 1_int64
      t = t_try
      call env%simlog%end_section('step-end status=accepted')
      if (stat == 0) exit
    end do

    if (stat == 0) return

    stat = -2
    errmsg = 'unable to take a coupled time step'
  end subroutine


  subroutine select_step_cause(this, thermal_hnext, hnext, cause)
    class(ns_ht_2d_solver), intent(in) :: this
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
    hlimit = this%flow%courant_time_step()
    if (hlimit < hnext) then
      hnext = hlimit
      cause = 'cfl'
    end if
  end subroutine


  real(r8) function last_time(this)
    class(ns_ht_2d_solver), intent(in) :: this

    last_time = this%thermal%last_time()
  end function


  real(r8) function initial_time_step(this)
    class(ns_ht_2d_solver), intent(in) :: this

    initial_time_step = this%dt_init
  end function


  integer(int64) function num_steps(this)
    class(ns_ht_2d_solver), intent(in) :: this

    num_steps = this%nstep
  end function


  !! Declare the temporal scalar fields published by this coupled solver.
  !! These fields are updated at each requested solution output and written
  !! by the simulation's output writer.
  subroutine init_temporal_output(this, data)

    class(ns_ht_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Set the current values of the temporal scalar fields published by this
  !! coupled solver.
  subroutine set_temporal_output(this, data)

    class(ns_ht_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Returns the current local cell pressure and velocity, including ghosts.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(ns_ht_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%flow%get_cell_flow_soln(pressure, velocity)
  end subroutine


  !! Returns the current face-normal velocity, including ghost faces.
  subroutine get_face_velocity(this, velocity)
    class(ns_ht_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%flow%get_face_velocity(velocity)
  end subroutine


  !! Returns a no-copy view of the full-local mask used to distinguish genuine
  !! flow equations from dummy equations.
  subroutine get_cell_flow_active(this, active)
    class(ns_ht_2d_solver), target, intent(in) :: this
    logical, pointer, intent(out) :: active(:)

    call this%flow%get_cell_flow_active(active)
  end subroutine


  subroutine get_cell_heat_soln(this, enth)
    class(ns_ht_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: enth(:)

    call this%thermal%get_cell_heat_soln(enth)
  end subroutine


  subroutine get_cell_temp_soln(this, temp)
    class(ns_ht_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: temp(:)

    call this%thermal%get_cell_temp_soln(temp)
  end subroutine

end module ns_ht_2d_solver_type
