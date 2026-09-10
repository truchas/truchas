!!
!! T2D_FLOW_THERMAL_SOLVER_TYPE
!!
!! This module defines T2D_FLOW_THERMAL_SOLVER, which coordinates one attempted
!! two-dimensional incompressible Navier--Stokes/thermal-transport step. It
!! advances material transport, converts its fluxes to an enthalpy rate,
!! attempts thermal transport, and then advances flow momentum and pressure.
!! The mesh, material distribution, and models remain sim-owned. The flow
!! solver owns the coupled step state and provides accessors for data needed
!! by the coupled time integrator.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flow_thermal_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use t2d_unstr_mesh_type
  use t2d_flow_model_type
  use t2d_flow_material_mapping_type
  use t2d_flow_material_transport_type
  use t2d_flow_mechanics_type
  use t2d_thermal_model_type
  use t2d_thermal_solver_type
  use t2d_flow_thermal_enthalpy_advector_type
  implicit none
  private

  type, public :: t2d_flow_thermal_solver
    private
    type(t2d_unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(material_distribution), pointer :: matl_dist => null() ! unowned reference
    type(t2d_flow_model), pointer :: flow_model => null() ! unowned reference
    type(t2d_flow_material_mapping) :: matl_map
    type(t2d_flow_mechanics) :: mechanics
    type(t2d_flow_material_transport) :: material_transport
    type(t2d_flow_thermal_enthalpy_advector) :: enthalpy_advector
    type(t2d_thermal_solver), pointer :: thermal => null()
    real(r8), allocatable :: temp(:), enthalpy_increment(:), flow_vfrac(:,:), flow_vfrac_old(:,:), &
        matl_vfrac_old(:,:)
    integer :: ncell_onP
    integer(int64) :: nstep = 0_int64
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: step
    procedure :: num_steps
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: get_cell_flow_active
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
    procedure :: courant_time_step
    final :: delete
  end type

contains

  subroutine init(this, env, flow_model, ht_model, matl_model, matl_dist, &
      params, stat, errmsg)

    class(t2d_flow_thermal_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_flow_model), target, intent(inout) :: flow_model
    type(t2d_thermal_model), target, intent(in) :: ht_model
    type(material_model), intent(in) :: matl_model
    type(material_distribution), target, intent(in) :: matl_dist
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: tracking_subcycles
    real(r8) :: courant_number, tracking_cutoff
    type(parameter_list), pointer :: flow_params, momentum_params, projection_params, thermal_params
    type(parameter_list), pointer :: tracking_params => null()
    character(:), allocatable :: tracking_algorithm
    integer, allocatable :: flow_pids(:), priority(:)
    logical :: simple_default

    stat = 0
    this%flow_model => flow_model
    ASSERT(size(matl_dist%vfrac,1) == matl_model%nmatl)
    if (.not.params%is_sublist('flow') .or. .not.params%is_sublist('thermal')) then
      stat = 1
      errmsg = 'solver requires flow and thermal sublists'
      return
    end if
    flow_params => params%sublist('flow')
    thermal_params => params%sublist('thermal')
    simple_default = .false.
    if (matl_model%nmatl_real == 1 .and. matl_model%nphase_real == 1 .and. .not.matl_model%have_void) &
      simple_default = matl_model%is_fluid(1)
    tracking_algorithm = 'geometric'
    tracking_cutoff = 1.0e-6_r8
    tracking_subcycles = 4
    if (simple_default) tracking_algorithm = 'simple'
    if (flow_params%is_sublist('volume-tracking')) then
      tracking_params => flow_params%sublist('volume-tracking')
      call tracking_params%get('algorithm', tracking_algorithm, default=tracking_algorithm, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call tracking_params%get('cutoff', tracking_cutoff, default=tracking_cutoff, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (tracking_cutoff <= 0.0_r8 .or. tracking_cutoff >= 1.0_r8) then
        stat = 1
        errmsg = 'solver.flow.volume-tracking.cutoff must be in (0,1)'
        return
      end if
      call tracking_params%get('subcycles', tracking_subcycles, default=tracking_subcycles, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (tracking_subcycles < 1) then
        stat = 1
        errmsg = 'solver.flow.volume-tracking.subcycles must be at least one'
        return
      end if
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
    call flow_params%get('courant-number', courant_number, default=0.5_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (courant_number <= 0.0_r8 .or. courant_number > 1.0_r8) then
      stat = 1
      errmsg = 'solver.flow.courant-number must be in (0,1]'
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
    allocate(this%temp(flow_model%mesh%ncell_onP), this%enthalpy_increment(flow_model%mesh%ncell_onP), &
        this%matl_vfrac_old(matl_model%nmatl,flow_model%mesh%ncell_onP), &
        this%flow_vfrac(this%matl_map%num_material(),flow_model%mesh%ncell), &
        this%flow_vfrac_old(this%matl_map%num_material(),flow_model%mesh%ncell))
    call this%matl_map%get_reduced_volume_fractions(matl_dist, this%flow_vfrac)
    call flow_model%mesh%cell_imap%gather_offp(this%flow_vfrac)
    this%flow_vfrac_old = this%flow_vfrac
    call flow_model%set_volume_fractions(this%flow_vfrac)
    if (flow_model%inviscid) then
      call this%mechanics%init(env, flow_model, projection_params=projection_params, courant_number=courant_number, &
          stat=stat, errmsg=errmsg)
    else
      if (.not.flow_params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a momentum-solver sublist'
        return
      end if
      momentum_params => flow_params%sublist('momentum-solver')
      call this%mechanics%init(env, flow_model, momentum_params, projection_params, courant_number, stat, errmsg)
    end if
    if (stat /= 0) return
    allocate(priority(this%matl_map%num_material()))
    call this%matl_map%get_priority(priority)
    call this%material_transport%init(env, flow_model%mesh, this%matl_map%num_real_fluid(), &
        this%matl_map%num_fluid(), this%matl_map%num_material(), algorithm=tracking_algorithm, &
        priority=priority, cutoff=tracking_cutoff, subcycles=tracking_subcycles)
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
    type(t2d_flow_thermal_solver), intent(inout) :: this

    if (associated(this%thermal)) deallocate(this%thermal)
  end subroutine


  !! Set the initial flow and thermal states at TIME using DT for their
  !! respective initial-condition procedures.
  subroutine set_initial_state(this, env, matl_model, time, dt, velocity, temp, stat, errmsg)

    class(t2d_flow_thermal_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: time, dt, velocity(:,:), temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%thermal%set_initial_state(env, time, dt, temp, stat, errmsg)
    if (stat /= 0) return
    call this%thermal%get_cell_temp_soln(this%temp)
    call this%matl_map%get_phase_volume_fractions(matl_model, this%matl_dist, this%temp, this%flow_vfrac)
    call this%mesh%cell_imap%gather_offp(this%flow_vfrac)
    call this%mechanics%set_initial_material_state(this%flow_vfrac, this%temp)
    call this%mechanics%set_buoyancy_temperature(this%temp)
    call this%mechanics%set_initial_state(env, time, dt, velocity, stat)
    if (stat /= 0) then
      errmsg = 'initializing flow state failed'
      return
    end if
    this%nstep = 0_int64
  end subroutine


  !! Attempt one coupled step from T_N to T_NP1. Thermal failure is reported
  !! as recoverable; a flow failure is non-recoverable.
  subroutine step(this, env, matl_model, t_n, t_np1, stat, errmsg, hnext)

    class(t2d_flow_thermal_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(out) :: hnext

    real(r8), pointer :: face_velocity(:)

    ASSERT(t_np1 > t_n)
    ASSERT(this%thermal%last_time() == t_n)
    stat = 0
    call env%timer%start('flow/material-transport')
    this%flow_vfrac_old = this%flow_vfrac
    call this%mechanics%get_face_velocity(face_velocity)
    call this%material_transport%advance(env, t_n, t_np1, face_velocity, this%flow_vfrac)
    call env%timer%stop('flow/material-transport')
    call this%thermal%get_cell_temp_soln(this%temp)
    this%matl_vfrac_old = this%matl_dist%vfrac
    call this%matl_map%apply_phase_fluxes(this%mesh, this%material_transport%flux_volumes, this%matl_dist)
    call this%matl_map%get_phase_volume_fractions(matl_model, this%matl_dist, this%temp, this%flow_vfrac)
    call this%mesh%cell_imap%gather_offp(this%flow_vfrac)
    call this%mechanics%set_volume_fractions(this%flow_vfrac)
    call this%mechanics%set_pre_solidification_state()
    call env%timer%start('thermal/advection')
    call this%enthalpy_advector%get_advected_enthalpy(t_n, this%temp, &
        this%material_transport%flux_volumes, this%enthalpy_increment)
    call env%timer%stop('thermal/advection')
    !! TODO: A future adaptive BDF1 thermal error estimate should measure
    !! the conduction update relative to this advected enthalpy state, so
    !! material-front motion does not masquerade as thermal truncation error.
    call this%thermal%set_ext_enthalpy_rate(this%enthalpy_increment / (t_np1 - t_n))
    call this%thermal%step(env, t_n, t_np1, stat, errmsg, hnext=hnext)
    if (stat /= 0) then
      this%matl_dist%vfrac = this%matl_vfrac_old
      this%flow_vfrac = this%flow_vfrac_old
      call this%mechanics%set_volume_fractions(this%flow_vfrac_old)
      call this%mechanics%set_pre_solidification_state()
      stat = 1
      return
    end if
    call this%thermal%get_cell_temp_soln(this%temp)
    call this%matl_map%get_phase_volume_fractions(matl_model, this%matl_dist, this%temp, this%flow_vfrac)
    call this%mesh%cell_imap%gather_offp(this%flow_vfrac)
    call this%mechanics%set_volume_fractions(this%flow_vfrac)
    call this%mechanics%set_buoyancy_temperature(this%temp)
    if (.not.this%flow_model%unsteady_stokes) then
      call this%mechanics%advance_momentum(env, t_n, t_np1, stat, errmsg, this%material_transport%flux_volumes)
    else
      call this%mechanics%advance_momentum(env, t_n, t_np1, stat, errmsg)
    end if
    if (stat /= 0) then
      call this%mechanics%reject_step()
      call this%thermal%reject_step()
      call this%thermal%get_cell_temp_soln(this%temp)
      this%matl_dist%vfrac = this%matl_vfrac_old
      this%flow_vfrac = this%flow_vfrac_old
      call this%mechanics%set_volume_fractions(this%flow_vfrac_old)
      call this%mechanics%set_pre_solidification_state()
      call this%mechanics%set_buoyancy_temperature(this%temp)
      stat = -3
      if (.not.allocated(errmsg)) errmsg = 'flow momentum update failed'
      return
    end if
    call this%mechanics%commit_step()
    call this%thermal%commit_step
    this%nstep = this%nstep + 1_int64
  end subroutine


  integer(int64) function num_steps(this)
    class(t2d_flow_thermal_solver), intent(in) :: this

    num_steps = this%nstep
  end function


  !! Declare the temporal scalar fields published by this coupled solver.
  !! These fields are updated at each requested solution output and written
  !! by the simulation's output writer.
  subroutine init_temporal_output(this, data)

    class(t2d_flow_thermal_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Set the current values of the temporal scalar fields published by this
  !! coupled solver.
  subroutine set_temporal_output(this, data)

    class(t2d_flow_thermal_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Returns the current local cell pressure and velocity, including ghosts.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(t2d_flow_thermal_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%mechanics%get_cell_flow_soln(pressure, velocity)
  end subroutine


  !! Returns the current face-normal velocity, including ghost faces.
  subroutine get_face_velocity(this, velocity)
    class(t2d_flow_thermal_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%mechanics%get_face_velocity(velocity)
  end subroutine


  !! Returns a no-copy view of the full-local mask used to distinguish genuine
  !! flow equations from dummy equations.
  subroutine get_cell_flow_active(this, active)
    class(t2d_flow_thermal_solver), target, intent(in) :: this
    logical, pointer, intent(out) :: active(:)

    call this%mechanics%get_cell_flow_active(active)
  end subroutine


  subroutine get_cell_heat_soln(this, enth)
    class(t2d_flow_thermal_solver), intent(in) :: this
    real(r8), intent(inout) :: enth(:)

    call this%thermal%get_cell_heat_soln(enth)
  end subroutine


  subroutine get_cell_temp_soln(this, temp)
    class(t2d_flow_thermal_solver), intent(in) :: this
    real(r8), intent(inout) :: temp(:)

    call this%thermal%get_cell_temp_soln(temp)
  end subroutine

  real(r8) function courant_time_step(this)
    class(t2d_flow_thermal_solver), intent(in) :: this

    courant_time_step = this%mechanics%courant_time_step()
  end function

end module t2d_flow_thermal_solver_type
