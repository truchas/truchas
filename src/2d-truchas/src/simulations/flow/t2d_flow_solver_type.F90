!!
!! T2D_FLOW_SOLVER_TYPE
!!
!! This module defines T2D_FLOW_SOLVER, the isothermal incompressible
!! Navier--Stokes step solver. It advances material transport and flow
!! mechanics as one transaction, retaining the material-resolved donor-cell
!! fluxes needed by the advective momentum update. Time-step selection and
!! target-time integration are provided by T2D_FLOW_INTEGRATOR.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flow_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use t2d_flow_model_type
  use t2d_flow_bc_type
  use t2d_flow_mechanics_type
  use t2d_flow_material_mapping_type
  use t2d_flow_material_transport_type
  implicit none
  private

  type, public :: t2d_flow_solver
    private
    type(t2d_flow_material_mapping) :: matl_map
    type(t2d_flow_model), pointer :: model => null()  ! unowned reference
    type(t2d_flow_mechanics) :: mechanics
    type(t2d_flow_material_transport) :: material_transport
    real(r8), allocatable :: vfrac(:,:)
    integer(int64) :: nstep = 0_int64
  contains
    procedure :: init
    procedure :: set_volume_fractions
    procedure :: set_initial_material_state
    procedure :: get_reduced_volume_fractions
    procedure :: update_material_distribution
    procedure :: set_buoyancy_temperature
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_cell_flow_active
    procedure :: get_face_velocity
    procedure :: step
    procedure :: num_steps
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: courant_time_step
  end type

contains

  subroutine init(this, env, model, matl_model, params, stat, errmsg)
    class(t2d_flow_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_flow_model), target, intent(inout) :: model
    type(material_model), intent(in) :: matl_model
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: nrealfluid, nfluid, nmat
    integer, allocatable :: priority(:), phase_ids(:)
    character(:), allocatable :: algorithm
    type(parameter_list), pointer :: tracking_params => null(), momentum_params, projection_params
    real(r8) :: courant_number, tracking_cutoff
    integer :: tracking_subcycles
    character(96) :: message
    logical :: simple_default

    stat = 0
    this%model => model
    if (matl_model%nphase_real /= matl_model%nmatl_real) then
      stat = 1
      errmsg = 'isothermal flow requires single-phase materials'
      return
    end if
    call this%matl_map%init(matl_model, stat, errmsg)
    if (stat /= 0) return

    simple_default = .false.
    if (matl_model%nmatl_real == 1 .and. matl_model%nphase_real == 1 .and. .not.matl_model%have_void) &
      simple_default = matl_model%is_fluid(1)
    algorithm = 'geometric'
    tracking_cutoff = 1.0e-6_r8
    tracking_subcycles = 4
    if (simple_default) algorithm = 'simple'
    if (params%is_sublist('volume-tracking')) then
      tracking_params => params%sublist('volume-tracking')
      call tracking_params%get('algorithm', algorithm, default=algorithm, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // tracking_params%path() // ': ' // errmsg
        return
      end if
      call tracking_params%get('cutoff', tracking_cutoff, default=tracking_cutoff, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // tracking_params%path() // ': ' // errmsg
        return
      end if
      if (tracking_cutoff <= 0.0_r8 .or. tracking_cutoff >= 1.0_r8) then
        stat = 1
        errmsg = 'processing ' // tracking_params%path() // ': "cutoff" must be in (0,1)'
        return
      end if
      call tracking_params%get('subcycles', tracking_subcycles, default=tracking_subcycles, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // tracking_params%path() // ': ' // errmsg
        return
      end if
      if (tracking_subcycles < 1) then
        stat = 1
        errmsg = 'processing ' // tracking_params%path() // ': "subcycles" must be at least one'
        return
      end if
      call this%matl_map%set_priority(tracking_params, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // tracking_params%path() // ': ' // errmsg
        return
      end if
    end if
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
    write(message, '(a,es11.4,a)') 'Using Courant number ', courant_number, '.'
    call env%simlog%info(trim(message))
    if (.not.params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'requires a "projection-solver" sublist'
      return
    end if
    projection_params => params%sublist('projection-solver')
    call env%simlog%info('Using ' // trim(algorithm) // ' volume tracking.')
    nrealfluid = this%matl_map%num_real_fluid()
    nfluid = this%matl_map%num_fluid()
    nmat = this%matl_map%num_material()
    allocate(priority(nmat))
    call this%matl_map%get_priority(priority)
    allocate(phase_ids(nrealfluid))
    call this%matl_map%get_real_fluid_phase_ids(phase_ids)
    call model%init_material(matl_model, phase_ids, stat, errmsg, nfluid=nfluid)
    if (stat /= 0) return
    if (model%inviscid) then
      call this%mechanics%init(env, model, projection_params=projection_params, courant_number=courant_number, &
          stat=stat, errmsg=errmsg)
    else
      if (.not.params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a "momentum-solver" sublist'
        return
      end if
      momentum_params => params%sublist('momentum-solver')
      call this%mechanics%init(env, model, momentum_params, projection_params, courant_number, stat, errmsg)
    end if
    if (stat /= 0) return
    allocate(this%vfrac(nmat,model%mesh%ncell))
    this%vfrac = 0.0_r8
    this%vfrac(1,:) = 1.0_r8
    call this%material_transport%init(env, model%mesh, nrealfluid, nfluid, nmat, algorithm, priority, tracking_cutoff, &
        tracking_subcycles)
    call configure_inflow_material(this, model%bc, stat, errmsg)

  contains

    subroutine configure_inflow_material(this, bc, stat, errmsg)
      class(t2d_flow_solver), intent(inout) :: this
      type(t2d_flow_bc), intent(in) :: bc
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      integer :: i, slot
      integer, allocatable :: assigned(:)

      stat = 0
      if (.not.allocated(bc%inflow_material)) then
        errmsg = ''
        return
      end if
      allocate(assigned(model%mesh%nface_onP), source=0)
      do i = 1, size(bc%inflow_material)
        slot = this%matl_map%slot_index(bc%inflow_material(i)%name)
        if (slot == 0 .or. slot > this%matl_map%num_fluid()) then
          stat = 1
          errmsg = 'invalid flow inflow material: "' // bc%inflow_material(i)%name // '"'
          return
        end if
        if (any(assigned(bc%inflow_material(i)%face) /= 0 .and. &
            assigned(bc%inflow_material(i)%face) /= slot)) then
          stat = 1
          errmsg = 'conflicting flow inflow materials on a boundary face'
          return
        end if
        assigned(bc%inflow_material(i)%face) = slot
        call this%material_transport%set_inflow_material(slot, bc%inflow_material(i)%face)
      end do
      errmsg = ''
    end subroutine
  end subroutine


  subroutine set_volume_fractions(this, vfrac)
    class(t2d_flow_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    ASSERT(size(vfrac,1) == size(this%vfrac,1))
    ASSERT(size(vfrac,2) == size(this%vfrac,2))
    this%vfrac = vfrac
    call this%mechanics%set_volume_fractions(vfrac)
  end subroutine


  subroutine set_initial_material_state(this, vfrac, temperature)
    class(t2d_flow_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%mechanics%set_initial_material_state(vfrac, temperature)
    ASSERT(size(vfrac,1) == size(this%vfrac,1))
    ASSERT(size(vfrac,2) == size(this%vfrac,2))
    this%vfrac = vfrac
  end subroutine


  subroutine get_reduced_volume_fractions(this, matl_dist, vfrac)
    class(t2d_flow_solver), intent(in) :: this
    type(material_distribution), intent(in) :: matl_dist
    real(r8), allocatable, intent(out) :: vfrac(:,:)

    allocate(vfrac(size(this%vfrac,1),size(this%vfrac,2)))
    call this%matl_map%get_reduced_volume_fractions(matl_dist, vfrac)
  end subroutine


  !! Update the simulation-owned material distribution from the current flow
  !! distribution before it is used for output or by another physics model.
  subroutine update_material_distribution(this, matl_dist)
    class(t2d_flow_solver), intent(in) :: this
    type(material_distribution), intent(inout) :: matl_dist

    call this%matl_map%put_reduced_volume_fractions(this%vfrac, matl_dist)
  end subroutine


  subroutine set_buoyancy_temperature(this, temperature)
    class(t2d_flow_solver), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    call this%mechanics%set_buoyancy_temperature(temperature)
  end subroutine

  !! Set STATE from an input velocity. The common initial-condition solver
  !! projects the velocity and computes an initial pressure with its temporary
  !! Stokes step, as mainline does when it omits initial momentum transport.
  subroutine set_initial_state(this, env, time, dt, velocity, stat)
    class(t2d_flow_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    call this%mechanics%set_initial_state(env, time, dt, velocity, stat)
    if (stat /= 0) return
    this%nstep = 0_int64
  end subroutine


  !! Return no-copy views of the accepted cell-centered pressure and velocity.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(t2d_flow_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%mechanics%get_cell_flow_soln(pressure, velocity)
  end subroutine


  !! Return a no-copy view of the full-local flow-equation mask.
  subroutine get_cell_flow_active(this, active)
    class(t2d_flow_solver), target, intent(in) :: this
    logical, pointer, intent(out) :: active(:)

    call this%mechanics%get_cell_flow_active(active)
  end subroutine


  !! Return a no-copy view of the accepted face-normal velocity.
  subroutine get_face_velocity(this, velocity)
    class(t2d_flow_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%mechanics%get_face_velocity(velocity)
  end subroutine


  !! Advance STATE from T_N to T_NP1. The time step is derived from the two
  !! endpoint times so callers retain exact target times. This is the
  !! isothermal wrapper: it first obtains material transport from the old face
  !! velocity and then advances momentum and pressure.
  subroutine step(this, env, t_n, t_np1, stat, errmsg)
    class(t2d_flow_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg

    real(r8), pointer :: vfrac_trial(:,:), face_velocity(:)

    call env%timer%start('flow/material-transport')
    call this%mechanics%get_face_velocity(face_velocity)
    call this%material_transport%advance(env, t_n, t_np1, face_velocity, this%vfrac)
    call this%material_transport%get_trial_volume_fractions(vfrac_trial)
    call env%timer%stop('flow/material-transport')
    call this%mechanics%set_volume_fractions(vfrac_trial)
    if (.not.this%model%unsteady_stokes) then
      call this%mechanics%advance_momentum(env, t_n, t_np1, stat, errmsg, &
          this%material_transport%flux_volumes(:this%matl_map%num_real_fluid(),:))
    else
      call this%mechanics%advance_momentum(env, t_n, t_np1, stat, errmsg)
    end if
    if (stat /= 0) then
      call this%mechanics%reject_step()
      call this%mechanics%set_volume_fractions(this%vfrac)
      return
    end if
    this%vfrac = vfrac_trial
    call this%mechanics%commit_step()
    this%nstep = this%nstep + 1_int64
  end subroutine


  integer(int64) function num_steps(this)
    class(t2d_flow_solver), intent(in) :: this

    num_steps = this%nstep
  end function

  !! Declare the temporal scalar fields published by this solver.
  !! These fields are updated at each requested solution output and written
  !! by the simulation's output writer.
  subroutine init_temporal_output(this, data)
    class(t2d_flow_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Set the current values of the temporal scalar fields published by this
  !! solver.
  subroutine set_temporal_output(this, data)
    class(t2d_flow_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Return the maximum step size requested by the flow mechanics for the old
  !! face-normal velocity.
  function courant_time_step(this) result(dt)
    class(t2d_flow_solver), intent(in) :: this
    real(r8) :: dt

    dt = this%mechanics%courant_time_step()
  end function

end module t2d_flow_solver_type
