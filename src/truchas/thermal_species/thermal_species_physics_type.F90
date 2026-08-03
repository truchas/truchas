!!
!! THERMAL_SPECIES_PHYSICS_TYPE
!!
!! This module defines the package-level object for the heat/species transport
!! physics. It owns the common mesh/material state, selects the concrete
!! thermal/species solver variant needed by the active physics, and provides
!! the public operations used by the driver to advance, restart, and query the
!! solution.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module thermal_species_physics_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_logging_services
  use matl_mesh_func_type
  use mesh_interop
  use truchas_timers
  use truchas_logging_services
  use unstr_mesh_type
  use enthalpy_advector_class
  use conc_advector_type
  use parameter_list_type
  use thermal_species_solver_class
  use htsd_solver_type, only: htsd_solver
  use alloy_solver_type, only: alloy_solver
  use fht_solver_type, only: fht_solver
  use ht_solver_type, only: ht_solver
  use sd_solver_type, only: sd_solver
  implicit none
  private

  type, public :: thermal_species_physics
    !! Problem characteristics.
    logical :: have_heat_transport = .false.
    logical :: have_species_transport = .false.
    logical :: have_em_heat = .false.
    logical :: have_fluid_flow = .false.
    logical :: have_phase_change = .false.
    logical :: have_void = .false.
    logical :: have_advected_heat = .false.
    integer :: num_species = 0
    !! The mesh and material mesh function.
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(matl_mesh_func) :: mmf
    !! Selected solver.
    class(thermal_species_solver), allocatable :: solver
    class(enthalpy_advector), allocatable :: hadv
    type(conc_advector), allocatable :: cadv(:)
    real(r8) :: cutvof
  contains
    procedure :: init
    procedure :: has_heat_transport
    procedure :: has_species_transport
    procedure :: get_num_species
    procedure :: step
    procedure :: commit_step
    procedure :: restart
    procedure :: set_initial_state

    procedure :: get_temp
    procedure :: get_enthalpy
    procedure :: get_phi
    procedure :: get_temp_grad
    procedure :: get_face_temp
    procedure :: get_face_temp_view

    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
  end type

  enum, bind(c)
    enumerator :: SOLVER_UNSET = 0
    enumerator :: SOLVER_HTSD
    enumerator :: SOLVER_FHT
    enumerator :: SOLVER_HT
    enumerator :: SOLVER_SD
    enumerator :: SOLVER_ALLOY
  end enum

contains

  logical function has_heat_transport(this)
    class(thermal_species_physics), intent(in) :: this
    has_heat_transport = this%have_heat_transport
  end function

  logical function has_species_transport(this)
    class(thermal_species_physics), intent(in) :: this
    has_species_transport = this%have_species_transport
  end function

  integer function get_num_species(this)
    class(thermal_species_physics), intent(in) :: this
    get_num_species = this%num_species
  end function

  integer function solver_from_name(name)
    character(*), intent(in) :: name
    select case (name)
    case ('ht')
      solver_from_name = SOLVER_HT
    case ('sd')
      solver_from_name = SOLVER_SD
    case ('htsd')
      solver_from_name = SOLVER_HTSD
    case ('alloy')
      solver_from_name = SOLVER_ALLOY
    case default
      INSIST(.false.)
    end select
  end function

  integer function runtime_solver(this, declared_solver, integrator)

    class(thermal_species_physics), intent(in) :: this
    integer, intent(in) :: declared_solver
    character(*), intent(in) :: integrator

    runtime_solver = declared_solver
    if (declared_solver == SOLVER_ALLOY) then
      if (this%have_void) then
        call TLS_fatal ('INIT: alloy solidification does not support VOID material.')
      end if
      if (integrator /= 'adaptive-bdf2') then
        call TLS_fatal ('INIT: alloy solidification requires adaptive-bdf2 stepping.')
      end if
    else if (this%have_void .and. this%have_fluid_flow) then
      !! Transient void uses the special FHT solver, which is only defined for
      !! heat transfer without species transport.
      if (declared_solver /= SOLVER_HT) then
        INSIST(.false.)
      end if
      runtime_solver = SOLVER_FHT
      if (integrator /= 'nonadaptive-bdf1') then
        call TLS_fatal ('INIT: thermal/species system characteristics are incompatible with STEPPING_METHOD choice.')
      end if
    else if (integrator /= 'adaptive-bdf2') then
      call TLS_fatal ('INIT: thermal/species system characteristics are incompatible with STEPPING_METHOD choice.')
    end if

  end function runtime_solver

  subroutine set_solver_capabilities(this, solver)
    class(thermal_species_physics), intent(inout) :: this
    integer, intent(in) :: solver
    select case (solver)
    case (SOLVER_HTSD)
      this%have_heat_transport = .true.
      this%have_species_transport = .true.
    case (SOLVER_FHT, SOLVER_HT)
      this%have_heat_transport = .true.
      this%have_species_transport = .false.
    case (SOLVER_SD)
      this%have_heat_transport = .false.
      this%have_species_transport = .true.
    case (SOLVER_ALLOY)
      this%have_heat_transport = .true.
      !! Alloy solute is an internal phase-change state, not the package's
      !! conventional SPECIES_TRANSPORT field.  Keeping this false avoids the
      !! ordinary species initialization and material-fraction update paths.
      this%have_species_transport = .false.
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine set_species_count(this, solver, solvers_params)

    class(thermal_species_physics), intent(inout) :: this
    integer, intent(in) :: solver
    type(parameter_list), pointer, intent(in) :: solvers_params

    type(parameter_list), pointer :: htsd_params, htsd_model_params, alloy_params
    type(parameter_list), pointer :: sd_params, species_model_params

    select case (solver)
    case (SOLVER_HTSD)
      htsd_params => solvers_params%sublist('htsd')
      htsd_model_params => htsd_params%sublist('model')
      species_model_params => htsd_model_params%sublist('species')
      call species_model_params%get('num-species', this%num_species)
    case (SOLVER_SD)
      sd_params => solvers_params%sublist('species')
      species_model_params => sd_params%sublist('model')
      call species_model_params%get('num-species', this%num_species)
    case (SOLVER_ALLOY)
      alloy_params => solvers_params%sublist('alloy')
      call alloy_params%get('num-comp', this%num_species)
    case default
      this%num_species = 0
    end select
    this%have_species_transport = (this%num_species > 0 .and. solver /= SOLVER_ALLOY)

  end subroutine set_species_count

  subroutine get_declared_solver_params(declared_solver, solvers_params, solver_params)

    integer, intent(in) :: declared_solver
    type(parameter_list), pointer, intent(in) :: solvers_params
    type(parameter_list), pointer, intent(out) :: solver_params

    type(parameter_list), pointer :: ht_params, sd_params, htsd_params, alloy_params

    select case (declared_solver)
    case (SOLVER_HT)
      ht_params => solvers_params%sublist('heat')
      solver_params => ht_params%sublist('solver')
    case (SOLVER_SD)
      sd_params => solvers_params%sublist('species')
      solver_params => sd_params%sublist('solver')
    case (SOLVER_HTSD)
      htsd_params => solvers_params%sublist('htsd')
      solver_params => htsd_params%sublist('solver')
    case (SOLVER_ALLOY)
      alloy_params => solvers_params%sublist('alloy')
      solver_params => alloy_params
    case default
      INSIST(.false.)
    end select

  end subroutine get_declared_solver_params

  subroutine init_solver(this, solver_kind, solvers_params, stat, errmsg)
    !! TARGET starts the lifetime chain for solver components that are later
    !! retained by pointer in the IDAESOL adapter objects.
    class(thermal_species_physics), intent(inout), target :: this
    integer, intent(in) :: solver_kind
    type(parameter_list), pointer, intent(in) :: solvers_params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: ht_params, sd_params, fht_params, htsd_params, alloy_params

    select case (solver_kind)
    case (SOLVER_HTSD)
      allocate(htsd_solver :: this%solver)
      select type (solver => this%solver)
      type is (htsd_solver)
        htsd_params => solvers_params%sublist('htsd')
        call solver%init(this%mesh, this%mmf, htsd_params, stat, errmsg)
      end select

    case (SOLVER_SD)
      allocate(sd_solver :: this%solver)
      select type (solver => this%solver)
      type is (sd_solver)
        sd_params => solvers_params%sublist('species')
        call solver%init(this%mesh, this%mmf, sd_params, stat, errmsg)
      end select

    case (SOLVER_FHT)
      allocate(fht_solver :: this%solver)
      select type (solver => this%solver)
      type is (fht_solver)
        fht_params => solvers_params%sublist('heat')
        call solver%init(this%mesh, this%mmf, fht_params, this%hadv, stat, errmsg)
      end select

    case (SOLVER_HT)
      allocate(ht_solver :: this%solver)
      select type (solver => this%solver)
      type is (ht_solver)
        ht_params => solvers_params%sublist('heat')
        call solver%init(this%mesh, this%mmf, ht_params, stat, errmsg)
      end select

    case (SOLVER_ALLOY)
      allocate(alloy_solver :: this%solver)
      select type (solver => this%solver)
      type is (alloy_solver)
        alloy_params => solvers_params%sublist('alloy')
        call solver%init(this%mesh, this%mmf, alloy_params, stat, errmsg)
      end select

    case default
      INSIST(.false.)
    end select

  end subroutine init_solver

  subroutine init_heat_advector(this, stat, errmsg)

    use enthalpy_advector1_type
    use physics_module, only: flow

    class(thermal_species_physics), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(enthalpy_advector1), allocatable :: hadv1

    this%have_fluid_flow = .false.
    if (flow) then
      this%have_fluid_flow = .true.
      allocate(hadv1)
      call hadv1%init(this%mesh, stat, errmsg)
      if (stat /= 0) return
      call move_alloc(hadv1, this%hadv)
    else
      stat = 0
    end if

  end subroutine init_heat_advector

  subroutine init_species_advectors(this, stat, errmsg)

    class(thermal_species_physics), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i

    stat = 0
    if (this%have_fluid_flow .and. this%have_species_transport) then
      allocate(this%cadv(this%num_species))
      do i = 1, this%num_species
        call this%cadv(i)%init(this%mesh, i, stat, errmsg)
        if (stat /= 0) return
      end do
    end if

  end subroutine init_species_advectors


  subroutine init(this, mesh, params)

    use em_heat_driver, only: em_heat_enabled

    class(thermal_species_physics), intent(inout), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), pointer, intent(in) :: params

    integer :: stat
    integer :: declared_solver
    integer :: solver_type = SOLVER_UNSET
    character(:), allocatable :: declared_solver_name, integrator, errmsg
    type(parameter_list), pointer :: solvers_params, solver_params

    call TLS_info ('')
    call TLS_info ('Initializing thermal/species physics ...')

    !! Common initialization.
    this%mesh => mesh

    call mmf_init (this%mesh, this%mmf, stat, errmsg)
    if (stat /= 0) call TLS_fatal ('INIT: ' // errmsg)

    !TODO? reimplement this
    !call verify_material_compatibility (this%mmf, stat, errmsg)
    !if (stat /= 0) call TLS_fatal ('INIT: ' // trim(errmsg))

    !! Problem attributes
    this%have_em_heat = em_heat_enabled()
    this%have_phase_change = multiphase_problem(this%mmf)
    this%have_void = void_is_present()
    this%have_advected_heat = .false.

    call params%get('declared-solver', declared_solver_name)
    declared_solver = solver_from_name(declared_solver_name)
    if (declared_solver == SOLVER_ALLOY) this%have_phase_change = .false.
    solvers_params => params%sublist('solvers')
    call get_declared_solver_params(declared_solver, solvers_params, solver_params)

    call init_heat_advector(this, stat, errmsg)
    if (stat /= 0) call TLS_fatal('INIT: error initializing enthalpy advection: ' // errmsg)

    call solver_params%get('cutvof', this%cutvof, default=0.0_r8)
    call solver_params%get('integrator', integrator)

    solver_type = runtime_solver(this, declared_solver, integrator)
    call set_solver_capabilities(this, solver_type)
    this%have_advected_heat = this%have_heat_transport .and. this%have_fluid_flow .and. solver_type /= SOLVER_FHT
    call set_species_count(this, solver_type, solvers_params)
    call init_species_advectors(this, stat, errmsg)
    if (stat /= 0) call TLS_fatal('INIT: error initializing scalar advection: ' // errmsg)

    call init_solver(this, solver_type, solvers_params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('INIT: ' // errmsg)

    call TLS_info ('  thermal/species physics initialized')

  contains

    logical function multiphase_problem (mmf)
      use material_model_driver, only: matl_model
      type(matl_mesh_func), intent(in) :: mmf
      integer, allocatable :: list(:)
      integer :: i, p1, p2
      call mmf%get_all_matl(list, drop_void=.true.)
      do i = size(list), 1, -1
        call matl_model%get_matl_phase_index_range(list(i), p1, p2)
        if (p2 > p1) exit
      end do
      multiphase_problem = (i /= 0)
    end function

  end subroutine init


  subroutine step(this, t, hnext, errc)

    class(thermal_species_physics), intent(inout) :: this
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc

    real(r8) :: h, tlast

    call start_timer('Thermal/Species Physics')

    tlast = this%solver%last_time()
    h = t - tlast

    if (this%have_fluid_flow) then
      call update_mmf_from_matl(this%mmf)
      if (this%have_void) call cull_material_fragments(this%mmf, this%cutvof)
    end if

    if (this%have_heat_transport) call update_adv_heat
    if (allocated(this%cadv)) call update_adv_conc

    call this%solver%step(t, hnext, errc)

    call stop_timer('Thermal/Species Physics')

  contains

    subroutine update_adv_heat

      use em_heat_driver, only: em_heat_ptr

      real(r8), allocatable :: q_rate(:)

      if (this%have_em_heat .or. this%have_advected_heat) then
        associate (ncell_onP => this%mesh%ncell_onP)
        allocate(q_rate(ncell_onP))
        !! Joule heat source.
        if (this%have_em_heat) then
          q_rate(:) = em_heat_ptr() * this%mesh%volume(:ncell_onP)
        else
          q_rate = 0.0_r8
        end if
        !! Advected heat source.
        if (this%have_advected_heat) then
          block
            real(r8) :: tcell(ncell_onP), dQ(ncell_onP)
            call this%get_temp(tcell) !TODO: can be a view rather than a copy
            call this%hadv%get_advected_enthalpy(tcell, dQ)
            q_rate = q_rate + dQ/h
          end block
        end if
        call this%solver%set_ext_enthalpy_rate(q_rate)
        deallocate(q_rate)
        end associate
      end if

    end subroutine update_adv_heat

    subroutine update_adv_conc

      integer :: i
      real(r8) :: phi(this%mesh%ncell), dphi(this%mesh%ncell_onP)

      do i = 1, this%num_species
        call this%get_phi(i, phi)
        call this%mesh%cell_imap%gather_offp(phi)
        call this%cadv(i)%get_advected_scalar(tlast, phi, dphi)
        dphi = dphi/h ! turn into a rate
        call this%solver%set_ext_species_rate(i, dphi)
      end do

    end subroutine update_adv_conc

  end subroutine step

  subroutine commit_step(this)

    use zone_module, only: zone

    class(thermal_species_physics), intent(inout) :: this

    call start_timer('Thermal/Species Physics')

    call this%solver%commit_step

    !! Update MATL in contexts that can modify the phase distribution.
    if (this%have_phase_change) call update_phase_fractions
    select type (solver => this%solver)
    type is (alloy_solver)
      block
        real(r8), pointer :: liquid_frac(:)
        call solver%get_liq_frac_view(liquid_frac)
        call update_matl_from_alloy(this%mmf, liquid_frac)
      end block
    end select

    !! Update Truchas data structures to reflect the new HT/SD solution.
    if (this%have_heat_transport) then
      zone%temp_old = zone%temp ! store for dT/dt output
      call this%get_temp(zone%temp)
      call this%get_enthalpy(zone%enthalpy)
    end if

    !! Write some info about the step.
    call this%solver%log_step_stats

    call stop_timer('Thermal/Species Physics')

  contains

    subroutine update_phase_fractions
      real(r8), allocatable :: state(:,:)
      allocate(state(this%mesh%ncell_onP,1))
      call this%get_temp(state(:,1))
      call update_matl_from_mmf(this%mmf, state)
    end subroutine

  end subroutine commit_step

  !! The effect of calling this subroutine is to restart or reset the solver so
  !! that its subsequent numerical behavior is as if it was starting integration
  !! from an initial state equal to the current state.  This mainly means
  !! dropping any previous solution history in the BDF2 solver and recomputing
  !! an approximation to the initial state time derivative.  The argument DT is
  !! used to compute that time derivative and is best chosen equal to the next
  !! time step size; however it has no bearing on the next step size used.

  subroutine restart(this, dt)
    class(thermal_species_physics), intent(inout) :: this
    real(r8), intent(in) :: dt
    call this%solver%restart(dt)
  end subroutine

  subroutine set_initial_state(this, t, dt, temp, conc)

    class(thermal_species_physics), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)

    !! Permute the cell temperature array to the solver ordering.
    if (this%have_heat_transport) then
      ASSERT(present(temp))
      ASSERT(size(temp) == this%mesh%ncell_onP)
    end if

    !! Permute the cell concentration array to the solver ordering.
    if (this%have_species_transport) then
      INSIST(present(conc))
      ASSERT(size(conc,dim=1) == this%mesh%ncell_onP)
      ASSERT(size(conc,dim=2) == this%num_species)
    end if

    !! Set the initial state in the appropriate solver. Solver init builds
    !! structural data only; time-dependent state, including moving view-factor
    !! state, is established here before any moving-VF event can advance it.
    call this%solver%set_initial_state(t, dt, temp, conc)

  end subroutine set_initial_state

  subroutine get_temp(this, array)
    class(thermal_species_physics), intent(in) :: this
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_heat_transport)
    call this%solver%get_cell_temp_copy(array)
  end subroutine

  subroutine get_enthalpy(this, array)
    class(thermal_species_physics), intent(in) :: this
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_heat_transport)
    call this%solver%get_cell_heat_copy(array)
  end subroutine

  subroutine get_phi(this, n, array)
    class(thermal_species_physics), intent(in) :: this
    integer,  intent(in)  :: n
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_species_transport)
    call this%solver%get_cell_conc_copy(n, array)
  end subroutine get_phi

  subroutine get_temp_grad(this, array)
    class(thermal_species_physics), intent(inout) :: this
    real(r8), intent(inout) :: array(:,:)
    ASSERT(size(array,dim=2) >= this%mesh%ncell_onP)
    ASSERT(size(array,dim=1) == 3)
    ASSERT(this%have_heat_transport)
    call this%solver%get_cell_temp_grad(array(:,:this%mesh%ncell_onP))
  end subroutine get_temp_grad

  !! Copy the current face temperatures on the new distributed mesh.
  subroutine get_face_temp(this, array)
    class(thermal_species_physics), intent(in) :: this
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%nface_onP)
    ASSERT(this%have_heat_transport)
    call this%solver%get_face_temp_copy(array)
  end subroutine get_face_temp

  !! Get reference to the current face temperatures on the new distributed mesh.
  subroutine get_face_temp_view(this, view)
    class(thermal_species_physics), intent(in) :: this
    real(r8), pointer :: view(:)
    ASSERT(this%have_heat_transport)
    call this%solver%get_face_temp_view(view)
  end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CULL_MATERIAL_FRAGMENTS
 !!
 !! This procedure removes all material from cells with a small non-void
 !! volume fraction.  Where the void volume fraction exceeds 1 - THRESHOLD,
 !! the void volume fraction is set to 1 and the volume fraction for all other
 !! materials is set to 0.  Note that this doesn't remove materials with a
 !! small volume fraction unless the total non-void volume fraction is less
 !! than THRESHOLD.
 !!

  subroutine cull_material_fragments (mmf, threshold, culled)

    use string_utilities, only: i_to_c

    type(matl_mesh_func), intent(inout) :: mmf
    real(r8), intent(in) :: threshold
    logical, intent(out), optional :: culled

    integer :: n, j
    integer, pointer :: matID(:)
    real(r8), pointer :: vfrac(:,:)
    integer :: cell_count

    cell_count = 0
    do n = 1, mmf%num_reg()
      matID => mmf%reg_matl(n)
      if (matID(1) /= 0) cycle  ! no void in this region
      vfrac => mmf%reg_vol_frac(n)
      if (associated(vfrac)) then ! multi-material region
        do j = 1, size(vfrac,dim=1)
          if (vfrac(j,1) == 1.0_r8) then  ! ensure all other volume fractions are 0
            vfrac(j,2:) = 0.0_r8
          else if (vfrac(j,1) >= 1.0_r8 - threshold) then ! make the cell totally void
            cell_count = cell_count + 1
            vfrac(j,1) = 1.0_r8
            vfrac(j,2:) = 0.0_r8
          end if
        end do
      end if
    end do
    cell_count = global_sum(cell_count)
    if (cell_count > 0) &
        call TLS_info ('thermal/species: culled material fragments from ' // i_to_c(cell_count) // ' cells.')
    if (present(culled)) culled = (cell_count > 0)

  end subroutine cull_material_fragments

  !! Called in response to moving view-factor events from the global event
  !! queue.  This advances the loaded view-factor interval; ordinary
  !! time-dependent interpolation is handled inside the radiation problem.
  !! Eventually this action could be registered directly with the event.

  subroutine update_moving_vf(this)

    class(thermal_species_physics), intent(inout) :: this

    select type (solver => this%solver)
    type is (htsd_solver)
      call solver%update_moving_vf
    type is (fht_solver)
      call solver%update_moving_vf
    type is (ht_solver)
      call solver%update_moving_vf
    type is (sd_solver)
      !! Species-only transport has no view-factor radiation state.
    type is (alloy_solver)
      !! Alloy solidification has no view-factor radiation state.
    class default
      INSIST(.false.)
    end select

  end subroutine


  subroutine add_moving_vf_events(this, eventq, rank)

    use sim_event_queue_type

    class(thermal_species_physics), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank

    select type (solver => this%solver)
    type is (htsd_solver)
      call solver%add_moving_vf_events(eventq, rank)
    type is (fht_solver)
      call solver%add_moving_vf_events(eventq, rank)
    type is (ht_solver)
      call solver%add_moving_vf_events(eventq, rank)
    type is (sd_solver)
      !! Species-only transport has no view-factor radiation events.
    type is (alloy_solver)
      !! Alloy solidification has no view-factor radiation events.
    class default
      INSIST(.false.)
    end select

  end subroutine

end module thermal_species_physics_type
