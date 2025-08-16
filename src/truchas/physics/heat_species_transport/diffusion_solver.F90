!!
!! DIFFUSION_SOLVER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Markus Berndt <berndt@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module diffusion_solver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use diffusion_solver_data
  use mesh_manager
  use parallel_communication
  use truchas_logging_services
  use mfd_disc_type
  use matl_mesh_func_type
  use string_utilities, only: i_to_c
  use mesh_interop
  use ht_model_type
  use ht_solver_type
  use alloy_model_type
  use alloy_solver_type
  use FHT_model_type
  use FHT_solver_type
  use HTSD_model_type, only: htsd_model, htsd_model_delete
  use HTSD_solver_type
  use truchas_timers
  use truchas_logging_services
  use unstr_mesh_type
  use enthalpy_advector_class
  use conc_advector_type
  use parameter_list_type
  use rad_problem_type, only: vf_event
  implicit none
  private

  public :: ds_init, ds_step, ds_accept
  public :: read_ds_namelists
  public :: ds_set_initial_state, ds_restart
  public :: ds_delete
  public :: update_moving_vf, add_moving_vf_events, vf_event

  !! These return cell-centered results relative to the old Truchas mesh.
  public :: ds_get_temp, ds_get_enthalpy, ds_get_phi
  public :: ds_get_temp_grad
  public :: ds_get_alloy_conc_view, ds_get_alloy_liquid_conc, ds_get_alloy_solid_conc

  !! These return results relative to the new distributed mesh.
  public :: ds_get_cell_temp, ds_get_face_temp
  public :: ds_get_face_temp_view

  type :: ds_driver
    !! Problem characteristics.
    logical :: have_heat_transfer = .false.
    logical :: have_species_diffusion = .false.
    logical :: have_joule_heat = .false.
    logical :: have_fluid_flow = .false.
    logical :: have_phase_change = .false.
    logical :: have_void = .false.
    !! The mesh, discretization, and material mesh function.
    type(unstr_mesh), pointer :: mesh => null()
    type(mfd_disc), pointer :: disc => null()
    type(matl_mesh_func), pointer :: mmf => null()
    !! Solver selections
    integer :: solver_type = 0
    !! The standard model and solver.
    type(HTSD_model),  pointer :: mod1 => null()
    type(HTSD_solver), pointer :: sol1 => null()
    !! The special model and solver that works with transient void.
    type(FHT_model),  pointer :: mod2  => null()
    type(FHT_solver), pointer :: sol2 => null()
    !! New prototype HT-only solver
    type(ht_model),  pointer :: mod3 => null()
    type(ht_solver), pointer :: sol3 => null()
    type(alloy_model),  pointer :: mod4 => null()
    type(alloy_solver), pointer :: sol4 => null()
    class(enthalpy_advector), allocatable :: hadv
    type(conc_advector), allocatable :: cadv(:)
    real(r8) :: cutvof
    type(parameter_list) :: ds_params, bc_params, thermal_source_params, &
        species_bc_params, species_source_params, alloy_params
  end type ds_driver
  type(ds_driver), save, target :: this

  integer, parameter :: SOLVER1 = 1 ! the standard solver
  integer, parameter :: SOLVER2 = 2 ! special solver that works with transient void
  integer, parameter :: SOLVER3 = 3 ! new HT prototype
  integer, parameter :: SOLVER4 = 4 ! prototype alloy solidification

contains

  subroutine read_ds_namelists (lun)

    use thermal_bc_namelist
    use thermal_source_namelist
    use species_bc_namelist
    use species_source_namelist
    use enclosure_radiation_namelist
    use diffusion_solver_namelist
    use diffusion_solver_data, only: ds_sys_type, num_species, void_temperature
    use physics_module, only: alloy_solidification
    use alloy_namelist

    integer, intent(in) :: lun
    type(parameter_list), pointer :: plist

    call ds_data_init
    if (alloy_solidification) then
      call read_diffusion_solver_namelist(lun, ds_sys_type, num_species, this%alloy_params)
      call read_alloy_namelist(lun, this%alloy_params)
      plist => this%alloy_params%sublist('sources')
      call read_thermal_source_namelists(lun, plist)
      plist => this%alloy_params%sublist('bc')
      call read_thermal_bc_namelists(lun, plist)
    else
      call read_diffusion_solver_namelist(lun, ds_sys_type, num_species, this%ds_params)
      call read_thermal_bc_namelists(lun, this%bc_params)
      call read_thermal_source_namelists(lun, this%thermal_source_params)
      call read_species_bc_namelists(lun, this%species_bc_params)
      call read_species_source_namelists(lun, this%species_source_params)
      call read_enclosure_radiation_namelists(lun)
      call this%ds_params%get('void-temperature', void_temperature)
    end if

  end subroutine read_ds_namelists


  subroutine ds_step(h, hnext, errc)

    real(r8), intent(in)  :: h
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc

    real(r8) :: t

    call start_timer('Diffusion Solver')

    if (this%have_fluid_flow) then
      call update_mmf_from_matl(this%mmf)
      if (this%have_void) call cull_material_fragments(this%mmf, this%cutvof)
    end if

    if (this%have_heat_transfer) call update_adv_heat
    if (allocated(this%cadv)) then
      select case (this%solver_type)
      case (SOLVER4)
        t = this%sol4%last_time()
      case default ! necessarily using SOLVER1
        t = HTSD_solver_last_time(this%sol1)
      end select
      call update_adv_conc
    end if

    select case (this%solver_type)
    case (SOLVER1)  ! HT/SD solver with static void
      t = h + HTSD_solver_last_time(this%sol1)
      call HTSD_solver_advance_state(this%sol1, t, hnext, errc)
    case (SOLVER2)  ! HT solver with transient void
      t = h + FHT_solver_last_time(this%sol2)
      call FHT_solver_advance_state(this%sol2, t, errc)
      hnext = merge(h/2, huge(1.0_r8), (errc /= 0))
    case (SOLVER3)  ! new HT prototype
      t = h + this%sol3%last_time()
      call this%sol3%step(t, hnext, errc)
    case (SOLVER4)  ! alloy solidification
      t = h + this%sol4%last_time()
      call this%sol4%step(t, hnext, errc)
!block
!  real(r8), pointer :: view(:)
!  call this%sol4%get_cell_temp_view(view)
!  print '(a,*(f8.4))', 'temp:', view
!  call this%sol4%get_liq_frac_view(view)
!  print '(a,*(f8.4))', 'lfrac:', view
!end block
    case default
      INSIST(.false.)
    end select

    call stop_timer('Diffusion Solver')

  contains

    subroutine update_adv_heat

      use ih_driver, only: joule_power_density

      real(r8), allocatable :: q_ds(:)

      if (this%have_joule_heat .or. (this%have_fluid_flow .and. this%solver_type /= SOLVER2)) then
        allocate(q_ds(this%mesh%ncell))
        !! Joule heat source.
        if (this%have_joule_heat) then
          q_ds(:this%mesh%ncell_onP) = joule_power_density()
          call this%mesh%cell_imap%gather_offp(q_ds)
        else
          q_ds = 0.0_r8
        end if
        !! Advected heat source.
        if (this%have_fluid_flow .and. this%solver_type /= SOLVER2) then
          block
            real(r8) :: tcell(this%mesh%ncell_onP), dQ(this%mesh%ncell)
            call ds_get_temp(tcell) !TODO: can be a view rather than a copy
            call this%hadv%get_advected_enthalpy(tcell, dQ(:this%mesh%ncell_onP))
            call this%mesh%cell_imap%gather_offp(dQ)
            q_ds = q_ds + (dQ / (h * this%mesh%volume))
          end block
        end if
        select case (this%solver_type)
        case (SOLVER1)  ! HT/SD solver with static void
          call this%mod1%set_ht_adv_source(q_ds)
        case (SOLVER2)  ! HT solver with transient void
          call this%mod2%set_ht_adv_source(q_ds)
        case (SOLVER3)
          call this%mod3%set_ht_adv_source(q_ds)
        case (SOLVER4)
          call this%sol4%set_adv_heat(q_ds)
        end select
        deallocate(q_ds)
      end if

    end subroutine update_adv_heat

    subroutine update_adv_conc

      integer :: i
      real(r8) :: phi(this%mesh%ncell), dphi(this%mesh%ncell)

      select case (this%solver_type)
      case (SOLVER1)
        do i = 1, num_species
          call ds_get_phi(i, phi)
          call this%mesh%cell_imap%gather_offp(phi)
          call this%cadv(i)%get_advected_scalar(t, phi, dphi)
          call this%mesh%cell_imap%gather_offp(dphi)
          dphi = dphi / (h*this%mesh%volume) ! turn into a source (per unit volume-time)
          call this%mod1%set_sd_adv_source(i, dphi)
        end do
      case (SOLVER4)
        do i = 1, this%sol4%num_comp
          call this%sol4%get_C_liq(i, phi)
          call this%mesh%cell_imap%gather_offp(phi)
          call this%cadv(i)%get_advected_scalar(t, phi, dphi)
          call this%mesh%cell_imap%gather_offp(dphi)
          dphi = dphi / (h*this%mesh%volume) ! turn into a source (per unit volume-time)
          call this%sol4%set_cdot(i, dphi)
        end do
      case default
        INSIST(.false.)
      end select

    end subroutine update_adv_conc

  end subroutine ds_step

  subroutine ds_accept

    use zone_module, only: zone

    real(r8) :: hlast
    integer :: counters(6)
    character(len=80) :: message(2)
    real(r8), pointer :: state(:,:)

    call start_timer('Diffusion Solver')

    select case (this%solver_type)
    case (SOLVER1)  ! HT/SD solver with static void
      call HTSD_solver_commit_pending_state(this%sol1)
    case (SOLVER2)  ! HT solver with transient void
      call FHT_solver_commit_pending_state(this%sol2)
    case (SOLVER3)  ! new HT prototype
      call this%sol3%commit_pending_state
    case (SOLVER4)  ! new HT prototype
      call this%sol4%commit_pending_state
    case default
      INSIST(.false.)
    end select

    !! Update MATL in contexts that can modify the phase distribution.
    if (this%solver_type == SOLVER4) then
      !FIXME: AWFUL UGLY HACK -- ASSUMES A SINGLE 2_PHASE MATERIAL
      block
        use legacy_matl_api, only: define_matl
        real(r8) :: vof(2,this%mesh%ncell_onp)
        real(r8), pointer :: lfrac(:)
        integer ::i
        call this%sol4%get_liq_frac_view(lfrac)
        vof(1,:) = 1 - lfrac(:this%mesh%ncell_onp)
        vof(2,:) = lfrac(:this%mesh%ncell_onp)
        call define_matl(vof)
        INSIST(allocated(zone%phi))
        INSIST(size(zone%phi,dim=1) == this%sol4%num_comp)
        do i = 1, this%sol4%num_comp
          call this%sol4%get_C_liq(i, zone%phi(i,:))
        end do
      end block
    else if (this%have_phase_change) then
      call create_state_array(state)
      call update_matl_from_mmf(this%mmf, state)
      deallocate(state)
    end if

    !! Update Truchas data structures to reflect the new HT/SD solution.
    if (this%have_heat_transfer) then
      zone%temp_old = zone%temp ! store for dT/dt output
      call ds_get_temp(zone%temp)
      call ds_get_enthalpy(zone%enthalpy)
    end if

    !! Write some info about the step.
    select case (this%solver_type)
    case (SOLVER1)
      hlast = HTSD_solver_last_step_size(this%sol1)
      call HTSD_solver_get_stepping_stats (this%sol1, counters)
      write(message,1) hlast, counters(1:2), counters(4:6)
    case (SOLVER2)
      hlast = FHT_solver_last_step_size(this%sol2)
      call FHT_solver_get_stepping_stats (this%sol2, counters)
      write(message,2) hlast, counters(2:4)
    case (SOLVER3)
      hlast = this%sol3%last_step_size()
      call this%sol3%get_stepping_stats(counters)
      write(message,1) hlast, counters(1:2), counters(4:6)
    case (SOLVER4)
      hlast = this%sol4%last_step_size()
      call this%sol4%get_stepping_stats(counters)
      write(message,1) hlast, counters(1:2), counters(4:6)
    case default
      INSIST(.false.)
    end select
    1 format(/,'DS: dt=',es9.3,', NFUN:NPC=',i7.7,':',i5.5,', NNR:NNF:NSR=',3(i4.4,:,':'))
    2 format(/,'DS: dt=',es9.3,', NFUN:NPC:NPA=',3(i7.7,:,':'))
    call TLS_info (message)

    call stop_timer('Diffusion Solver')

  contains

    subroutine create_state_array(state)
      real(r8), pointer :: state(:,:)
      integer :: n, i
      n = 1; if (this%have_heat_transfer) n = 0
      allocate(state(this%mesh%ncell_onP,n:num_species))
      select case (this%solver_type)
      case (SOLVER1)
        !TODO! HTSD_model_new_state_array returns a ncells state, whereas the user
        !TODO! of this routine (update_matl_from_mmf) expects a ncells_onP state.
        !TODO! Is it worthwhile changing that routine, and doing unnecessary comm,
        !TODO! to use the HTSD_model_new_state_array routine instead?
        if (this%have_heat_transfer) call HTSD_solver_get_cell_temp_copy(this%sol1, state(:,0))
        do i = 1, num_species
          call HTSD_solver_get_cell_conc_copy(this%sol1, i, state(:,i))
        end do
      case (SOLVER2)
        call FHT_solver_get_cell_temp_copy(this%sol2, state(:,0))
      case (SOLVER3)
        call this%sol3%get_cell_temp_copy(state(:,0))
      case (SOLVER4)
        call this%sol4%get_cell_temp_copy(state(:,0))
      case default
        INSIST(.false.)
      end select
    end subroutine

  end subroutine ds_accept

  subroutine ds_get_temp (array)
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_temp_copy (this%sol1, array)
    case (SOLVER2)
      call FHT_solver_get_cell_temp_copy (this%sol2, array)
    case (SOLVER3)
      call this%sol3%get_cell_temp_copy(array)
    case (SOLVER4)
      call this%sol4%get_cell_temp_copy(array)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_enthalpy (array)
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_heat_copy (this%sol1, array)
    case (SOLVER2)
      call FHT_solver_get_cell_heat_copy (this%sol2, array)
    case (SOLVER3)
      call this%sol3%get_cell_heat_copy(array)
    case (SOLVER4)
      call this%sol4%get_cell_heat_copy(array)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_phi (n, array)
    integer,  intent(in)  :: n
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_species_diffusion)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_conc_copy (this%sol1, n, array)
    case default
      INSIST(.false.)
    end select
  end subroutine ds_get_phi

  subroutine ds_get_alloy_conc_view(view)
    real(r8), pointer, intent(out) :: view(:,:)
    select case (this%solver_type)
    case (SOLVER4)
      call this%sol4%get_cell_conc_view(view)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_alloy_liquid_conc(n, array)
    integer, intent(in) :: n
    real(r8), intent(inout) :: array(:)
    select case (this%solver_type)
    case (SOLVER4)
      call this%sol4%get_C_liq(n, array)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_alloy_solid_conc(n, array)
    integer, intent(in) :: n
    real(r8),intent(inout) :: array(:)
    select case (this%solver_type)
    case (SOLVER4)
      call this%sol4%get_C_sol(n, array)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_temp_grad (array)
    real(r8), intent(inout) :: array(:,:)
    ASSERT(size(array,dim=2) >= this%mesh%ncell_onP)
    ASSERT(size(array,dim=1) == 3)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_temp_grad (this%sol1, array(:,:this%mesh%ncell_onP))
    case (SOLVER2)
      call FHT_solver_get_cell_temp_grad (this%sol2, array(:,:this%mesh%ncell_onP))
    case (SOLVER3)
      call this%sol3%get_cell_temp_grad(array(:,:this%mesh%ncell_onP))
    case (SOLVER4)
      call this%sol4%get_cell_temp_grad(array(:,:this%mesh%ncell_onP))
    case default
      INSIST(.false.)
    end select
    array(:,this%mesh%ncell_onP+1:) = 0.0_r8  ! gap elements
  end subroutine ds_get_temp_grad

  !! Get reference to the current cell temperatures on the new distributed mesh.
  subroutine ds_get_cell_temp (array)
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%ncell_onP)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_temp_copy (this%sol1, array)
    case (SOLVER2)
      call FHT_solver_get_cell_temp_copy (this%sol2, array)
    case (SOLVER3)
      call this%sol3%get_cell_temp_copy(array)
    case (SOLVER4)
      call this%sol4%get_cell_temp_copy(array)
    case default
      INSIST(.false.)
    end select
  end subroutine ds_get_cell_temp

  !! Get reference to the current face temperatures on the new distributed mesh.
  subroutine ds_get_face_temp (array)
    real(r8), intent(inout) :: array(:)
    ASSERT(size(array) >= this%mesh%nface_onP)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_face_temp_copy (this%sol1, array)
    case (SOLVER2)
      call FHT_solver_get_face_temp_copy (this%sol2, array)
    case (SOLVER3)
      call this%sol3%get_face_temp_copy(array)
    case (SOLVER4)
      call this%sol4%get_face_temp_copy(array)
    case default
      INSIST(.false.)
    end select
  end subroutine ds_get_face_temp

  !! Get reference to the current face temperatures on the new distributed mesh.
  subroutine ds_get_face_temp_view (view)
    real(r8), pointer :: view(:)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_face_temp_view (this%sol1, view)
    case (SOLVER2)
      call FHT_solver_get_face_temp_view (this%sol2, view)
    case (SOLVER3)
      call this%sol3%get_face_temp_view(view)
    case (SOLVER4)
      call this%sol4%get_face_temp_view(view)
    case default
      INSIST(.false.)
    end select
  end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DS_INIT
 !!
  subroutine ds_init (tinit)

    use ih_driver, only: ih_enabled
    use FHT_model_factory
    use FHT_solver_factory
    use HTSD_model_factory
    use HTSD_solver_factory
    use ht_model_factory
    use ht_solver_factory
    use alloy_solver_factory
    use physics_module, only: flow, alloy_solidification
    use enthalpy_advector1_type
    use thermal_bc_factory1_type
    use species_bc_factory1_type
    use thermal_source_factory_type
    use species_source_factory_type
    use physical_constants, only: stefan_boltzmann, absolute_zero

    real(r8), intent(in) :: tinit

    integer :: stat, i, n
    character(len=200) :: errmsg
    type(enthalpy_advector1), allocatable :: hadv1
    character(:), allocatable :: integrator, errmsg2

    call TLS_info ('')
    call TLS_info ('Initializing diffusion solver ...')

    !! Common initialization.
    this%mesh => unstr_mesh_ptr(mesh_name)
    INSIST(associated(this%mesh))

    allocate(this%disc)
    call this%disc%init (this%mesh, use_new_mfd=.true.)

    allocate(this%mmf)
    call mmf_init (this%mesh, this%mmf, stat, errmsg2)
    if (stat /= 0) call TLS_fatal ('DS_INIT: ' // errmsg2)

    !TODO? reimplement this
    !call verify_material_compatibility (this%mmf, stat, errmsg)
    !if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))

    !! Problem attributes
    this%have_heat_transfer = heat_eqn
    this%have_species_diffusion = (num_species > 0)
    this%have_joule_heat = ih_enabled()
    this%have_fluid_flow = .false.
    this%have_phase_change = multiphase_problem(this%mmf)
    this%have_void = void_is_present()

    if (flow) then
      this%have_fluid_flow = .true.
      allocate(hadv1)
      call hadv1%init(this%mesh, stat, errmsg2)
      if (stat /= 0) call TLS_fatal('DS_INIT: error initializing enthalpy advection: ' // errmsg2)
      call move_alloc(hadv1, this%hadv)
      if (num_species > 0) then
        allocate(this%cadv(num_species))
        do i = 1, num_species
          call this%cadv(i)%init(this%mesh, i, stat, errmsg2)
          if (stat /= 0) call TLS_fatal('DS_INIT: error initializing concentration advection: ' // errmsg2)
        end do
      else if (alloy_solidification) then
        call this%alloy_params%get('num-comp', n)
        allocate(this%cadv(n))
        do i = 1, n
          call this%cadv(i)%init(this%mesh, i, stat, errmsg2)
          if (stat /= 0) call TLS_fatal('DS_INIT: error initializing concentration advection: ' // errmsg2)
        end do
      end if
    end if

    call this%ds_params%get('cutvof', this%cutvof, default=0.0_r8)
    if (alloy_solidification) then
      call this%alloy_params%get('integrator', integrator)
    else
      call this%ds_params%get('integrator', integrator)
    end if

    if (alloy_solidification .and. this%have_void) &
        call TLS_fatal('DS_INIT: alloy solidification does not support presence of VOID')

    !! Figure out which diffusion solver we should be running, and ensure
    !! that the user has selected a compatible integration method.
    if (this%have_void .and. this%have_fluid_flow) then
      !! Transient void; use special solver.
      if (this%have_species_diffusion) then
        !! Only implemented for HT.
        INSIST(.false.)
      end if
      this%solver_type = SOLVER2
      if (integrator /= 'nonadaptive-bdf1') then
        call TLS_fatal ('DS_INIT: diffusion system characteristics are incompatible with STEPPING_METHOD choice.')
      end if
    else
      !! Void (if any) is fixed; use standard solver
      if (this%have_species_diffusion) then
        this%solver_type = SOLVER1
      else if (alloy_solidification) then
        this%solver_Type = SOLVER4
      else
        this%solver_Type = SOLVER3
      end if
      if (integrator /= 'adaptive-bdf2') then
        call TLS_fatal ('DS_INIT: diffusion system characteristics are incompatible with STEPPING_METHOD choice.')
      end if
    end if

    select case (this%solver_type)
    case (SOLVER1)
      block
        type(thermal_bc_factory1)    :: tbc_fac
        type(species_bc_factory1)    :: sbc_fac
        type(thermal_source_factory) :: tsrc_fac
        type(species_source_factory) :: ssrc_fac
        call tbc_fac%init(this%mesh, stefan_boltzmann, absolute_zero, this%bc_params)
        call sbc_fac%init(this%mesh, this%species_bc_params)
        call tsrc_fac%init(this%mesh, this%thermal_source_params)
        call ssrc_fac%init(this%mesh, this%species_source_params)
        this%mod1 => create_HTSD_model(tinit, this%disc, this%mmf, tbc_fac, sbc_fac, tsrc_fac, ssrc_fac, stat, errmsg)
      end block
      if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))
      this%sol1 => create_HTSD_solver (this%mmf, this%mod1, this%ds_params)

    case (SOLVER2)
      block
        type(thermal_bc_factory1)    :: bc_fac
        type(thermal_source_factory) :: tsrc_fac
        call bc_fac%init(this%mesh, stefan_boltzmann, absolute_zero, this%bc_params)
        call tsrc_fac%init(this%mesh, this%thermal_source_params)
        this%mod2 => create_FHT_model (tinit, this%disc, this%mmf, bc_fac, tsrc_fac, stat, errmsg)
      end block
      if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))
      this%sol2 => create_FHT_solver(this%mmf, this%mod2, this%ds_params)
      call move_alloc(this%hadv, this%sol2%hadv)

    case (SOLVER3)
      block
        type(thermal_bc_factory1)    :: tbc_fac
        type(thermal_source_factory) :: tsrc_fac
        call tbc_fac%init(this%mesh, stefan_boltzmann, absolute_zero, this%bc_params)
        call tsrc_fac%init(this%mesh, this%thermal_source_params)
        this%mod3 => create_ht_model(tinit, this%disc, this%mmf, tbc_fac, tsrc_fac, stat, errmsg2)
      end block
      if (stat /= 0) call TLS_fatal('DS_INIT: ' // trim(errmsg2))
      this%sol3 => create_ht_solver(this%mmf, this%mod3, this%ds_params, stat, errmsg2)
      if (stat /= 0) call TLS_fatal('DS_INIT: ' // trim(errmsg2))

    case (SOLVER4)
      block
        use material_model_driver, only: matl_model
        use material_class
        integer :: n
        class(material), pointer :: matl
        character(:), allocatable :: name
        call this%alloy_params%get('material', name, stat, errmsg2)
        if (stat /= 0) call TLS_fatal('DS_INIT: ' // errmsg2)
        n = matl_model%matl_index(name)
        if (n == 0) call TLS_fatal('DS_INIT: unknown material: ' // name)
        call matl_model%get_matl_ref(n, matl)
        allocate(this%mod4)
        call this%mod4%init(this%mesh, matl, this%alloy_params, stat, errmsg2)
        if (stat /= 0) call TLS_fatal('DS_INIT: ' // errmsg2)
        this%sol4 => create_alloy_solver(this%mod4, this%alloy_params, stat, errmsg2)
        if (stat /= 0) call TLS_fatal('DS_INIT: ' // trim(errmsg2))
      end block
    case default
      INSIST(.false.)
    end select

    call TLS_info ('  diffusion solver initialized')

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
    end function multiphase_problem

  end subroutine ds_init

  subroutine ds_delete ()
    if (associated(this%sol1)) then
      call HTSD_solver_delete (this%sol1)
      deallocate(this%sol1)
    end if
    if (associated(this%sol2)) then
      call FHT_solver_delete (this%sol2)
      deallocate(this%sol2)
    end if
    if (associated(this%mod1)) then
      call HTSD_model_delete (this%mod1)
      deallocate(this%mod1)
    end if
    if (associated(this%mod2)) then
      call FHT_model_delete (this%mod2)
      deallocate(this%mod2)
    end if
    if (associated(this%mmf)) deallocate(this%mmf)
    if (associated(this%disc)) deallocate(this%disc)
    this%mesh => null()
  end subroutine ds_delete

  subroutine ds_set_initial_state (t, dt, temp, conc)

    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)

    !! Permute the cell temperature array to the DS ordering.
    if (this%have_heat_transfer) then
      ASSERT(present(temp))
      ASSERT(size(temp) == this%mesh%ncell_onP)
    end if

    !! Permute the cell concentration array to the DS ordering.
    if (this%have_species_diffusion) then
      INSIST(present(conc))
      ASSERT(size(conc,dim=1) == this%mesh%ncell_onP)
      ASSERT(size(conc,dim=2) == num_species)
    end if

    !! Set the initial state in the appropriate solver.
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_set_initial_state (this%sol1, t, dt, temp, conc)
    case (SOLVER2)
      call FHT_solver_set_initial_state (this%sol2, t, temp)
    case (SOLVER3)
      call this%sol3%set_initial_state(t, temp, dt)
    case (SOLVER4)
      call this%sol4%set_initial_state(t, temp, dt)
    case default
      INSIST(.false.)
    end select

    if (this%solver_type == SOLVER4) then
      block
        use zone_module
        integer :: i
        INSIST(.not.allocated(zone%phi))
        allocate(zone%phi(this%sol4%num_comp, this%mesh%ncell_onp))
        do i = 1, size(zone%phi,dim=1)
          call this%sol4%get_C_liq(i, zone%phi(i,:))
        end do
      end block
    end if

  end subroutine ds_set_initial_state

  !! The effect of calling this subroutine is to restart or reset the solver so
  !! that its subsequent numerical behavior is as if it was starting integration
  !! from an initial state equal to the current state.  This mainly means
  !! dropping any previous solution history in the BDF2 solver and recomputing
  !! an approximation to the initial state time derivative.  The argument DT is
  !! used to compute that time derivative and is best chosen equal to the next
  !! time step size; however it has no bearing on the next step size used.

  subroutine ds_restart (dt)

    real(r8), intent(in) :: dt

    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_restart (this%sol1, dt)
    case (SOLVER2)
      ! nothing to do here yet
    case (SOLVER3)
      call this%sol3%restart(dt)
    case (SOLVER4)
      call this%sol4%restart(dt)
    case default
      INSIST(.false.)
    end select

  end subroutine ds_restart

!TODO: Replace this with something equivalent? The new material data type has no
! notion of being temperature dependent nor of being multi-component
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!
! !! VERIFY_MATERIAL_COMPATIBILITY
! !!
! !! This routine verifies that the attributes of the material systems
! !! referenced by the material mesh function MMF are compatible with the
! !! requirements of the type of diffusion system being solved.   If an
! !! incompatibility is detected STAT returns a nonzero value and an
! !! explanatory error message is returned in ERRMSG.
! !!
! !! The material system requirements are as follows.
! !!
! !! Species diffusion only:
! !!   * temperature-independent material systems
! !!   * number of species equal to the number of material components less 1
! !! Heat conduction only:
! !!   * temperature-dependent material systems
! !!   * single-component material systems
! !! Heat conduction/species diffusion:
! !!   * temperature-dependent material systems
! !!   * number of species equal to the number of material components less 1
! !!
!
!  subroutine verify_material_compatibility (mmf, stat, errmsg)
!
!    use material_system
!    use material_table
!
!    type(matl_mesh_func), intent(in) :: mmf
!    integer, intent(out) :: stat
!    character(len=*), intent(out) :: errmsg
!
!    integer :: i
!    integer, allocatable :: matid(:)
!    type(mat_system), pointer :: ms
!
!    !! Retrieve a list of all the material IDs that may be encountered.
!    call mmf%get_all_matl(matid, drop_void=.true.)
!
!    !! Verify that the material system attributes are compatible
!    !! with the constraints imposed by the type of diffusion system.
!    select case (ds_sys_type)
!    case (DS_SPEC_SYS)      ! Species diffusion
!      do i = 1, size(matid)
!        ms => mt_get_material(matid(i))
!        ASSERT(associated(ms))
!        if (ms_temp_dep(ms)) then
!          stat = -1
!          errmsg = 'diffusion system type requires temperature-independent material systems'
!          return
!        end if
!        if (ms_num_component(ms) /= num_species+1) then
!          stat = -1
!          errmsg = 'diffusion system type requires ' // i_to_c(num_species+1) // &
!                   '-component material systems'
!          return
!        end if
!      end do
!    case (DS_TEMP_SYS)      ! Heat conduction
!      do i = 1, size(matid)
!        ms => mt_get_material(matid(i))
!        ASSERT(associated(ms))
!        if (.not.ms_temp_dep(ms)) then
!          stat = -1
!          errmsg = 'diffusion system type requires temperature-dependent material systems'
!          return
!        end if
!        if (ms_num_component(ms) /= 1) then
!          stat = -1
!          errmsg = 'diffusion system type requires single-component material systems'
!          return
!        end if
!      end do
!    case (DS_TEMP_SPEC_SYS) ! Heat conduction and species diffusion
!      do i = 1, size(matid)
!        ms => mt_get_material(matid(i))
!        ASSERT(associated(ms))
!        if (.not.ms_temp_dep(ms)) then
!          stat = -1
!          errmsg = 'diffusion system type requires temperature-dependent material systems'
!          return
!        end if
!        if (ms_num_component(ms) /= num_species+1) then
!          stat = -1
!          errmsg = 'diffusion system type requires ' // i_to_c(num_species+1) // &
!                   '-component material systems'
!          return
!        end if
!      end do
!    case default
!      INSIST(.false.)
!    end select
!
!    stat = 0
!    errmsg = ''
!
!  end subroutine verify_material_compatibility

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
        call TLS_info ('DS: culled material fragments from ' // i_to_c(cell_count) // ' cells.')
    if (present(culled)) culled = (cell_count > 0)

  end subroutine cull_material_fragments

  subroutine update_moving_vf
    select case (this%solver_type)
    case (SOLVER1)
      call this%mod1%update_moving_vf
    case (SOLVER2)
      call this%mod2%update_moving_vf
    case (SOLVER3)
      call this%mod3%update_moving_vf
    case (SOLVER4)
      ! nothing
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine add_moving_vf_events(eventq)
    use sim_event_queue_type
    type(sim_event_queue), intent(inout) :: eventq
    select case (this%solver_type)
    case (SOLVER1)
      call this%mod1%add_moving_vf_events(eventq)
    case (SOLVER2)
      call this%mod2%add_moving_vf_events(eventq)
    case (SOLVER3)
      call this%mod3%add_moving_vf_events(eventq)
    case (SOLVER4)
      ! nothing
    case default
      INSIST(.false.)
    end select
  end subroutine

end module diffusion_solver
