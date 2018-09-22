!!
!! DIFFUSION_SOLVER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Markus Berndt <berndt@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module diffusion_solver

  use kinds
  use diffusion_solver_data
  use mesh_manager
  use parallel_communication
  use truchas_logging_services
  use mfd_disc_type
  use matl_mesh_func_type
  use string_utilities, only: i_to_c
  use source_mesh_function
  use mesh_interop
  use FHT_model_type
  use FHT_solver_type
  use HTSD_model_type
  use HTSD_solver_type
  use truchas_timers
  use truchas_logging_services
  use unstr_mesh_type
  implicit none
  private

  public :: ds_init, ds_step
  public :: read_ds_namelists
  public :: ds_set_initial_state, ds_restart
  public :: ds_delete

  !! These return cell-centered results relative to the old Truchas mesh.
  public :: ds_get_temp, ds_get_enthalpy, ds_get_phi
  public :: ds_get_temp_grad

  !! These return results relative to the new distributed mesh.
  public :: ds_get_cell_temp, ds_get_face_temp

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
    !! Saved references to the model sources.
    type(source_mf), pointer :: ht_source => null()
    type(source_mf), pointer :: sd_source(:) => null()
    !! Solver selections
    integer :: solver_type = 0
    !! The standard model and solver.
    type(HTSD_model),  pointer :: mod1 => null()
    type(HTSD_solver), pointer :: sol1 => null()
    !! The special model and solver that works with transient void.
    type(FHT_model),  pointer :: mod2  => null()
    type(FHT_solver), pointer :: sol2 => null()
  end type ds_driver
  type(ds_driver), save :: this
  
  integer, parameter :: SOLVER1 = 1 ! the standard solver
  integer, parameter :: SOLVER2 = 2 ! special solver that works with transient void
  
contains

  subroutine read_ds_namelists (lun)
  
    use ds_boundary_condition_input
    use ds_interface_condition_input
    use ds_source_input, only: read_ds_source
    use ER_input
  
    integer, intent(in) :: lun
    
    call read_ds_namelist (lun)
    call read_ds_boundary_condition (lun)
    call read_ds_interface_condition (lun)
    call read_ds_source (lun)
    
    call ERI_read_enclosure_radiation (lun)
    call ERI_read_enclosure_surface (lun)
    
  end subroutine read_ds_namelists

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DS_STEP
 !!

  subroutine ds_step (h, hnext, errc)

    use zone_module, only: zone
    use flow_phase_change, only: set_reference_fluid_density, set_solidified_density
    use cutoffs_module, only: cutvof

    real(r8), intent(inout) :: h
    real(r8), intent(out)   :: hnext
    integer,  intent(out)   :: errc

    real(r8) :: hlast, t
    integer :: counters(6)
    character(len=80) :: message(2)
    real(r8), pointer :: state(:,:)
    
    call start_timer ('Diffusion Solver')

    if (this%have_heat_transfer .and. this%have_phase_change) call set_reference_fluid_density
    
    if (this%have_fluid_flow) then
      call update_mmf_from_matl (this%mmf)
      if (this%have_void) call cull_material_fragments (this%mmf, cutvof)
    end if
    
    if (this%have_heat_transfer) call update_adv_heat
    if (this%have_species_diffusion) call update_adv_conc
    
    select case (this%solver_type)
    case (SOLVER1)  ! HT/SD solver with static void
      t = h + HTSD_solver_last_time(this%sol1)
      call HTSD_solver_advance_state (this%sol1, h, hnext, errc)
      if (errc == 0) call HTSD_solver_commit_pending_state (this%sol1)
    case (SOLVER2)  ! HT solver with transient void
      t = h + FHT_solver_last_time(this%sol2)
      call FHT_solver_advance_state (this%sol2, t, errc)
      if (errc == 0) call FHT_solver_commit_pending_state (this%sol2)
      hnext = huge(1.0_r8)  ! cede time step control to someone else
    case default
      INSIST(.false.)
    end select
    if (errc /= 0) return
    
    !! Update MATL in contexts that can modify the phase distribution.
    if (this%have_phase_change) then
      call create_state_array (state)
      call update_matl_from_mmf (this%mmf, state)
      deallocate(state)
    end if

    !! Update Truchas data structures to reflect the new HT/SD solution.
    if (this%have_heat_transfer) then
      call ds_get_temp (zone%temp)
      call ds_get_enthalpy (zone%enthalpy)
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
    case default
      INSIST(.false.)
    end select
    1 format(/,'DS: dt=',es9.3,', NFUN:NPC=',i7.7,':',i5.5,', NNR:NNF:NSR=',3(i4.4,:,':'))
    2 format(/,'DS: dt=',es9.3,', NFUN:NPC:NPA=',3(i7.7,:,':'))
    call TLS_info (message)
    
    if (this%have_heat_transfer .and. this%have_phase_change) call set_solidified_density
    
    call stop_timer ('Diffusion Solver')

  contains
  
    subroutine create_state_array (state)
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
          call HTSD_solver_get_cell_conc_copy (this%sol1, i, state(:,i))
        end do
      case (SOLVER2)
        call FHT_solver_get_cell_temp_copy (this%sol2, state(:,0))
      case default
        INSIST(.false.)
      end select
    end subroutine create_state_array
    
    subroutine update_adv_heat

      use legacy_mesh_api, only: ncells
      use advection_module,   only: compute_advected_enthalpy
      use EM_data_proxy,      only: joule_power_density
      use index_partitioning, only: gather_boundary

      real(r8), allocatable :: q_t(:), q_ds(:), dQ_t(:), dQ_ds(:), T_t(:)
      
      if (this%have_joule_heat .or. (this%have_fluid_flow .and. this%solver_type /= SOLVER2)) then
        allocate(q_ds(this%mesh%ncell))
        !! Joule heat source.
        if (this%have_joule_heat) then
          allocate(q_t(ncells))
          q_t = joule_power_density()
          q_ds(:this%mesh%ncell_onP) = q_t(:this%mesh%ncell_onP)
          call gather_boundary (this%mesh%cell_ip, q_ds)
          deallocate(q_t)
        else
          q_ds = 0.0_r8
        end if
        !! Advected heat source.
        if (this%have_fluid_flow .and. this%solver_type /= SOLVER2) then
          allocate(T_t(ncells), dQ_t(ncells))
          call ds_get_temp (T_t)
          call compute_advected_enthalpy (T_t, dQ_t)
          deallocate(T_t)
          allocate(dQ_ds(this%mesh%ncell))
          dQ_ds(:this%mesh%ncell_onP) = dQ_t(:this%mesh%ncell_onP)
          call gather_boundary (this%mesh%cell_ip, dQ_ds)
          q_ds = q_ds + (dQ_ds / (h * this%mesh%volume))
          deallocate(dQ_ds, dQ_t)
        end if
        call smf_set_extra_source (this%ht_source, q_ds)
        deallocate(q_ds)
      end if

    end subroutine update_adv_heat
    
    subroutine update_adv_conc

      use legacy_mesh_api, only: ncells
      use advection_module,   only: advected_phi
      use index_partitioning, only: gather_boundary

      integer :: i
      real(r8), allocatable :: q_t(:), q_ds(:), phi_t(:)
      
      if (this%have_fluid_flow) then
        allocate(q_t(ncells), q_ds(this%mesh%ncell), phi_t(ncells))
        do i = 1, num_species
          call ds_get_phi (i, phi_t)
          call advected_phi(phi_t, q_t)    ! species deltas for the time step
          q_ds(:this%mesh%ncell_onP) = q_t(:this%mesh%ncell_onP)
          call gather_boundary (this%mesh%cell_ip, q_ds)
          q_ds = q_ds / (h * this%mesh%volume) ! convert to a rate density
          call smf_set_extra_source (this%sd_source(i), q_ds)
        end do
        deallocate(q_t, q_ds, phi_t)
      end if

    end subroutine update_adv_conc

  end subroutine ds_step
  
  subroutine ds_get_temp (array)
    use legacy_mesh_api, only: ncells
    real(r8), intent(inout) :: array(:)
    ASSERT(this%have_heat_transfer)
    ASSERT(size(array) == ncells)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_temp_copy (this%sol1, array)
    case (SOLVER2)
      call FHT_solver_get_cell_temp_copy (this%sol2, array)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_enthalpy (array)
    use legacy_mesh_api, only: ncells
    real(r8), intent(inout) :: array(:)
    ASSERT(this%have_heat_transfer)
    ASSERT(size(array) == ncells)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_heat_copy (this%sol1, array)
    case (SOLVER2)
      call FHT_solver_get_cell_heat_copy (this%sol2, array)
    case default
      INSIST(.false.)
    end select
  end subroutine

  subroutine ds_get_phi (n, array)
    use legacy_mesh_api, only: ncells
    integer,  intent(in)  :: n
    real(r8), intent(inout) :: array(:)
    ASSERT(this%have_species_diffusion)
    ASSERT(size(array) == ncells)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_conc_copy (this%sol1, n, array)
    case default
      INSIST(.false.)
    end select
  end subroutine ds_get_phi

  subroutine ds_get_temp_grad (array)
    use legacy_mesh_api, only: ncells
    real(r8), intent(inout) :: array(:,:)
    ASSERT(size(array,dim=2) == ncells)
    ASSERT(size(array,dim=1) == 3)
    ASSERT(this%have_heat_transfer)
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_get_cell_temp_grad (this%sol1, array(:,:this%mesh%ncell_onP))
    case (SOLVER2)
      call FHT_solver_get_cell_temp_grad (this%sol2, array(:,:this%mesh%ncell_onP))
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
    case default
      INSIST(.false.)
    end select
  end subroutine ds_get_face_temp
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DS_INIT
 !!
  subroutine ds_init ()

    use EM_data_proxy, only: EM_is_on
    use fluid_data_module,  only: fluid_flow
    use FHT_model_factory
    use FHT_solver_factory
    use HTSD_model_factory
    use HTSD_solver_factory

    integer :: stat
    character(len=200) :: errmsg
    character(:), allocatable :: errmsg2

    call TLS_info ('')
    call TLS_info ('Initializing diffusion solver ...')
    
    !! Common initialization.
    this%mesh => unstr_mesh_ptr(mesh_name)
    INSIST(associated(this%mesh))
    
    allocate(this%disc)
    call this%disc%init (this%mesh, use_new_mfd)
    
    allocate(this%mmf)
    call mmf_init (this%mesh, this%mmf, stat, errmsg2)
    if (stat /= 0) call TLS_fatal ('DS_INIT: ' // errmsg2)
    
    call verify_material_compatibility (this%mmf, stat, errmsg)
    if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))
    
    !! Problem attributes
    this%have_heat_transfer = heat_eqn
    this%have_species_diffusion = (num_species > 0)
    this%have_joule_heat = EM_is_on()
    this%have_fluid_flow = fluid_flow
    this%have_phase_change = multiphase_problem(this%mmf)
    this%have_void = void_is_present()
    
    !! Figure out which diffusion solver we should be running, and ensure
    !! that the user has selected a compatible integration method.
    if (this%have_void .and. this%have_fluid_flow) then
      !! Transient void; use special solver.
      if (this%have_species_diffusion) then
        !! Only implemented for HT.
        INSIST(.false.)
      end if
      this%solver_type = SOLVER2
      if (integrator /= DS_NONADAPTIVE_BDF1) then
        call TLS_fatal ('DS_INIT: diffusion system characteristics are incompatible with STEPPING_METHOD choice.')
      end if
    else
      !! Void (if any) is fixed; use standard solver
      this%solver_type = SOLVER1
      if (integrator /= DS_ADAPTIVE_BDF2) then
        call TLS_fatal ('DS_INIT: diffusion system characteristics are incompatible with STEPPING_METHOD choice.')
      end if
    end if
    
    select case (this%solver_type)
    case (SOLVER1)
      this%mod1 => create_HTSD_model (this%disc, this%mmf, stat, errmsg)
      if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))
      if (this%have_heat_transfer) this%ht_source => this%mod1%ht%source
      if (this%have_species_diffusion) this%sd_source => this%mod1%sd%source
      this%sol1 => create_HTSD_solver (this%mmf, this%mod1, stat, errmsg)
      if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))
      
    case (SOLVER2)
      this%mod2 => create_FHT_model (this%disc, this%mmf, stat, errmsg)
      if (stat /= 0) call TLS_fatal ('DS_INIT: ' // trim(errmsg))
      this%ht_source => this%mod2%q ! we need this to set the advected heat at each step
      this%sol2 => create_FHT_solver(this%mmf, this%mod2, stat, errmsg)
      
    case default
      INSIST(.false.)
    end select
    
    call TLS_info ('  Diffusion solver initialized.')

  contains
  
    logical function multiphase_problem (mmf)
      use material_system
      use material_table
      type(matl_mesh_func), intent(in) :: mmf
      integer, allocatable :: list(:)
      type(mat_system), pointer :: ms
      integer :: i
      call mmf%get_all_matl(list, drop_void=.true.)
      do i = size(list), 1, -1
        ms => mt_get_material(list(i))
        ASSERT(associated(ms))
        if (ms_num_phase(ms) > 1) exit
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
    this%ht_source => null()
    this%sd_source => null()
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
  
    use legacy_mesh_api, only: ncells

    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional, target :: temp(:), conc(:,:)
    
    real(r8), pointer :: temp_ds(:) => null(), conc_ds(:,:) => null()

    !! Permute the cell temperature array to the DS ordering.
    if (this%have_heat_transfer) then
      ASSERT(present(temp))
      ASSERT(size(temp) == ncells)
      temp_ds => temp(:this%mesh%ncell_onP)
    end if
    
    !! Permute the cell concentration array to the DS ordering.
    if (this%have_species_diffusion) then
      INSIST(present(conc))
      ASSERT(size(conc,dim=1) == ncells)
      ASSERT(size(conc,dim=2) == num_species)
      conc_ds => conc(:this%mesh%ncell_onP,:)
    end if
    
    !! Set the initial state in the appropriate solver.
    select case (this%solver_type)
    case (SOLVER1)
      call HTSD_solver_set_initial_state (this%sol1, t, temp_ds, conc_ds, dt)
    case (SOLVER2)
      call FHT_solver_set_initial_state (this%sol2, t, temp_ds)
    case default
      INSIST(.false.)
    end select
    
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
    case default
      INSIST(.false.)
    end select

  end subroutine ds_restart

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VERIFY_MATERIAL_COMPATIBILITY
 !!
 !! This routine verifies that the attributes of the material systems
 !! referenced by the material mesh function MMF are compatible with the
 !! requirements of the type of diffusion system being solved.   If an
 !! incompatibility is detected STAT returns a nonzero value and an
 !! explanatory error message is returned in ERRMSG.
 !!
 !! The material system requirements are as follows.
 !!
 !! Species diffusion only:
 !!   * temperature-independent material systems
 !!   * number of species equal to the number of material components less 1
 !! Heat conduction only:
 !!   * temperature-dependent material systems
 !!   * single-component material systems
 !! Heat conduction/species diffusion:
 !!   * temperature-dependent material systems
 !!   * number of species equal to the number of material components less 1
 !!

  subroutine verify_material_compatibility (mmf, stat, errmsg)

    use material_system
    use material_table

    type(matl_mesh_func), intent(in) :: mmf
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i
    integer, allocatable :: matid(:)
    type(mat_system), pointer :: ms

    !! Retrieve a list of all the material IDs that may be encountered.
    call mmf%get_all_matl(matid, drop_void=.true.)

    !! Verify that the material system attributes are compatible
    !! with the constraints imposed by the type of diffusion system.
    select case (ds_sys_type)
    case (DS_SPEC_SYS)      ! Species diffusion
      do i = 1, size(matid)
        ms => mt_get_material(matid(i))
        ASSERT(associated(ms))
        if (ms_temp_dep(ms)) then
          stat = -1
          errmsg = 'diffusion system type requires temperature-independent material systems'
          return
        end if
        if (ms_num_component(ms) /= num_species+1) then
          stat = -1
          errmsg = 'diffusion system type requires ' // i_to_c(num_species+1) // &
                   '-component material systems'
          return
        end if
      end do
    case (DS_TEMP_SYS)      ! Heat conduction
      do i = 1, size(matid)
        ms => mt_get_material(matid(i))
        ASSERT(associated(ms))
        if (.not.ms_temp_dep(ms)) then
          stat = -1
          errmsg = 'diffusion system type requires temperature-dependent material systems'
          return
        end if
        if (ms_num_component(ms) /= 1) then
          stat = -1
          errmsg = 'diffusion system type requires single-component material systems'
          return
        end if
      end do
    case (DS_TEMP_SPEC_SYS) ! Heat conduction and species diffusion
      do i = 1, size(matid)
        ms => mt_get_material(matid(i))
        ASSERT(associated(ms))
        if (.not.ms_temp_dep(ms)) then
          stat = -1
          errmsg = 'diffusion system type requires temperature-dependent material systems'
          return
        end if
        if (ms_num_component(ms) /= num_species+1) then
          stat = -1
          errmsg = 'diffusion system type requires ' // i_to_c(num_species+1) // &
                   '-component material systems'
          return
        end if
      end do
    case default
      INSIST(.false.)
    end select

    stat = 0
    errmsg = ''

  end subroutine verify_material_compatibility

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

end module diffusion_solver
