!!
!! EM_HEAT_DRIVER
!!
!! Procedures for driving the computation of the EM heat source for EM
!! heating simulations.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Refactored February 2024, July 2025
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module em_heat_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use simpl_mesh_type
  use data_mapper_class
  use parameter_list_type
  use parallel_communication
  use truchas_logging_services
  use truchas_timers
  use induction_heat_solver_type
  use microwave_heat_solver_type
  implicit none
  private

  public :: em_heat_enabled, em_heat_driver_init, em_heat_driver_final
  public :: update_em_heat, em_heat_ptr, get_em_heat_event_times
  public :: read_em_heat_namelists, read_em_heat_restart_data, skip_em_heat_restart_data

  type :: em_heat_driver_data
    class(base_mesh), pointer :: ht_mesh => null() ! unowned reference
    type(simpl_mesh), pointer :: em_mesh => null() ! unowned reference
    class(data_mapper), allocatable :: ht2em
    logical :: const_eps=.false., const_mu=.false., const_sigma=.false.
    real(r8), allocatable :: eps(:), epsi(:), mu(:), sigma(:) ! on EM mesh
    real(r8), allocatable :: q_ht(:) ! heat source on HT mesh
    real(r8), allocatable :: q(:) ! computed heat source on EM mesh
    logical :: use_mw_solver
    type(induction_heat_solver), allocatable :: ih_solver
    type(microwave_heat_solver), allocatable :: mh_solver
  end type

  type(em_heat_driver_data), allocatable, target :: this
  type(parameter_list), pointer :: params => null()

contains

  logical function em_heat_enabled()
    use physics_module, only: em_heating
    em_heat_enabled = em_heating
  end function

  subroutine em_heat_driver_final
    if (allocated(this)) deallocate(this)
    if (associated(params)) deallocate(params)
  end subroutine

  function em_heat_ptr() result(ptr)
    real(r8), pointer :: ptr(:)
    ptr => this%q_ht
  end function

  !! This driver subroutine is called by the Truchas input driver and it
  !! coordinates the reading of all the namelists pertaining to EM heating,
  !! and assembles the input into the PARAMS parameter list held as a
  !! private module variable.

  subroutine read_em_heat_namelists(lun)
    use physics_module, only: induction_heating, microwave_heating
    use induction_heating_namelist
    use microwave_heating_namelist
    use induction_source_field_namelist
    use electromagnetic_bc_namelist
    use mesh_manager, only: enable_mesh
    integer, intent(in) :: lun
    logical :: exists
    type(parameter_list), pointer :: plist
    allocate(params)
    if (induction_heating) then
      call read_induction_heating_namelist(lun, params)
      plist => params%sublist('external-field')
      call read_induction_source_field_namelist(lun, plist)
    else if (microwave_heating) then
      call read_microwave_heating_namelist(lun, params)
    end if
    plist => params%sublist('bc')
    call read_electromagnetic_bc_namelists(lun, plist)
    call enable_mesh('em', exists)
    if (.not.exists) call TLS_fatal('EM_MESH namelist was not specified')
  end subroutine

  !! Return a list of times where the computed EM heat source changes abruptly.
  !! This info can be used to control time stepping to hit those times exactly
  !! and reduce the step size and restart the HT ODE integrator in response to
  !! the discontinuous source.

  subroutine get_em_heat_event_times(times)
    real(r8), allocatable, intent(out) :: times(:)
    if (this%use_mw_solver) then
      call this%mh_solver%get_event_times(times)
    else
      call this%ih_solver%get_event_times(times)
    end if
  end subroutine

  !! Initialize the EM heat driver. This is directly responsible for setting up
  !! the mapping between the HT and EM meshes and invoking the initialization
  !! of the user-specified type of EM heat solver (induction, microwave).

  subroutine em_heat_driver_init

    use physics_module, only: induction_heating, microwave_heating
    use mesh_manager, only: named_mesh_ptr, simpl_mesh_ptr
    use em_properties

    integer :: n, stat
    character(:), allocatable :: errmsg, data_mapper_kind

    ASSERT(.not.allocated(this))

    call start_timer('EM heat')

    call tls_info('')
    call tls_info('Initializing the EM heat solver ...')

    allocate(this)

    this%ht_mesh => named_mesh_ptr('main')
    this%em_mesh => simpl_mesh_ptr('em')

    !! Generate the mapping between the HT and EM meshes
    call params%get('data-mapper-kind', data_mapper_kind, stat, errmsg, default='default')
    if (stat /= 0) call TLS_fatal('EM_HEAT_DRIVER_INIT: ' // errmsg)
    select case (data_mapper_kind)
    case ('default', 'kuprat')
      block
        use kuprat_mapper_type
        call tls_info('creating kuprat mesh-to-mesh mapper ...')
        allocate(kuprat_mapper :: this%ht2em)
      end block
    case ('portage')
#ifdef USE_PORTAGE
      block
        use portage_mapper_type
        call tls_info('creating portage mesh-to-mesh mapper ...')
        allocate(portage_mapper :: this%ht2em)
      end block
#else
      call tls_fatal('EM_HEAT_DRIVER_INIT: this Truchas build does not support "portage" for data-mapper-kind')
#endif
    case default
      call tls_fatal('EM_HEAT_DRIVER_INIT: unknown data-mapper-kind: ' // data_mapper_kind)
    end select
    call this%ht2em%init(this%ht_mesh, this%em_mesh)

    n = this%em_mesh%ncell
    allocate(this%eps(n), this%epsi(n), this%mu(n), this%sigma(n))

    allocate(this%q_ht(this%ht_mesh%ncell_onp))
    allocate(this%q(this%em_mesh%ncell_onp))

    call define_default_em_properties

    block
      use physical_constants
      call params%set('vacuum-permittivity', vacuum_permittivity)
      call params%set('vacuum-permeability', vacuum_permeability)
    end block

    INSIST(induction_heating .neqv. microwave_heating)
    this%use_mw_solver = microwave_heating
    if (this%use_mw_solver) then
      allocate(this%mh_solver)
      call this%mh_solver%init(this%em_mesh, params, stat, errmsg)
      if (stat /= 0) call tls_fatal('EM_HEAT_DRIVER_INIT: error creating microwave solver: ' // errmsg)
    else
      allocate(this%ih_solver)
      call this%ih_solver%init(this%em_mesh, params, stat, errmsg)
      if (stat /= 0) call tls_fatal('EM_HEAT_DRIVER_INIT: error creating induction solver: ' // errmsg)
    end if

    call stop_timer('EM heat')

  end subroutine em_heat_driver_init

  !! Update the EM heat source array referenced by EM_HEAT_PTR() for the given
  !! time T and HT mesh cell temperatures TEMP. This is intended to be called
  !! every time step, and T should be the start of the step. The subroutine only
  !! updates the heat source when necessary, due to changes in the forcing of
  !! the EM system or significant changes to temperature-dependent EM material
  !! properties.

  subroutine update_em_heat(t, temp)

    real(r8), intent(in) :: t, temp(:)

    integer :: stat
    character(:), allocatable :: errmsg
    character(80) :: string

    ASSERT(allocated(this))

    call start_timer('EM heat')

    call set_em_properties(this, temp)
    if (this%use_mw_solver) then
      call this%mh_solver%update_em_heat(t, this%eps, this%epsi, this%mu, this%sigma, this%q, stat, errmsg)
      if (stat < 0) call tls_fatal('UPDATE_EM_HEAT: microwave solver error: ' // errmsg)
    else
      call this%ih_solver%update_em_heat(t, this%eps, this%epsi, this%mu, this%sigma, this%q, stat, errmsg)
      if (stat < 0) call tls_fatal('UPDATE_EM_HEAT: induction solver error: ' // errmsg)
    end if

    if (stat == 2) call write_em_heat_restart_data(t)

    if (stat > 0) then ! Q was updated
      call set_joule_power_density(this%q)
      write(string,fmt='(2(a,es11.4))') 'EM_HEAT_DRIVER: |Q|_max=', global_maxval(this%q), &
          ', Q_total=', global_dot_product(this%q, abs(this%em_mesh%volume(:size(this%q))))
      call tls_info(trim(string))
    end if

    call stop_timer('EM heat')

  contains

    subroutine set_joule_power_density(values)
      real(r8), intent(in) :: values(:)
      !call start_timer('mesh-to-mesh mapping')
      call this%ht2em%map_field(values(:this%em_mesh%ncell_onP), this%q_ht, defval=0.0_r8, &
                           map_type=GLOBALLY_CONSERVATIVE, pullback=.true.)
      !call stop_timer('mesh-to-mesh mapping')
    end subroutine

  end subroutine update_em_heat

  !! This auxiliary subroutine sets the values of the EPS, MU, and SIGMA
  !! property arrays. The property values are computed on the heat transfer
  !! mesh using the current value of temperature, and mapped to the EM mesh.
  !!
  !! NB: At the outset, the CONST_* flags should be false, so that the values
  !! are set, but if a property is constant (not varying with temperature)
  !! subsequent calls will skip computing and setting the value.

  subroutine set_em_properties(this, temp)
    use em_properties
    class(em_heat_driver_data), intent(inout) :: this
    real(r8), intent(in) :: temp(:)
    real(r8), allocatable :: value(:) !TODO persistent workspace
    allocate(value(this%ht_mesh%ncell_onP))
    if (.not.this%const_eps) then
      call get_permittivity(temp, value)
      call set_em_mesh_prop(value, this%eps, defval=1.0_r8)
      call get_permittivity_im(temp, value)
      call set_em_mesh_prop(value, this%epsi, defval=0.0_r8)
      this%const_eps = permittivity_is_const() .and. permittivity_im_is_const()
    end if
    if (.not.this%const_mu) then
      call get_permeability(temp, value)
      call set_em_mesh_prop(value, this%mu, defval=1.0_r8)
      this%const_mu = permeability_is_const()
    end if
    if (.not.this%const_sigma) then
      call get_conductivity(temp, value)
      call set_em_mesh_prop(value, this%sigma, defval=0.0_r8)
      this%const_sigma = conductivity_is_const()
    end if

  contains

    subroutine set_em_mesh_prop(ht_mesh_prop, em_mesh_prop, defval)
      real(r8), intent(in)  :: ht_mesh_prop(:), defval
      real(r8), intent(out) :: em_mesh_prop(:)
      call this%ht2em%map_field(ht_mesh_prop, em_mesh_prop(:this%em_mesh%ncell_onP), defval=defval, map_type=LOCALLY_BOUNDED)
      call this%em_mesh%cell_imap%gather_offp(em_mesh_prop)
    end subroutine

  end subroutine set_em_properties

  !! Write the EM heat data and inputs used to compute it to the HDF output
  !! file. This is checkpoint data used to generate restart files. Output for
  !! visualization is written elsewhere.

  subroutine write_em_heat_restart_data(t)

    use truchas_h5_outfile, only: th5_sim_group
    use truchas_danu_output_data, only: outfile, io_group_size
    use truchas_env, only: output_file_name

    real(r8), intent(in) :: t

    integer, save :: sim_num = 0
    character(32) :: sim_name
    type(th5_sim_group) :: sim

    ASSERT(allocated(this))

    call outfile%reopen(output_file_name('h5'), io_group_size, is_IOP)

    !! Write the data.
    sim_num = sim_num + 1
    write(sim_name,'(a,i3.3)') 'EM', sim_num
    call outfile%add_sim_group(trim(sim_name), sim)
    call sim%write_attr('TIME', t)
    call TLS_info('EM_HEAT_DRIVER: writing EM restart data for ' // trim(sim_name))

    if (this%use_mw_solver) then
      call sim%write_attr('EM-KIND', 2) ! microwave heating
      call this%mh_solver%write_restart_data(sim)
    else
      call sim%write_attr('EM-KIND', 1) ! induction heating
      call this%ih_solver%write_restart_data(sim)
    end if

    call outfile%close

  end subroutine write_em_heat_restart_data

  !! Read the EM heat segment of the restart file opened (and prepositioned)
  !! on UNIT, and initialize the Q_HT component and parameters specific to
  !! the kind of EM heat simulation. If any errors or incompatibilities with
  !! data read from the input file are encountered, execution is halted.

  subroutine read_em_heat_restart_data(unit, version)

    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c
    use parallel_communication, only: global_sum

    integer, intent(in) :: unit,  version

    integer :: n

    ASSERT(allocated(this))

    call tls_info('Reading the EM heat data from the restart file')

    call read_var(unit, n, 'READ_EM_HEAT_RESTART_DATA: error reading EM_HEAT_TYPE record')
    select case (n)
    case (1) ! induction heating
      if (this%use_mw_solver) then
        call tls_fatal('READ_EM_HEAT_RESTART_DATA: incompatible solver data')
      else
        call this%ih_solver%read_restart_data(unit, version)
      end if
    case (2) ! microwave heating
      if (this%use_mw_solver) then
        call this%mh_solver%read_restart_data(unit, version)
      else
        call tls_fatal('READ_EM_HEAT_RESTART_DATA: incompatible solver data')
      endif
    case default
      INSIST(.false.)
    end select

  end subroutine read_em_heat_restart_data

  !! Advance over the EM heat segment of the restart file opened (and
  !! prepositioned) on UNIT. This is used when the restart file contains
  !! EM heat data but em heating is not enabled.

  subroutine skip_em_heat_restart_data(unit, version)
    use restart_utilities, only: skip_records, read_var
    integer, intent(in) :: unit, version
    integer :: n
    call read_var(unit, n, 'SKIP_EM_HEAT_RESTART_DATA: error reading EM_HEAT_TYPE record')
    select case (n)
    !NB: skip_restart_data is nopass so can be invoked when the variable is not allocated!
    case (1)
      call this%ih_solver%skip_restart_data(unit, version)
    case (2)
      call this%mh_solver%skip_restart_data(unit, version)
    end select
  end subroutine

end module em_heat_driver
