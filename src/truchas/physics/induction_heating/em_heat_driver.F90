!!
!! EM_HEAT_DRIVER
!!
!! Procedures for driving the computation of the EM heat source for EM
!! heating simulations.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Refactored February 2024
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
  implicit none
  private

  public :: em_heat_enabled, em_heat_driver_init, em_heat_driver_final
  public :: update_em_heat, em_heat_ptr
  public :: read_em_heat_namelists, read_em_heat_restart_data, skip_em_heat_restart_data

  type :: em_heat_driver_data
    class(base_mesh), pointer :: ht_mesh => null() ! unowned reference
    type(simpl_mesh), pointer :: em_mesh => null() ! unowned reference
    class(data_mapper), allocatable :: ht2em
    ! EM properties
    logical :: const_eps=.false., const_mu=.false., const_sigma=.false.
    real(r8), allocatable :: eps(:), epsi(:), mu(:), sigma(:)  ! on EM mesh
    real(r8), allocatable :: q_ht(:)
    ! Computed EM heat and the variable inputs used
    real(r8), allocatable :: q(:)
    logical :: q_written = .false.
    type(induction_heat_solver), allocatable :: ih_solver
  end type

  type(em_heat_driver_data), allocatable, target :: this
  type(parameter_list), pointer :: params

contains

  logical function em_heat_enabled()
    use physics_module, only: electromagnetics
    em_heat_enabled = electromagnetics
  end function

  subroutine em_heat_driver_final
    if (allocated(this)) deallocate(this)
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
    use electromagnetics_namelist
    use induction_source_field_namelist
    use electromagnetic_bc_namelist
    use mesh_manager, only: enable_mesh
    integer, intent(in) :: lun
    logical :: exists
    type(parameter_list), pointer :: plist
    allocate(params)
    call read_electromagnetics_namelist(lun, params)
    plist => params%sublist('external-field')
    call read_induction_source_field_namelist(lun, plist)
    plist => params%sublist('bc')
    call read_electromagnetic_bc_namelists(lun, plist)
    call enable_mesh('em', exists)
    if (.not.exists) call TLS_fatal('EM_MESH namelist was not specified')
  end subroutine

  subroutine em_heat_driver_init

    use mesh_manager, only: named_mesh_ptr, simpl_mesh_ptr
    use em_properties

    integer :: n, stat
    character(:), allocatable :: errmsg, data_mapper_kind
    type(parameter_list), pointer :: plist

    ASSERT(.not.allocated(this))

    call TLS_info('')
    call TLS_info('Initializing electromagnetic heat solver ...')

    allocate(this)

    this%ht_mesh => named_mesh_ptr('main')
    this%em_mesh => simpl_mesh_ptr('em')

    !! Generate the mapping between the HT and EM meshes
    call params%get('data-mapper-kind', data_mapper_kind, stat, errmsg, default='default')
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    select case (data_mapper_kind)
    case ('default', 'kuprat')
      block
        use kuprat_mapper_type
        call TLS_info('creating kuprat mesh-to-mesh mapper ...')
        allocate(kuprat_mapper :: this%ht2em)
      end block
    case ('portage')
#ifdef USE_PORTAGE
      block
        use portage_mapper_type
        call TLS_info('creating portage mesh-to-mesh mapper ...')
        allocate(portage_mapper :: this%ht2em)
      end block
#else
      call TLS_fatal('IH_DRIVER_INIT: this Truchas build does not support "portage" for data-mapper-kind')
#endif
    case default
      call TLS_fatal('IH_DRIVER_INIT: unknown data-mapper-kind: ' // data_mapper_kind)
    end select
    call this%ht2em%init(this%ht_mesh, this%em_mesh)

    n = this%em_mesh%ncell
    allocate(this%eps(n), this%epsi(n), this%mu(n), this%sigma(n))

    allocate(this%q_ht(this%ht_mesh%ncell_onP))
    !TODO: NEED TO ADD Q_EPSI AND INCLUDE WITH RESTART DATA
    allocate(this%q(this%em_mesh%ncell_onP))
    !NB: the allocation status of Q_DATA indicates whether Q contains data

    call define_default_em_properties

    allocate(this%ih_solver)
    call this%ih_solver%init(this%em_mesh, params)

  end subroutine em_heat_driver_init


  subroutine update_em_heat(t, temp)

    real(r8), intent(in) :: t, temp(:)

    integer :: stat

    ASSERT(allocated(this))

    call start_timer('EM heat')

    call set_em_properties(this, temp)
    call this%ih_solver%update_em_heat(t, this%eps, this%epsi, this%mu, this%sigma, this%q, stat)
    ! OR call update_dielectric_heat
    if (stat > 0) then ! Q was updated
      call set_joule_power_density(this%q)
      this%q_written = .false.
    end if

    if (.not.this%q_written) then
      call write_em_heat_restart_data(t)
      this%q_written = .true.
    end if

    call stop_timer('EM heat')

  contains

    subroutine set_joule_power_density(values)
      real(r8), intent(in) :: values(:)
      call start_timer('mesh-to-mesh mapping')
      call this%ht2em%map_field(values(:this%em_mesh%ncell_onP), this%q_ht, defval=0.0_r8, &
                           map_type=GLOBALLY_CONSERVATIVE, pullback=.true.)
      call stop_timer('mesh-to-mesh mapping')
    end subroutine

  end subroutine update_em_heat

  !! These auxiliary subroutines set the values of the EPS, MU, and SIGMA
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
      call set_permittivity(this, value)
      call get_permittivity_im(temp, value)
      call set_permittivity_im(this, value)
      this%const_eps = permittivity_is_const() .and. permittivity_im_is_const()
    end if
    if (.not.this%const_mu) then
      call get_permeability(temp, value)
      call set_permeability(this, value)
      this%const_mu = permeability_is_const()
    end if
    if (.not.this%const_sigma) then
      call get_conductivity(temp, value)
      call set_conductivity(this, value)
      this%const_sigma = conductivity_is_const()
    end if
  end subroutine

  subroutine set_permittivity(this, values)
    class(em_heat_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%eps(:this%em_mesh%ncell_onP), defval=1.0_r8, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%eps)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

  subroutine set_permittivity_im(this, values)
    class(em_heat_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%epsi(:this%em_mesh%ncell_onP), defval=0.0_r8, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%epsi)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

  subroutine set_permeability(this, values)
    class(em_heat_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%mu(:this%em_mesh%ncell_onP), defval=1.0_r8, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%mu)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

  subroutine set_conductivity(this, values)
    class(em_heat_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%sigma(:this%em_mesh%ncell_onP), defval=0.0_r8, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%sigma)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

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

    call TLS_info ('  reading the EM heat data from the restart file')

    call read_var(unit, n, 'READ_EM_DATA: error reading EM_HEAT_TYPE record')
    select case (n)
    case (1) ! induction heating
      call this%ih_solver%read_restart_data(unit, version)
    case (2) ! microwave heating
      !call read_mh_data(unit, version)
    case default
    end select

    !TODO: joule size is only on-process cells?
    call read_var(unit, n, 'READ_EM_DATA: error reading NHEAT record')
    if (n /= global_sum(size(this%q_ht))) &
        call halt('READ_EM_DATA: incompatible NHEAT value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_ht, &
        this%ht_mesh%xcell(:this%ht_mesh%ncell_onp), &
        'READ_EM_DATA: error reading HEAT record')
    this%q_written = .false.

  end subroutine read_em_heat_restart_data

  !! Advance over the Joule heat segment of the restart file opened (and
  !! prepositioned) on UNIT. This is used when the restart file contains Joule
  !! heat data but induction heating is not enabled. NB: the number of records
  !! skipped must match the number of records read by READ_JOULE_DATA.

  subroutine skip_em_heat_restart_data(unit, version)
    use restart_utilities, only: skip_records, read_var
    integer, intent(in) :: unit, version
    integer :: n
    call read_var(unit, n, 'SKIP_EM_DATA: error reading EM_HEAT_TYPE record')
    select case (n)
    case (1)
      call this%ih_solver%skip_restart_data(unit, version)
    !case (2)
    !  call skip_mh_data(unit, version)
    end select
    call skip_records(unit, 2, 'SKIP_EM_DATA: error skipping the EM heat data')
  end subroutine

  !! Write the Joule heat data and inputs used to compute it to the HDF output
  !! file. This is checkpoint data used to generate restart files. Output for
  !! visualization is written elsewhere.

  subroutine write_em_heat_restart_data(t)

    use parallel_communication
    use permutations
    use truchas_h5_outfile, only: th5_sim_group
    use truchas_danu_output_data, only: outfile, io_group_size
    use truchas_env, only: output_file_name

    real(r8), intent(in) :: t

    integer :: n
    integer, allocatable :: cell_perm(:)
    real(r8), allocatable :: col_mu(:), col_sigma(:)

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
    call TLS_info('  writing EM restart data for ' // trim(sim_name))

    call sim%write_attr('EM-KIND', 1) ! induction heating
    call this%ih_solver%write_restart_data(sim)

    call sim%write_dist_array('HEAT', this%q_ht, global_sum(size(this%q_ht)))

    call outfile%close

  end subroutine write_em_heat_restart_data

end module em_heat_driver
