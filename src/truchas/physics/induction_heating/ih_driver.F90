!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ih_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use simpl_mesh_type
  use data_mapper_class
  use parallel_communication
  use ih_source_factory_type
  use truchas_logging_services
  use truchas_timers
  use scalar_func_class
  implicit none
  private

  public :: ih_enabled, ih_driver_init, update_joule_heat
  public :: joule_power_density
  public :: set_EM_simulation_on_or_off, EM_is_on
  public :: read_joule_data, skip_joule_data

  !! EM simulation state variable
  logical, save :: EM_enabled = .false.

  type :: ih_driver_data
    class(base_mesh), pointer :: ht_mesh => null() ! unowned reference
    type(simpl_mesh), pointer :: em_mesh => null() ! unowned reference
    class(data_mapper), allocatable :: ht2em
    type(ih_source_factory) :: src_fac
    character(32) :: coil_md5sum
    logical :: const_eps=.false., const_mu=.false., const_sigma=.false.
    real(r8) :: matl_change_threshold
    real(r8), allocatable :: eps(:), mu(:), sigma(:)  ! on EM mesh
    ! Joule heat and the varying inputs used to compute it
    real(r8), allocatable :: q(:) ! on HT mesh
    real(r8), allocatable :: q_data(:), q_eps(:), q_mu(:), q_sigma(:)
    logical :: q_written = .false.
    logical :: q_restart = .false. ! true when q was read from a restart file
  end type
  type(ih_driver_data), allocatable, target :: this

contains

  subroutine ih_driver_init(t)

    use mesh_manager, only: named_mesh_ptr, simpl_mesh_ptr
    use electromagnetics_namelist, only: params
    use EM_properties

    real(r8), intent(in) :: t

    integer :: n, stat
    character(:), allocatable :: errmsg, data_mapper_kind

    !TODO: barrier if EM is not enabled?

    call TLS_info('')
    call TLS_info('Initializing induction heating solver ...')

    allocate(this)

    this%ht_mesh => named_mesh_ptr('main')
    this%em_mesh => simpl_mesh_ptr('em')

    !! Generate the mapping between the HT and EM meshes
    call params%get('data-mapper-kind', data_mapper_kind, stat=stat, errmsg=errmsg, default='default')
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
    allocate(this%eps(n), this%mu(n), this%sigma(n))

    allocate(this%q(this%ht_mesh%ncell_onP))
    allocate(this%q_eps(n), this%q_mu(n), this%q_sigma(n))
    !NB: the allocation status of Q_DATA indicates whether Q contains data

    !! Initialize the external H-field source factory
    call this%src_fac%init(params, t, stat, errmsg)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    this%coil_md5sum = this%src_fac%coil_geom_fingerprint()

    call define_default_em_properties

    call params%get('matl-change-threshold', this%matl_change_threshold, stat=stat, errmsg=errmsg, default=0.3_r8)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    if (this%matl_change_threshold <= 0.0_r8) call TLS_fatal('IH_DRIVER_INIT: matl-change-threshold must be > 0.0')

  end subroutine ih_driver_init

  !! These auxiliary subroutines set the values of the EPS, MU, and SIGMA
  !! property arrays. The property values are computed on the heat transfer
  !! mesh using the current value of temperature, and mapped to the EM mesh.
  !! NB: At the outset, the CONST_* flags should be false, so that the values
  !! are set, but if a property is constant (not varying with temperature)
  !! subsequent calls will skip computing and setting the value.

  subroutine set_em_properties(this)
    use em_properties
    class(ih_driver_data), intent(inout) :: this
    real(r8), allocatable :: value(:) !TODO persistent workspace
    allocate(value(this%ht_mesh%ncell_onP))
    if (.not.this%const_eps) then
      call get_permittivity(value)
      call set_permittivity(this, value)
      this%const_eps = permittivity_is_const()
    end if
    if (.not.this%const_mu) then
      call get_permeability(value)
      call set_permeability(this, value)
      this%const_mu = permeability_is_const()
    end if
    if (.not.this%const_sigma) then
      call get_conductivity(value)
      call set_conductivity(this, value)
      this%const_sigma = conductivity_is_const()
    end if
  end subroutine

  subroutine set_permittivity(this, values)
    use physical_constants, only: vacuum_permittivity
    class(ih_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%eps(:this%em_mesh%ncell_onP), defval=vacuum_permittivity, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%eps)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

  subroutine set_permeability(this, values)
    use physical_constants, only: vacuum_permeability
    class(ih_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%mu(:this%em_mesh%ncell_onP), defval=vacuum_permeability, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%mu)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

  subroutine set_conductivity(this, values)
    class(ih_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%sigma(:this%em_mesh%ncell_onP), defval=0.0_r8, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%sigma)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine

  subroutine update_joule_heat(t)

    use EM_properties
    use electromagnetics_namelist, only: params

    real(r8), intent(in) :: t
    real(r8) :: s
    character(10) :: ss

    !TODO: barrier if EM is not enabled?

    call start_timer('joule heat')

    call set_em_properties(this)

    if (allocated(this%q_data)) then ! there is an existing Q that may be usable

      !! Determine if the restart Joule heat is usable; if not, compute it.
      if (this%src_fac%source_is_zero(t)) then
        if (this%src_fac%source_differs(this%q_data, t)) then
          call TLS_info('magnetic source field changed to zero; Joule heat set to zero.')
          call zero_joule_power_density(this, t)
        end if
      else if (matl_prop_differ(this)) then
        call TLS_info('EM properties have changed; computing the Joule heat ...')
        call compute_joule_heat(this, t, params)
      else if (this%src_fac%source_differs(this%q_data, t)) then
        if (this%src_fac%source_is_scaled(this%q_data, t, s)) then
          write(ss,fmt='(es9.3)') s
          call TLS_info('magnetic source field was scaled by ' // trim(ss) // '; Joule heat scaled accordingly.')
          call scale_joule_power_density(this, t, s)
        else
          call TLS_info('magnetic source field changed; computing the Joule heat ...')
          call compute_joule_heat(this, t, params)
        end if
      end if

    else ! just compute the Joule heat

      if (this%src_fac%source_is_zero(t)) then
        call TLS_info('no magnetic source field; setting the Joule heat to zero.')
        call zero_joule_power_density(this, t)
      else
        call TLS_info('computing the Joule heat ...')
        call compute_joule_heat(this, t, params)
      end if

    end if

    if (.not.this%q_written) then
      call danu_write_joule(t)
      this%q_written = .true.
    end if

    call stop_timer('joule heat')

  end subroutine update_joule_heat

  !! Compute the Joule heat by solving Maxwell's equations.

  subroutine compute_joule_heat(this, t, params)

    use parameter_list_type
    use em_bc_factory_type

    class(ih_driver_data), intent(inout), target :: this
    real(r8), intent(in) :: t
    type(parameter_list), intent(inout) :: params

    real(r8), allocatable :: q(:)
    type(em_bc_factory) :: bc_fac
    real(r8) :: freq

    call start_timer('simulation')

    allocate(q(this%em_mesh%ncell))
    call this%src_fac%set_time(t)
    freq = this%src_fac%source_frequency()
    call bc_fac%init(this%em_mesh, this%src_fac, params)
    call compute_joule_heat_emtd(this%em_mesh, freq, this%eps, this%mu, this%sigma, bc_fac, params, q)
    call set_joule_power_density(q)
    this%q_data = this%src_fac%source_data(t)
    this%q_eps = this%eps
    this%q_mu = this%mu
    this%q_sigma = this%sigma
    this%q_written = .false.
    this%q_restart = .false.

    call stop_timer('simulation')

  contains

    subroutine set_joule_power_density(values)
      real(r8), intent(in) :: values(:)
      call start_timer('mesh-to-mesh mapping')
      call this%ht2em%map_field(values(:this%em_mesh%ncell_onP), this%q, defval=0.0_r8, &
                           map_type=GLOBALLY_CONSERVATIVE, pullback=.true.)
      call stop_timer('mesh-to-mesh mapping')
    end subroutine

  end subroutine compute_joule_heat

  !! Computes the Joule heat by solving the time domain Maxwell equations over
  !! several cycles of the external magnetic source field until a periodic
  !! steady state is achieved, and calculating the time-averaged Joule heat
  !! over one period.

  subroutine compute_joule_heat_emtd(mesh, freq, eps, mu, sigma, bc_fac, params, q)

    use simpl_mesh_type
    use parameter_list_type
    use tdme_joule_heat_sim_type
    use em_bc_factory_type

    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: freq, eps(:), mu(:), sigma(:)
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    real(r8), intent(out) :: q(:)

    type(tdme_joule_heat_sim) :: sim
    integer :: stat
    character(:), allocatable :: errmsg

    call sim%init(mesh, freq, eps, mu, sigma, bc_fac, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    call sim%compute(q, stat, errmsg)
    if (stat < 0) then
      call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    else if (stat > 0) then
      call TLS_warn('COMPUTE_JOULE_HEAT: ' // errmsg // '; continuing anyway.')
    end if

  end subroutine compute_joule_heat_emtd

  !! Set the Joule heat to 0. Used when the magnetic source is 0 and
  !! the actual computation is unnecessary.

  subroutine zero_joule_power_density(this, t)
    class(ih_driver_data), intent(inout) :: this
    real(r8), intent(in) :: t
    this%q = 0.0_r8
    this%q_data = this%src_fac%source_data(t)
    this%q_eps = this%eps
    this%q_mu = this%mu
    this%q_sigma = this%sigma
    this%q_written = .false.
    this%q_restart = .false.
  end subroutine

  !! Scales the existing Joule heat by the factor S**2. Used when the strength
  !! of the magnetic source is scaled by the factor S, and there are no other
  !! differences to the input (EM properties, frequency). In this case the new
  !! solution is known and actual computation is unnecessary.

  subroutine scale_joule_power_density(this, t, s)
    class(ih_driver_data), intent(inout) :: this
    real(r8), intent(in) :: t, s
    this%q = s**2 * this%q
    this%q_data = this%src_fac%source_data(t)
    !NB: we do not overwrite q_eps, q_mu, and q_sigma
    this%q_written = .false.
    this%q_restart = .false.
  end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_EM_SIMULATION_ON_OR_OFF, EM_IS_ON, EM_IS_OFF
 !!
 !! Set and query the state of the EM simulation.  When setting the state,
 !! only the value of the argument on the IO processor is relevant.
 !!

  subroutine set_EM_simulation_on_or_off (on)
    use parallel_communication, only : broadcast
    logical, intent(in) :: on
    EM_enabled = on
    call broadcast (EM_enabled)
  end subroutine

  logical function EM_is_on ()
    EM_is_on = EM_enabled
  end function

  logical function ih_enabled()
    ih_enabled = allocated(this)
  end function

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Truchas-side data retrieval procedures
 !!
 !! These functions return pointers to previously stored module data, which are
 !! cell-based arrays on the hex mesh.  Functions may be used in expressions,
 !! or as the target of a pointer assignment, however the target of the result
 !! should never be deallocated.
 !!

  function joule_power_density() result(ptr)
    real(r8), pointer :: ptr(:)
    ptr => this%q
  end function

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MATERIAL_HAS_CHANGED
 !!
 !! This procedure returns the value true if the EM material parameters
 !! returned by CONDUCTIVITY and PERMEABILITY differ significantly from the
 !! values used to compute the Joule heat returned by JOULE_POWER_DENSITY;
 !! otherwise the procedure returns the value false.
 !!
 !! The maximum relative change is taken as the difference measure. Only
 !! when this difference exceeds MATERIAL_CHANGE_THRESHOLD is the difference
 !! considered significant.
 !!
 !! For conductivity, only the conducting region (sigma>0) is considered.
 !! An underlying assumption is that the conducting region remains fixed
 !! throughout the simulation.
 !!
 !! The permittivity only enters the equations through the displacement
 !! current term, which is exceedingly small in this magnetostatic regime
 !! and ought to be dropped entirely.  Thus the solution is essentially
 !! independent of the permittivity and so ignore any changes in its value.
 !!

  !! This function returns true if the EM properties (permeability and
  !! conductivity) have changed significantly since that last time the
  !! Joule heat was computed; otherwise it returns false. The maximum
  !! relative difference in a property over the cells is taken as the
  !! measure of the change, and is considered significant if it exceeds
  !! a user-specified threshold.
  !!
  !! NB: permittivity is not considered since the displacement current
  !! term is an insignificant perturbation to the system in the low
  !! frequency, quasi-magnetostatic regime where induction heating occurs.
  !! In fact the permittivity may be modified for numerical purposes.
  !!
  !! NB: At initialization const_prop is false, ensuring

  logical function matl_prop_differ(this) result(differ)

    use parallel_communication, only: global_maxval

    class(ih_driver_data), intent(inout) :: this

    integer :: j
    real(r8) :: dmu, dsigma
    character(80) :: string

    differ = .false.

    if (.not.this%const_mu .or. this%q_restart) then
      dmu = global_maxval(abs(this%mu-this%q_mu)/this%mu)
      differ = .true.
    else
      dmu = 0.0_r8
    end if

    if (.not.this%const_sigma .or. this%q_restart) then
      dsigma = 0.0_r8
      do j = 1, size(this%sigma)
        if (this%sigma(j) > 0.0_r8) dsigma = max(dsigma, abs(this%sigma(j)-this%q_sigma(j))/this%sigma(j))
      end do
      dsigma = global_maxval(dsigma)
      differ = .true.
    else
      dsigma = 0.0_r8
    end if

    if (.not.differ) return ! nothing to check

    if (this%q_restart) then
      differ = (max(dmu, dsigma) > 0.0_r8)
      if (.not.differ) this%q_restart = .false. ! properties confirmed to be the same
    else
      differ = (max(dmu, dsigma) > this%matl_change_threshold)
    end if

    write(string,fmt='(3x,2(a,es10.3))') 'max relative change: sigma=', dsigma, ', mu=', dmu
    call TLS_info(string)

  end function matl_prop_differ

  !! Read the Joule heat segment of the restart file opened (and prepositioned)
  !! on UNIT, and initialize the Q, Q_DATA, Q_MU, AND Q_SIGMA components.
  !! If any errors or incompatibilities with data read from the input file
  !! are encountered, execution is halted.

  subroutine read_joule_data(unit, version)

    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c
    use parallel_communication, only: global_sum

    integer, intent(in) :: unit,  version

    integer :: n
    character(32) :: coil_md5sum

    call TLS_info ('  reading the Joule heat data from the restart file')

    call read_var(unit, coil_md5sum, 'READ_JOULE_DATA: error reading COIL_MD5SUM record')
    if (coil_md5sum /= this%coil_md5sum) call halt('READ_JOULE_DATA: incompatible coil geometry')
    call read_var(unit, n, 'READ_JOULE_DATA: error reading NQDATA record')
    allocate(this%q_data(n)) ! certainly the correct size because of matching MD5 sums
    call read_var(unit, this%q_data, 'READ_JOULE_DATA: error reading QDATA record')

    call read_var(unit, n, 'READ_JOULE_DATA: error reading NMU record')
    if (n /= this%em_mesh%cell_imap%global_size) &
        call halt('READ_JOULE_DATA: incompatible NMU value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_mu(:this%em_mesh%ncell_onP), &
        this%em_mesh%xcell(:this%em_mesh%ncell_onP), &
        'READ_JOULE_DATA: error reading MU record')
    call this%em_mesh%cell_imap%gather_offp(this%q_mu)

    call read_var(unit, n, 'READ_JOULE_DATA: error reading NSIGMA record')
    if (n /= this%em_mesh%cell_imap%global_size) &
        call halt('READ_JOULE_DATA: incompatible NSIGMA value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_sigma(:this%em_mesh%ncell_onP), &
        this%em_mesh%xcell(:this%em_mesh%ncell_onP), &
        'READ_JOULE_DATA: error reading SIGMA record')
    call this%em_mesh%cell_imap%gather_offp(this%q_sigma)

    !TODO: joule size is only on-process cells?
    call read_var(unit, n, 'READ_JOULE_DATA: error reading NJOULE record')
    if (n /= global_sum(size(this%q))) &
        call halt('READ_JOULE_DATA: incompatible NJOULE value: ' // i_to_c(n))
    call read_dist_array(unit, this%q, &
        this%ht_mesh%xcell(:this%ht_mesh%ncell_onp), &
        'READ_JOULE_DATA: error reading JOULE record')
    this%q_written = .false.
    this%q_restart = .true.

  end subroutine read_joule_data

  !! Advance over the Joule heat segment of the restart file opened (and
  !! prepositioned) on UNIT. This is used when the restart file contains Joule
  !! heat data but induction heating is not enabled. NB: the number of records
  !! skipped must match the number of records read by READ_JOULE_DATA.

  subroutine skip_joule_data(unit, version)
    use restart_utilities, only: skip_records
    integer, intent(in) :: unit, version
    call skip_records(unit, 9, 'SKIP_JOULE_DATA: error skipping the Joule heat data')
  end subroutine

  !! Write the Joule heat data and inputs used to compute it to the HDF output
  !! file. This is checkpoint data used to generate restart files. Output for
  !! visualization is written elsewhere.

  subroutine danu_write_joule(t)

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

    call outfile%reopen(output_file_name('h5'), io_group_size, is_IOP)

    n = global_sum(this%em_mesh%ncell_onP)

    !! NNC, 2/2017. Not entirely sure why the reordering for MU and SIGMA is
    !! necessary.  I believe if we wrote the cell map for the tet mesh, then
    !! post-processing tools could do it when needed, as is done for the main
    !! mesh, and we could dispense with the collation here and truly write in
    !! parallel.  FIXME

    !! Collate the cell permutation array.
    allocate(cell_perm(merge(n,0,is_IOP)))
    call gather(this%em_mesh%xcell(:this%em_mesh%ncell_onP), cell_perm)

    !! Collate the cell-based MU array on the tet mesh, and restore it to the external order.
    allocate(col_mu(merge(n,0,is_IOP)))
    call gather(this%q_mu(:this%em_mesh%ncell_onP), col_mu)
    if (is_IOP) call reorder(col_mu, cell_perm, forward=.true.)

    !! Collate the cell-based SIGMA array on the tet mesh, and restore it to the external order.
    allocate(col_sigma(merge(n,0,is_IOP)))
    call gather(this%q_sigma(:this%em_mesh%ncell_onP), col_sigma)
    if (is_IOP) call reorder(col_sigma, cell_perm, forward=.true.)

    !! Write the data.
    sim_num = sim_num + 1
    write(sim_name,'(a,i3.3)') 'EM', sim_num
    call outfile%add_sim_group(trim(sim_name), sim)
    call sim%write_attr('TIME', t)
    call TLS_info('  writing EM restart data for ' // trim(sim_name))
    call sim%write_attr('COIL_MD5SUM', this%coil_md5sum)
    call sim%write_repl_data('SOURCE_DATA', this%q_data)
    call sim%write_repl_data('MU', col_mu)
    call sim%write_repl_data('SIGMA', col_sigma)
    call sim%write_dist_array('JOULE', this%q, global_sum(size(this%q)))

    call outfile%close

  end subroutine danu_write_joule

end module ih_driver
