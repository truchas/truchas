!!
!! IH_DRIVER
!!
!! Procedures for driving the computation of the Joule heat for induction
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

module ih_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use simpl_mesh_type
  use data_mapper_class
  use ih_source_factory_type
  use parameter_list_type
  use parallel_communication
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  public :: ih_enabled, ih_driver_init, ih_driver_final
  public :: update_joule_heat, joule_power_density
  public :: read_ih_namelists, read_joule_data, skip_joule_data

  type :: ih_driver_data
    class(base_mesh), pointer :: ht_mesh => null() ! unowned reference
    type(simpl_mesh), pointer :: em_mesh => null() ! unowned reference
    class(data_mapper), allocatable :: ht2em
    type(ih_source_factory) :: src_fac
    logical :: use_emfd_solver
    character(32) :: coil_md5sum  ! fingerprint of the fixed coil geometry
    ! EM properties
    logical :: const_eps=.false., const_mu=.false., const_sigma=.false.
    real(r8) :: matl_change_threshold
    real(r8), allocatable :: eps(:), epsi(:), mu(:), sigma(:)  ! on EM mesh
    ! Computed Joule heat and the variable inputs used
    real(r8), allocatable :: q(:) ! on HT mesh
    real(r8), allocatable :: q_data(:), q_eps(:), q_mu(:), q_sigma(:)
    logical :: q_written = .false.
    logical :: q_restart = .false. ! true when q was read from a restart file
  end type

  type(ih_driver_data), allocatable, target :: this
  type(parameter_list), pointer :: params

contains

  logical function ih_enabled()
    use physics_module, only: electromagnetics
    ih_enabled = electromagnetics
  end function

  subroutine ih_driver_final
    if (allocated(this)) deallocate(this)
  end subroutine

  function joule_power_density() result(ptr)
    real(r8), pointer :: ptr(:)
    ptr => this%q
  end function

  !! This driver subroutine is called by the Truchas input driver and it
  !! coordinates the reading of all the namelists pertaining to induction
  !! heating, and assembles the input into the PARAMS parameter list held
  !! as a private module variable.

  subroutine read_ih_namelists(lun)
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

  subroutine ih_driver_init(t)

    use mesh_manager, only: named_mesh_ptr, simpl_mesh_ptr
    use em_properties

    real(r8), intent(in) :: t

    integer :: n, stat
    character(:), allocatable :: errmsg, data_mapper_kind
    type(parameter_list), pointer :: plist

    ASSERT(.not.allocated(this))

    call TLS_info('')
    call TLS_info('Initializing induction heating solver ...')

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

    allocate(this%q(this%ht_mesh%ncell_onP))
    !TODO: NEED TO ADD Q_EPSI AND INCLUDE WITH RESTART DATA
    allocate(this%q_eps(n), this%q_mu(n), this%q_sigma(n))
    !NB: the allocation status of Q_DATA indicates whether Q contains data

    !! Initialize the external H-field source factory
    plist => params%sublist('external-field')
    call this%src_fac%init(plist, t, stat, errmsg)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    this%coil_md5sum = this%src_fac%coil_geom_fingerprint()

    call define_default_em_properties

    call params%get('matl-change-threshold', this%matl_change_threshold, stat, errmsg, default=0.3_r8)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    if (this%matl_change_threshold <= 0.0_r8) call TLS_fatal('IH_DRIVER_INIT: matl-change-threshold must be > 0.0')

    call params%get('frequency-domain-solver', this%use_emfd_solver, stat=stat, errmsg=errmsg, default=.false.)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)

  end subroutine ih_driver_init

  !! This primary driver subroutine is called every heat transfer time step
  !! to update the Joule heat stored by the driver as needed in response to
  !! evolving conditions (time, temperature). The driver seeks to avoid
  !! computing the Joule heat when possible by zeroing the Joule heat, scaling
  !! it, or leaving it as is. Only when necessary does it invoke the
  !! computation of the Joule heat.

  subroutine update_joule_heat(t)

    use em_properties

    real(r8), intent(in) :: t
    real(r8) :: s
    character(10) :: ss

    ASSERT(allocated(this))

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

  !! This auxiliary subroutine computes and stores the Joule heat using the
  !! current values of the EM properties and external magnetic source.

  subroutine compute_joule_heat(this, t, params)

    use em_bc_factory_type
    use string_utilities, only: i_to_c
    use truchas_env, only: output_dir
    use physical_constants, only: vacuum_permittivity, vacuum_permeability

    class(ih_driver_data), intent(inout), target :: this
    real(r8), intent(in) :: t
    type(parameter_list), intent(inout) :: params

    real(r8), allocatable :: q(:)
    type(em_bc_factory) :: bc_fac
    real(r8) :: freq
    type(parameter_list), pointer :: plist
    integer, save :: sim_num = 0  ! global counter for the number of calls
    real(r8), parameter :: PI = 3.1415926535897932385_r8

    call start_timer('simulation')

    sim_num = sim_num + 1
    call params%set('graphics-file', trim(output_dir)//'tdme-'//i_to_c(sim_num)//'.vtkhdf')

    allocate(q(this%em_mesh%ncell_onP))
    call this%src_fac%set_time(t)
    freq = this%src_fac%H_freq()
    call params%set('omega', 2*PI*freq)
    call params%set('epsilon_0', vacuum_permittivity)
    call params%set('mu_0', vacuum_permeability)
    call bc_fac%init(this%em_mesh, this%src_fac, params)
    if (this%use_emfd_solver) then
      plist => params%sublist('emfd-solver')
      call plist%set('graphics-file', trim(output_dir)//'fdme-'//i_to_c(sim_num)//'.vtkhdf')
      call compute_joule_heat_emfd(this%em_mesh, freq, this%eps, this%epsi, this%mu, this%sigma, bc_fac, plist, q)
    else
      call compute_joule_heat_emtd(this%em_mesh, freq, this%eps, this%mu, this%sigma, bc_fac, params, q)
    end if
    call set_joule_power_density(q)
    this%q_data = this%src_fac%source_data(t)
    this%q_eps = this%eps
    this%q_mu = this%mu
    this%q_sigma = this%sigma
    this%q_written = .false.
    this%q_restart = .false.

    !TODO: make output subject to verbosity level
    block
      character(80) :: string
      write(string,fmt='(2(a,es11.4))') '|Q|_max=', global_maxval(q), &
          ', Q_total=', global_dot_product(q, abs(this%em_mesh%volume(:size(q))))
      call tls_info(trim(string))
    end block

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

  !! Compute Joule heat using the frequency domain EM solver

  subroutine compute_joule_heat_emfd(mesh, freq, eps, epsi, mu, sigma, bc_fac, params, q)

    use simpl_mesh_type
    use em_bc_factory_type
    use fdme_model_type
    use fdme_zvector_type
    use parameter_list_type
    use physical_constants, only: vacuum_permittivity, vacuum_permeability
    use emfd_nlsol_solver_type

    type(simpl_mesh), intent(inout), target :: mesh
    real(r8), intent(in) :: freq, eps(:), epsi(:), mu(:), sigma(:)
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    real(r8), intent(out) :: q(:)

    type(fdme_model), target :: model
    type(emfd_nlsol_solver) :: solver
    real(r8), parameter :: PI = 3.1415926535897932385_r8
    real(r8) :: t, omega
    complex(r8) :: efield(mesh%nedge)
    logical :: flag
    integer :: stat
    character(:), allocatable :: errmsg, filename

    call model%init(mesh, bc_fac, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)

    call solver%init(model, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)

    t = 0 ! dummy time
    omega = 2 * PI * freq

    !TODO? rework solver to use absolute eps and mu?
    call model%setup(t, eps/vacuum_permittivity, epsi/vacuum_permittivity, mu/vacuum_permeability, sigma, omega)

    efield = 0 ! initial guess
    call solver%solve(efield, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    call mesh%edge_imap%gather_offp(efield)

    call model%compute_heat_source(efield, q)

    !! Graphics output
    call params%get('graphics-output', flag, stat, errmsg, default=.false.)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    if (flag) then
      block
        complex(r8) :: bfield(mesh%nface), div_dfield(mesh%nnode)
        call model%compute_bfield(efield, bfield)
        call model%compute_div(efield, div_dfield)
        call params%get('graphics-file', filename, stat, errmsg)
        if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
        call emfd_vtk_graphics(filename, mesh, q, efield, bfield, mu, div_dfield, stat, errmsg)
        if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
      end block
    end if

  end subroutine compute_joule_heat_emfd

  subroutine emfd_vtk_graphics(filename, mesh, qfield, efield, bfield, mu, div_dfield, stat, errmsg)

    use vtkhdf_file_type
    use mimetic_discretization, only: w1_vector_on_cells, w2_vector_on_cells

    character(*), intent(in) :: filename
    type(simpl_mesh), intent(in) :: mesh
    real(r8), intent(in) :: qfield(:), mu(:)
    complex(r8), intent(inout) :: efield(:), bfield(:), div_dfield(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j
    type(vtkhdf_file) :: viz_file
    real(r8), allocatable :: g_scalar(:)
    complex(r8), allocatable :: g_vector(:,:), l_vector(:,:), g_zscalar(:)

    if (is_IOP) call viz_file%create(filename, stat, errmsg)
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast(errmsg)
      return
    end if

    call export_mesh

    allocate(g_scalar(merge(mesh%cell_imap%global_size, 0, is_IOP)))
    call gather(qfield(:mesh%ncell_onP), g_scalar)
    if (is_IOP) call viz_file%write_cell_dataset('Q_EM', g_scalar, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    allocate(g_vector(3,merge(mesh%cell_imap%global_size, 0, is_IOP)))
    allocate(l_vector(3,mesh%ncell))

    l_vector(:,:)%re = w1_vector_on_cells(mesh, efield%re)
    l_vector(:,:)%im = w1_vector_on_cells(mesh, efield%im)
    call gather(l_vector(:,:mesh%ncell_onP), g_vector)
    if (is_IOP) call viz_file%write_cell_dataset('E_re', g_vector%re, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

#ifdef GNU_PR117774
    if (is_IOP) call viz_file%write_cell_dataset('E_im', reshape([g_vector%im],shape(g_vector)), stat, errmsg)
#else
    if (is_IOP) call viz_file%write_cell_dataset('E_im', g_vector%im, stat, errmsg)
#endif
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call viz_file%write_cell_dataset('|E|', abs(g_vector), stat, errmsg)

    l_vector(:,:)%re = w2_vector_on_cells(mesh, bfield%re)
    l_vector(:,:)%im = w2_vector_on_cells(mesh, bfield%im)
    do j = 1, size(mu) ! convert cell-centered B to H
      l_vector(:,j) = l_vector(:,j) / mu(j)
    end do
    call gather(l_vector(:,:mesh%ncell_onP), g_vector)
    if (is_IOP) call viz_file%write_cell_dataset('H_re', g_vector%re, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

#ifdef GNU_PR117774
    if (is_IOP) call viz_file%write_cell_dataset('H_im', reshape([g_vector%im],shape(g_vector)), stat, errmsg)
#else
    if (is_IOP) call viz_file%write_cell_dataset('H_im', g_vector%im, stat, errmsg)
#endif
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call viz_file%write_cell_dataset('|H|', abs(g_vector), stat, errmsg)

    !! Output the mesh partition
    call gather(spread(real(this_PE,kind=r8), dim=1, ncopies=mesh%ncell_onP), g_scalar)
    if (is_IOP) call viz_file%write_cell_dataset('MPI rank', g_scalar, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    !! Divergence of the electric flux
    allocate(g_zscalar(merge(mesh%node_imap%global_size, 0, is_IOP)))
    call gather(div_dfield(:mesh%nnode_onP), g_zscalar)
    if (is_IOP) call viz_file%write_point_dataset('div_D_re', g_zscalar%re, stat, errmsg)
    call broadcast(stat)
#ifdef GNU_PR117774
    if (is_IOP) call viz_file%write_point_dataset('div_D_im', [g_zscalar%im], stat, errmsg)
#else
    if (is_IOP) call viz_file%write_point_dataset('div_D_im', g_zscalar%im, stat, errmsg)
#endif
    call broadcast(stat)
    if (is_IOP) call viz_file%write_point_dataset('|div_D|', abs(g_zscalar), stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call viz_file%close

  contains

    subroutine export_mesh

      use,intrinsic :: iso_fortran_env, only: int8

      integer, allocatable, target :: cnode(:,:)
      integer, allocatable :: xcnode(:)
      integer(int8), allocatable :: types(:)
      real(r8), allocatable :: x(:,:)
      integer, pointer :: connectivity(:)
      integer :: j, stat
      character(:), allocatable :: errmsg

      !! Collate the mesh data structure onto the IO process
      call mesh%get_global_cnode_array(cnode)
      call mesh%get_global_x_array(x)

      if (is_IOP) then
        xcnode = [(1+4*j, j=0, size(cnode,dim=2))]
        connectivity(1:size(cnode)) => cnode ! flattened view
        types = spread(VTK_TETRA, dim=1, ncopies=size(cnode,dim=2))
        call viz_file%write_mesh(x, connectivity, xcnode, types, stat, errmsg)
      end if
      call broadcast(stat)
      INSIST(stat == 0)

    end subroutine export_mesh

  end subroutine emfd_vtk_graphics

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

  !! These auxiliary subroutines set the values of the EPS, MU, and SIGMA
  !! property arrays. The property values are computed on the heat transfer
  !! mesh using the current value of temperature, and mapped to the EM mesh.
  !!
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
      call get_permittivity_im(value)
      call set_permittivity_im(this, value)
      this%const_eps = permittivity_is_const() .and. permittivity_im_is_const()
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

  subroutine set_permittivity_im(this, values)
    use physical_constants, only: vacuum_permittivity
    class(ih_driver_data), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    call start_timer('mesh-to-mesh mapping')
    call this%ht2em%map_field(values, this%epsi(:this%em_mesh%ncell_onP), defval=0.0_r8, map_type=LOCALLY_BOUNDED)
    call this%em_mesh%cell_imap%gather_offp(this%epsi)
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

  !! This auxiliary function returns true if the EM properties (permeability
  !! and conductivity) have changed significantly since that last time the
  !! Joule heat was computed; otherwise it returns false. The maximum relative
  !! difference in a property over the cells is taken as the measure of the
  !! change, and is considered significant if it exceeds a user-specified
  !! threshold.
  !!
  !! NB: Constant properties do not normally need to be checked. The exception
  !! is when the stored joule heat and its properties came from a restart file.
  !!
  !! NB: permittivity is not considered since the displacement current
  !! term is an insignificant perturbation to the system in the low
  !! frequency, quasi-magnetostatic regime where induction heating occurs.
  !! In fact the permittivity may be modified for numerical purposes.
  !!
  !! TODO: In the high frequency regime we need to consider changes to the
  !! complex permittivity for the frequency domain solver.

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

    ASSERT(allocated(this))

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

    ASSERT(allocated(this))

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
