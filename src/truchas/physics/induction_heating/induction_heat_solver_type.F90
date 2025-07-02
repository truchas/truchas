#include "f90_assert.fpp"

module induction_heat_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use ih_hfield_factory_type
  use parameter_list_type
  use parallel_communication
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: induction_heat_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
    type(ih_hfield_factory) :: hfield_fac
    real(r8) :: matl_change_threshold
    logical  :: const_mu=.false., const_sigma=.false.
    real(r8), allocatable :: q_data(:), q_eps(:), q_mu(:), q_sigma(:)
    logical :: use_fd_solver
    logical :: q_restart = .false. ! true when q was read from a restart file
    character(32) :: coil_md5sum  ! fingerprint of the fixed coil geometry
  contains
    procedure :: init
    procedure :: update_em_heat => update_joule_heat
    procedure :: write_restart_data, read_restart_data, skip_restart_data
  end type

contains

  subroutine init(this, mesh, params)

    use em_properties

    class(induction_heat_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout), target :: params

    integer :: n, stat
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist

    this%mesh => mesh
    this%params => params

    this%const_mu = permeability_is_const()
    this%const_sigma = conductivity_is_const()

    !TODO: NEED TO ADD Q_EPSI AND INCLUDE WITH RESTART DATA
    n = this%mesh%ncell
    allocate(this%q_eps(n), this%q_mu(n), this%q_sigma(n))
    !NB: the allocation status of Q_DATA indicates whether Q contains data

    !! Initialize the external H-field source factory
    plist => params%sublist('external-field')
    call this%hfield_fac%init(plist, stat, errmsg)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    this%coil_md5sum = this%hfield_fac%coil_geom_fingerprint()

    call params%get('matl-change-threshold', this%matl_change_threshold, stat, errmsg, default=0.3_r8)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)
    if (this%matl_change_threshold <= 0.0_r8) call TLS_fatal('IH_DRIVER_INIT: matl-change-threshold must be > 0.0')

    call params%get('use-fd-solver', this%use_fd_solver, stat=stat, errmsg=errmsg, default=.false.)
    if (stat /= 0) call TLS_fatal('IH_DRIVER_INIT: ' // errmsg)

  end subroutine init

  !! This primary driver subroutine is called every heat transfer time step
  !! to update the Joule heat stored by the driver as needed in response to
  !! evolving conditions (time, temperature). The driver seeks to avoid
  !! computing the Joule heat when possible by zeroing the Joule heat, scaling
  !! it, or leaving it as is. Only when necessary does it invoke the
  !! computation of the Joule heat.

  subroutine update_joule_heat(this, t, eps, epsi, mu, sigma, q, stat)

    class(induction_heat_solver), intent(inout) :: this
    real(r8), intent(in) :: t, eps(:), epsi(:), mu(:), sigma(:)
    real(r8), intent(inout) :: q(:)
    integer, intent(out) :: stat

    real(r8) :: s, freq
    character(10) :: ss

    enum, bind(c)
      enumerator :: KEEP=0, ZERO, SCALE, COMPUTE
    end enum

    if (allocated(this%q_data)) then ! there is an existing Q that may be usable
      if (this%hfield_fac%hfield_is_zero(t)) then
        stat = merge(ZERO, KEEP, this%hfield_fac%hfield_differs(this%q_data, t))
      else if (matl_prop_differ(this, eps, mu, sigma)) then
        stat = COMPUTE
      else if (this%hfield_fac%hfield_differs(this%q_data, t)) then
        stat = merge(SCALE, COMPUTE, this%hfield_fac%hfield_is_scaled(this%q_data, t, s))
      end if
    else ! no existing Q exists; set it
      stat = merge(ZERO, COMPUTE, this%hfield_fac%hfield_is_zero(t))
    end if

    select case (stat)
    case (KEEP)
      return
    case (ZERO)
      call TLS_info('no magnetic source field; setting the Joule heat to zero.')
      call zero_joule_power_density
    case (SCALE)
      write(ss,fmt='(es9.3)') s
      call TLS_info('magnetic source field was scaled by ' // trim(ss) // '; Joule heat scaled accordingly.')
      call scale_joule_power_density
    case (COMPUTE)
      call TLS_info('(re)computing the Joule heat ...')
      call this%hfield_fac%update_ih_hfield_func(t)
      freq = this%hfield_fac%frequency(t)
      call compute_joule_heat(this, freq, eps, epsi, mu, sigma, this%params, q)
      this%q_data = this%hfield_fac%hfield_data(t)
      this%q_eps = eps
      this%q_mu = mu
      this%q_sigma = sigma
      this%q_restart = .false.
    end select

  contains

    !! Set the Joule heat to 0. Used when the magnetic source is 0 and
    !! the actual computation is unnecessary.

    subroutine zero_joule_power_density
      q = 0.0_r8
      this%q_data = this%hfield_fac%hfield_data(t)
      this%q_eps = eps
      this%q_mu = mu
      this%q_sigma = sigma
      this%q_restart = .false.
    end subroutine

    !! Scales the existing Joule heat by the factor S**2. Used when the strength
    !! of the magnetic source is scaled by the factor S, and there are no other
    !! differences to the input (EM properties, frequency). In this case the new
    !! solution is known and actual computation is unnecessary.

    subroutine scale_joule_power_density
      q = s**2 * q
      this%q_data = this%hfield_fac%hfield_data(t)
      !NB: we do not overwrite q_eps, q_mu, and q_sigma
      this%q_restart = .false.
    end subroutine

  end subroutine update_joule_heat

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

  logical function matl_prop_differ(this, eps, mu, sigma) result(differ)

    use parallel_communication, only: global_maxval

    class(induction_heat_solver), intent(inout) :: this
    real(r8), intent(in) :: eps(:), mu(:), sigma(:)

    integer :: j
    real(r8) :: dmu, dsigma
    character(80) :: string

    differ = .false.

    if (.not.this%const_mu .or. this%q_restart) then
      dmu = global_maxval(abs(mu-this%q_mu)/mu)
      differ = .true.
    else
      dmu = 0.0_r8
    end if

    if (.not.this%const_sigma .or. this%q_restart) then
      dsigma = 0.0_r8
      do j = 1, size(sigma)
        if (sigma(j) > 0.0_r8) dsigma = max(dsigma, abs(sigma(j)-this%q_sigma(j))/sigma(j))
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

  !! This auxiliary subroutine computes and stores the Joule heat using the
  !! current values of the EM properties and external magnetic source.

  subroutine compute_joule_heat(this, freq, eps, epsi, mu, sigma, params, q)

    use em_bc_factory_type
    use string_utilities, only: i_to_c
    use truchas_env, only: output_dir

    class(induction_heat_solver), intent(inout), target :: this
    real(r8), intent(in) :: freq, eps(:), epsi(:), mu(:), sigma(:)
    type(parameter_list), intent(inout) :: params
    real(r8), intent(out) :: q(:)

    integer, save :: sim_num = 0  ! global counter for the number of calls

    call start_timer('simulation')

    sim_num = sim_num + 1

    if (this%use_fd_solver) then
      call params%set('graphics-file', trim(output_dir)//'fdme-'//i_to_c(sim_num)//'.vtkhdf')
      call compute_joule_heat_fdme(this%mesh, freq, eps, epsi, mu, sigma, params, q)
    else
      call params%set('graphics-file', trim(output_dir)//'tdme-'//i_to_c(sim_num)//'.vtkhdf')
      call compute_joule_heat_tdme(this%mesh, freq, eps, mu, sigma, params, q)
    end if

    !TODO: make output subject to verbosity level
    block
      character(80) :: string
      write(string,fmt='(2(a,es11.4))') '|Q|_max=', global_maxval(q), &
          ', Q_total=', global_dot_product(q, abs(this%mesh%volume(:size(q))))
      call tls_info(trim(string))
    end block

    call stop_timer('simulation')

  end subroutine compute_joule_heat

  !! Computes the Joule heat by solving the time domain Maxwell equations over
  !! several cycles of the external magnetic source field until a periodic
  !! steady state is achieved, and calculating the time-averaged Joule heat
  !! over one period.

  subroutine compute_joule_heat_tdme(mesh, freq, eps, mu, sigma, params, q)

    use simpl_mesh_type
    use parameter_list_type
    use tdme_joule_heat_sim_type
    use em_bc_factory_type
    use physical_constants, only: vacuum_permittivity, vacuum_permeability

    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: freq, eps(:), mu(:), sigma(:)
    type(parameter_list), intent(inout) :: params
    real(r8), intent(out) :: q(:)

    type(em_bc_factory) :: bc_fac
    type(tdme_joule_heat_sim) :: sim
    integer :: stat
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist
    real(r8) :: omega
    logical :: flag

    omega = 8*atan(1.0_r8)*freq
    call params%get('use-legacy-bc', flag, stat, errmsg, default=.false.)
    if (stat /= 0) return
    plist => params%sublist('bc')
    call bc_fac%init(mesh, omega, plist, use_legacy_bc=flag)

    call sim%init(mesh, freq, vacuum_permittivity*eps, vacuum_permeability*mu, sigma, bc_fac, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    call sim%compute(q, stat, errmsg)
    if (stat < 0) then
      call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    else if (stat > 0) then
      call TLS_warn('COMPUTE_JOULE_HEAT: ' // errmsg // '; continuing anyway.')
    end if

  end subroutine compute_joule_heat_tdme

  !! Compute Joule heat using the frequency domain EM solver

  subroutine compute_joule_heat_fdme(mesh, freq, eps, epsi, mu, sigma, params, q)

    use simpl_mesh_type
    use parameter_list_type
    use fdme_solver_type

    type(simpl_mesh), intent(inout), target :: mesh
    real(r8), intent(in) :: freq, eps(:), epsi(:), mu(:), sigma(:)
    type(parameter_list), intent(inout) :: params
    real(r8), intent(out) :: q(:)

    type(fdme_solver) :: solver
    real(r8) :: omega
    logical :: flag
    integer :: stat
    character(:), allocatable :: errmsg, filename

    omega = 8*atan(1.0_r8)*freq

    call solver%init(mesh, omega, eps, epsi, mu, sigma, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)

    call solver%solve(stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)

    call solver%get_heat_source(q)

    !! Graphics output
    call params%get('graphics-output', flag, stat, errmsg, default=.false.)
    if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    if (flag) then
      call params%get('graphics-file', filename, stat, errmsg)
      if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
      call fdme_vtk_graphics(solver, filename, mesh, q, stat, errmsg)
      if (stat /= 0) call TLS_fatal('COMPUTE_JOULE_HEAT: ' // errmsg)
    end if

  end subroutine compute_joule_heat_fdme


  subroutine write_restart_data(this, sim)

    use truchas_h5_outfile, only: th5_sim_group
    use parallel_communication
    use permutations

    class(induction_heat_solver), intent(in) :: this
    type(th5_sim_group), intent(in) :: sim

    integer :: n
    integer, allocatable :: cell_perm(:)
    real(r8), allocatable :: col_mu(:), col_sigma(:)

    n = global_sum(this%mesh%ncell_onP)

    !! NNC, 2/2017. Not entirely sure why the reordering for MU and SIGMA is
    !! necessary.  I believe if we wrote the cell map for the tet mesh, then
    !! post-processing tools could do it when needed, as is done for the main
    !! mesh, and we could dispense with the collation here and truly write in
    !! parallel.  FIXME

    !! Collate the cell permutation array.
    allocate(cell_perm(merge(n,0,is_IOP)))
    call gather(this%mesh%xcell(:this%mesh%ncell_onP), cell_perm)

    !! Collate the cell-based MU array on the tet mesh, and restore it to the external order.
    allocate(col_mu(merge(n,0,is_IOP)))
    call gather(this%q_mu(:this%mesh%ncell_onP), col_mu)
    if (is_IOP) call reorder(col_mu, cell_perm, forward=.true.)

    !! Collate the cell-based SIGMA array on the tet mesh, and restore it to the external order.
    allocate(col_sigma(merge(n,0,is_IOP)))
    call gather(this%q_sigma(:this%mesh%ncell_onP), col_sigma)
    if (is_IOP) call reorder(col_sigma, cell_perm, forward=.true.)

    call sim%write_attr('COIL_MD5SUM', this%coil_md5sum)
    call sim%write_repl_data('HFIELD_DATA', this%q_data)
    call sim%write_repl_data('MU', col_mu)
    call sim%write_repl_data('SIGMA', col_sigma)

  end subroutine write_restart_data

  !! Read the parameters specific to induction heating segment of the restart
  !! file opened (and prepositioned) on UNIT, and initialize the Q_DATA, Q_MU,
  !! AND Q_SIGMA components. If any errors or incompatibilities with data read
  !! from the input file are encountered, execution is halted.

  subroutine read_restart_data(this, unit, version)

    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c
    use parallel_communication, only: global_sum

    class(induction_heat_solver), intent(inout) :: this
    integer, intent(in) :: unit,  version

    integer :: n
    character(32) :: coil_md5sum

    this%q_restart = .true.

    call read_var(unit, coil_md5sum, 'INDUCTION_HEAT_SOLVER: error reading COIL_MD5SUM restart record')
    if (coil_md5sum /= this%coil_md5sum) call halt('INDUCTION_HEAT_SOLVER: incompatible restart coil geometry')
    call read_var(unit, n, 'INDUCTION_HEAT_SOLVER: error reading NQDATA restart record')
    allocate(this%q_data(n)) ! certainly the correct size because of matching MD5 sums
    call read_var(unit, this%q_data, 'INDUCTION_HEAT_SOLVER: error reading QDATA restart record')

    call read_var(unit, n, 'INDUCTION_HEAT_SOLVER: error reading NMU restart record')
    if (n /= this%mesh%cell_imap%global_size) &
        call halt('INDUCTION_HEAT_SOLVER: incompatible restart NMU value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_mu(:this%mesh%ncell_onP), &
        this%mesh%xcell(:this%mesh%ncell_onP), &
        'INDUCTION_HEAT_SOLVER: error reading MU restart record')
    call this%mesh%cell_imap%gather_offp(this%q_mu)

    call read_var(unit, n, 'INDUCTION_HEAT_SOLVER: error reading NSIGMA restart record')
    if (n /= this%mesh%cell_imap%global_size) &
        call halt('INDUCTION_HEAT_SOLVER: incompatible restart NSIGMA value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_sigma(:this%mesh%ncell_onP), &
        this%mesh%xcell(:this%mesh%ncell_onP), &
        'INDUCTION_HEAT_SOLVER: error reading SIGMA restart record')
    call this%mesh%cell_imap%gather_offp(this%q_sigma)

  end subroutine read_restart_data

  subroutine skip_restart_data(this, unit, version)
    use restart_utilities, only: skip_records
    class(induction_heat_solver), intent(in) :: this
    integer, intent(in) :: unit, version
    call skip_records(unit, 7, 'INDUCTION_HEAT_SOLVER: error skipping the restart data')
  end subroutine

  subroutine fdme_vtk_graphics(solver, filename, mesh, qfield, stat, errmsg)

    use fdme_solver_type
    use vtkhdf_file_type

    type(fdme_solver), intent(in) :: solver
    character(*), intent(in) :: filename
    type(simpl_mesh), intent(in) :: mesh
    real(r8), intent(in) :: qfield(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(vtkhdf_file) :: viz_file
    real(r8), allocatable :: g_scalar(:)
    complex(r8), allocatable :: g_vector(:,:), l_vector(:,:), g_zscalar(:), l_zscalar(:)

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

    call solver%get_cell_efield(l_vector)
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

    call solver%get_cell_hfield(l_vector)
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
    allocate(l_zscalar(mesh%nnode))
    call solver%get_div_dfield(l_zscalar)
    call gather(l_zscalar(:mesh%nnode_onP), g_zscalar)
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

  end subroutine fdme_vtk_graphics

end module induction_heat_solver_type
