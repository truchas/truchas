!!
!! MICROWAVE_HEAT_SOLVER_TYPE
!!
!! This module defines a derived type that drives the evaluation of the time
!! varying EM heat source (primarily dielectric loss) for microwave heating
!! simulations. Time dependence is due to stepwise changes in the waveguide
!! input powers, and changes in temperature-dependent EM properties. When
!! necessary, an external solver for the time harmonic Maxwell equations is
!! invoked to directly compute the EM heat source, but this is avoided when
!! possible by scaling a previously computed source, zeroing the source, or
!! leaving the source unchanged. Additional methods handle writing/reading
!! restart data.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! July 2025
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module microwave_heat_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use parameter_list_type
  use wg_port_bc_plist_factory_type
  use parallel_communication
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: microwave_heat_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
    type(wg_port_bc_plist_factory) :: wg_port_fac
    logical :: const_eps=.false., const_epsi=.false., const_mu=.false.
    logical :: prev_q_is_zero = .false.
    logical :: q_from_restart = .false., check_restart_q_props = .false.
    real(r8) :: omega, prop_change_threshold, prev_q_scf = 0.0_r8
    real(r8), allocatable :: q(:), q_eps(:), q_epsi(:), q_mu(:), q_data(:)
    logical :: graphics_output
    integer :: sim_num = 0
  contains
    procedure :: init
    procedure :: update_em_heat, get_event_times
    procedure :: write_restart_data, read_restart_data
    procedure, nopass :: skip_restart_data
    procedure, private :: em_prop_differ
  end type

contains

  subroutine init(this, mesh, params, stat, errmsg)

    use em_properties

    class(microwave_heat_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(parameter_list), pointer :: plist
    real(r8) :: freq

    this%mesh => mesh
    this%params => params

    this%const_eps  = permittivity_is_const()
    this%const_epsi = permittivity_im_is_const()
    this%const_mu   = permeability_is_const()

    n = this%mesh%ncell
    allocate(this%q_eps(n), this%q_epsi(n), this%q_mu(n))
    !NB: this%q allocated later when defined

    plist => params%sublist('bc')
    call this%wg_port_fac%init(params, plist, stat, errmsg)
    if (stat /= 0) return

    call params%get('prop-change-threshold', this%prop_change_threshold, stat, errmsg, default=0.3_r8)
    if (stat /= 0) return
    if (this%prop_change_threshold <= 0.0_r8) then
      stat = -1
      errmsg = 'prop-change-threshold is <= 0.0'
      return
    end if

    !! Single fixed frequency for all time
    call params%get('frequency', freq, stat, errmsg)
    if (stat /= 0) return
    if (freq <= 0) then
      stat = 1
      errmsg = 'non-positive frequency'
      return
    end if
    this%omega = 8*atan(1.0_r8)*freq

    call params%get('graphics-output', this%graphics_output, stat, errmsg, default=.false.)
    if (stat /= 0) return

  end subroutine

  !! Return a list of times when the waveguide input powers (and hence EM heat
  !! source) change abruptly. This info can be used to control time stepping
  !! to hit those times exactly and reduce the step size and restart the HT ODE
  !! integrator in response to the discontinuous heat source.

  subroutine get_event_times(this, times)
    class(microwave_heat_solver), intent(in) :: this
    real(r8), allocatable, intent(out) :: times(:)
    times = this%wg_port_fac%times
  end subroutine

  !! This primary subroutine is called every heat transfer time step to update
  !! the EM heat source Q for the given time T and cell-based property arrays.
  !! STAT returns the status: 0 if Q is unchanged; 1 if Q is updated; and 2 if
  !! Q is updated and restart data should be written (invoked by the caller).
  !! Otherwise, if an error occurs, STAT returns a negative value and ERRMSG
  !! returns an explanatory message.

  ! NB: Ideally only time and temperature would be passed, and the EM properties
  ! evaluated internally. However, materials and properties are only known on
  ! the heat transfer mesh currently, and so they must be evaluated there and
  ! then mapped to the EM mesh to be passed here.

  subroutine update_em_heat(this, t, eps, epsi, mu, sigma, q, stat, errmsg)

    class(microwave_heat_solver), intent(inout) :: this
    real(r8), intent(in) :: t, eps(:), epsi(:), mu(:), sigma(:)
    real(r8), intent(inout) :: q(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: scf

    integer :: task
    enum, bind(c) ! task values
      enumerator :: KEEP=0, ZERO, SCALE, COMPUTE
    end enum

    ASSERT(size(q) == this%mesh%ncell_onp)
    ASSERT(size(eps) == this%mesh%ncell)
    ASSERT(size(epsi) == this%mesh%ncell)
    ASSERT(size(mu) == this%mesh%ncell)
    ASSERT(size(sigma) == this%mesh%ncell)

    !! Determine what needs to be done, if anything
    if (allocated(this%q)) then ! there is an existing Q that may be usable
      if (this%wg_port_fac%power_is_zero(t)) then
        task = merge(KEEP, ZERO, this%prev_q_is_zero)
      else if (this%em_prop_differ(eps, epsi, mu)) then
        task = COMPUTE
      else if (this%wg_port_fac%power_differs(this%q_data, t)) then
        task = merge(SCALE, COMPUTE, this%wg_port_fac%power_is_scaled(this%q_data, t, scf))
        if (task == SCALE .and. scf == this%prev_q_scf) task = KEEP
      else
        task = KEEP
      end if
    else ! there is no existing Q and it needs to be set
      allocate(this%q(this%mesh%ncell_onP))
      task = merge(ZERO, COMPUTE, this%wg_port_fac%power_is_zero(t))
    end if

    select case (task)
    case (KEEP)
      if (this%q_from_restart) then
        this%q_from_restart = .false.
        q = this%q
        stat = 2 ! q is updated and write restart data
      else
        stat = 0 ! do nothing
      end if

    case (ZERO)
      q = 0.0_r8
      this%prev_q_is_zero = .true.
      this%prev_q_scf = 0
      if (this%q_from_restart) then
        this%q_from_restart = .false.
        stat = 2 ! q is updated and write restart data
      else
        stat = 1 ! q is updated
      end if

      call tls_info('')
      call tls_info('MICROWAVE_HEAT_SOLVER: no power; EM heat set to zero')

    case (SCALE)
      q = scf**2 * this%q
      this%prev_q_scf = scf
      this%prev_q_is_zero = .false.
      if (this%q_from_restart) then
        this%q_from_restart = .false.
        stat = 2 ! q is updated and write restart data
      else
        stat = 1 ! q is updated
      end if

      call tls_info('')
      call tls_info('MICROWAVE_HEAT_SOLVER: power scaled; EM heat scaled accordingly')

    case (COMPUTE)

      call this%wg_port_fac%set_plist_power(t)

      call start_timer('simulation')
      comp: block

        use fdme_solver_type
        use string_utilities, only: i_to_c
        use truchas_env, only: output_dir

        type(fdme_solver) :: solver
        character(:), allocatable :: filepath

        call tls_info('')
        call tls_info('MICROWAVE_HEAT_SOLVER: computing the EM heat')

        !TODO: use previous solution as initial guess (residual correction)
        call solver%init(this%mesh, this%omega, eps, epsi, mu, sigma, this%params, stat, errmsg)
        if (stat /= 0) exit comp
        call solver%solve(stat, errmsg)
        if (stat /= 0) exit comp

        call solver%get_heat_source(q)
        this%prev_q_is_zero = .false.
        this%prev_q_scf = 0

        this%q = q
        this%q_data = this%wg_port_fac%power_data(t)
        this%q_eps = eps
        this%q_epsi = epsi
        this%q_mu = mu

        if (this%graphics_output) then
          this%sim_num = this%sim_num + 1
          filepath = trim(output_dir) // 'fdme-' // i_to_c(this%sim_num) // '.vtkhdf'
          call fdme_vtk_graphics(solver, filepath, this%mesh, q, stat, errmsg)
          if (stat /= 0) exit comp
        end if

      end block comp
      call stop_timer('simulation')

      if (stat /= 0) then ! solver failure
        stat = -1
      else ! success
        stat = 2 ! updated value and write restart data
        this%q_from_restart = .false.
      end if

    end select

  end subroutine

  !! This auxiliary function returns true if the EM properties (permittivity
  !! and permeability) have changed significantly since that last time the
  !! dielectric heat was computed; otherwise it returns false. The maximum
  !! relative difference in a property over the cells is taken as the measure
  !! of the change, and is considered significant if it exceeds a user-specified
  !! threshold.
  !!
  !! NB: Constant properties do not normally need to be checked. The exception
  !! is when the stored EM heat and its properties came from a restart file.
  !!
  !! NB: conductivity is not considered as it can be lumped with the imaginary
  !! part of the permittivity if desired, and is generally insignificant at
  !! microwave frequencies.

  logical function em_prop_differ(this, eps, epsi, mu) result(differ)

    class(microwave_heat_solver), intent(inout) :: this
    real(r8), intent(in) :: eps(:), epsi(:), mu(:)

    integer :: j
    real(r8) :: deps, depsi, dmu
    character(80) :: string

    differ = .false.

    if (.not.this%const_eps .or. this%check_restart_q_props) then
      deps = global_maxval(abs(eps-this%q_eps)/eps)
      differ = .true.
    else
      deps = 0.0_r8
    end if

    if (.not.this%const_mu .or. this%check_restart_q_props) then
      dmu = global_maxval(abs(mu-this%q_mu)/mu)
      differ = .true.
    else
      dmu = 0.0_r8
    end if

    !TODO: Possibly not a great way to define the difference; something with loss tangent instead?
    if (.not.this%const_epsi .or. this%check_restart_q_props) then
      depsi = 0.0_r8
      do j = 1, size(epsi)
        if (epsi(j) > 0.0_r8) depsi = max(depsi, abs(epsi(j)-this%q_epsi(j))/epsi(j))
      end do
      depsi = global_maxval(depsi)
      differ = .true.
    else
      depsi = 0.0_r8
    end if

    if (.not.differ) return ! nothing to check

    if (this%check_restart_q_props) then
      this%check_restart_q_props = .false.
      differ = (max(deps, depsi, dmu) > 0.0_r8)
    else
      differ = (max(deps, depsi, dmu) > this%prop_change_threshold)
    end if

    write(string,fmt='(a,3(a,es8.2))') 'MICROWAVE_HEAT_SOLVER: max rel delta: ', &
        'eps=', deps, ', epsi=', depsi, ', mu=', dmu
    call tls_info(string)

  end function em_prop_differ


  subroutine write_restart_data(this, sim)

    use truchas_h5_outfile, only: th5_sim_group
    use permutations

    class(microwave_heat_solver), intent(in) :: this
    type(th5_sim_group), intent(in) :: sim

    integer :: n
    integer, allocatable :: cell_perm(:)
    real(r8), allocatable :: array(:)

    !! NNC, 2/2017. Not entirely sure why the reordering for EPS and MU are
    !! necessary.  I believe if we wrote the cell map for the tet mesh, then
    !! post-processing tools could do it when needed, as is done for the main
    !! mesh, and we could dispense with the collation here and truly write in
    !! parallel.  FIXME

    n = global_sum(this%mesh%ncell_onP)
    allocate(array(merge(n,0,is_IOP)))

    !! Collate the cell permutation array.
    allocate(cell_perm(merge(n,0,is_IOP)))
    call gather(this%mesh%xcell(:this%mesh%ncell_onP), cell_perm)

    !! Collate the cell-based Q array on the tet mesh, and restore it to the external order.
    call gather(this%q(:this%mesh%ncell_onP), array)
    if (is_IOP) call reorder(array, cell_perm, forward=.true.)
    call sim%write_repl_data('Q', array)

    !! Collate the cell-based EPS array on the tet mesh, and restore it to the external order.
    call gather(this%q_eps(:this%mesh%ncell_onP), array)
    if (is_IOP) call reorder(array, cell_perm, forward=.true.)
    call sim%write_repl_data('EPS', array)

    !! Collate the cell-based EPSI array on the tet mesh, and restore it to the external order.
    call gather(this%q_epsi(:this%mesh%ncell_onP), array)
    if (is_IOP) call reorder(array, cell_perm, forward=.true.)
    call sim%write_repl_data('EPSI', array)

    !! Collate the cell-based MU array on the tet mesh, and restore it to the external order.
    call gather(this%q_mu(:this%mesh%ncell_onP), array)
    if (is_IOP) call reorder(array, cell_perm, forward=.true.)
    call sim%write_repl_data('MU', array)

    call sim%write_repl_data('QDATA', this%q_data)

  end subroutine write_restart_data

  !! Read the parameters specific to microwave heating segment of the restart
  !! file opened (and prepositioned) on UNIT, and initialize the Q_DATA, Q_EPS,
  !! AND Q_MU components. If any errors or incompatibilities with data read
  !! from the input file are encountered, execution is halted.

  subroutine read_restart_data(this, unit, version)

    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c

    class(microwave_heat_solver), intent(inout) :: this
    integer, intent(in) :: unit,  version

    integer :: n

    this%check_restart_q_props = .true.
    this%q_from_restart = .true.

    call read_var(unit, n, 'MICROWAVE_HEAT_SOLVER: error reading NQ restart record')
    if (n /= this%mesh%cell_imap%global_size) &
        call halt('MICROWAVE_HEAT_SOLVER: incompatible restart NQ value: ' // i_to_c(n))
    allocate(this%q(this%mesh%ncell_onp))
    call read_dist_array(unit, this%q, this%mesh%xcell(:this%mesh%ncell_onP), &
        'MICROWAVE_HEAT_SOLVER: error reading Q restart record')
    call this%mesh%cell_imap%gather_offp(this%q_eps)

    call read_var(unit, n, 'MICROWAVE_HEAT_SOLVER: error reading NEPS restart record')
    if (n /= this%mesh%cell_imap%global_size) &
        call halt('MICROWAVE_HEAT_SOLVER: incompatible restart NEPS value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_eps(:this%mesh%ncell_onP), &
        this%mesh%xcell(:this%mesh%ncell_onP), &
        'MICROWAVE_HEAT_SOLVER: error reading EPS restart record')
    call this%mesh%cell_imap%gather_offp(this%q_eps)

    call read_var(unit, n, 'MICROWAVE_HEAT_SOLVER: error reading NEPSI restart record')
    if (n /= this%mesh%cell_imap%global_size) &
        call halt('MICROWAVE_HEAT_SOLVER: incompatible restart NEPSI value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_epsi(:this%mesh%ncell_onP), &
        this%mesh%xcell(:this%mesh%ncell_onP), &
        'MICROWAVE_HEAT_SOLVER: error reading EPSI restart record')
    call this%mesh%cell_imap%gather_offp(this%q_epsi)

    call read_var(unit, n, 'MICROWAVE_HEAT_SOLVER: error reading NMU restart record')
    if (n /= this%mesh%cell_imap%global_size) &
        call halt('MICROWAVE_HEAT_SOLVER: incompatible restart NMU value: ' // i_to_c(n))
    call read_dist_array(unit, this%q_mu(:this%mesh%ncell_onP), &
        this%mesh%xcell(:this%mesh%ncell_onP), &
        'MICROWAVE_HEAT_SOLVER: error reading MU restart record')
    call this%mesh%cell_imap%gather_offp(this%q_mu)

    call read_var(unit, n, 'MICROWAVE_HEAT_SOLVER: error reading NQDATA restart record')
    allocate(this%q_data(n))
    call read_var(unit, this%q_data, 'MICROWAVE_HEAT_SOLVER: error reading QDATA restart record')

  end subroutine read_restart_data

  subroutine skip_restart_data(unit, version)
    use restart_utilities, only: skip_records
    integer, intent(in) :: unit, version
    call skip_records(unit, 10, 'MICROWAVE_HEAT_SOLVER: error skipping the restart data')
  end subroutine

  !FIXME: following procedure is replicated in induction_heat_solver_type

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

end module microwave_heat_solver_type
