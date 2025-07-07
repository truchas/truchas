!!
!! TDME_JOULE_HEAT_SIM_TYPE
!!
!! This module defines a derived type that encapsulates the computation of
!! time-averaged Joule heat by integrating the time-domain Maxwell equations
!! over several cycles of the external magnetic field forcing until a periodic
!! steady state is reached. The computed Joule heat would typically be used
!! as the source term in an induction heating simulation. Note that time in
!! the Joule heat simulation is a "fast time" regarded as unfolding in an
!! instant of heat transfer's "slow time", and justifies averaging the fast
!! time Joule heat over a cycle of its periodic steady state.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Refactored February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Parameter list input. The PARAMS argument to INIT has the following form:
!!
!!  {
!!    "steps-per-cycle": INTEGER
!!    "ss-stopping-tolerance": FLOAT
!!    "maximum-source-cycles": INTEGER
!!    "graphics-output": BOOLEAN (default false)
!!    "maximum-cg-iterations": INTEGER
!!    "cg-stopping-tolerance": FLOAT
!!    "output-level": INTEGER
!!  }
!!


#include "f90_assert.fpp"

module tdme_joule_heat_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use tdme_model_type
  use tdme_solver_type
  use scalar_func_class
  use vtkhdf_file_type
  implicit none
  private

  type, public :: tdme_joule_heat_sim
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(tdme_model), pointer :: model => null()  ! owned
    type(tdme_solver) :: solver
    real(r8) :: escf, bscf, qscf
    integer :: steps_per_cycle, max_cycles
    real(r8) :: ss_tol
    logical :: graphics_output
    type(vtkhdf_file) :: viz_file
  contains
    procedure :: init
    procedure :: compute
    final :: tdme_joule_heat_delete
  end type

contains

  subroutine tdme_joule_heat_delete(this)
    type(tdme_joule_heat_sim), intent(inout) :: this
    if (associated(this%model)) deallocate(this%model)
  end subroutine

  subroutine init(this, mesh, freq, eps, mu, sigma, params, stat, errmsg)

    use em_bc_factory_type
    use parameter_list_type
    use bndry_func1_class

    class(tdme_joule_heat_sim), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: freq, eps(:), mu(:), sigma(:)
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(em_bc_factory) :: bc_fac
    class(bndry_func1), allocatable :: ebc, hbc
    real(r8) :: dt, eps_scf, sigma_scf, c_ratio, omega
    real(r8), allocatable :: model_eps(:), model_sigma(:)
    character(:), allocatable :: filename
    type(parameter_list), pointer :: plist
    logical :: flag

    ASSERT(size(eps) == mesh%ncell)
    ASSERT(size(mu)  == mesh%ncell)
    ASSERT(size(sigma) == mesh%ncell)
    !TODO? use property mesh function for eps, mu, sigma, like HT instead of bare arrays?

    this%mesh => mesh

    !! Scale time for a unit forcing period: t' = freq * t. This is handled
    !! as a change of time units which is reflected in a change of units of
    !! the equation coefficients and variables. Note that mu has no time unit.
    eps_scf = freq**2
    sigma_scf = freq
    ! Physical var = scale factor * computational var
    this%escf = freq**2
    this%bscf = freq
    this%qscf = freq**3

    !NB: scaling of nxE BC data moved inside the model.
    !FIXME: Get rid of the scaling entirely as it creates sticky problems
    !FIXME: with nxH BC handling.

    omega = 8*atan(1.0_r8)*freq
    call params%get('use-legacy-bc', flag, stat, errmsg, default=.false.)
    if (stat /= 0) return
    plist => params%sublist('bc')
    call bc_fac%init(this%mesh, omega, plist, use_legacy_bc=flag)

    call bc_fac%alloc_nxE_bc(ebc, stat, errmsg)
    if (stat /= 0) return !TODO: augment errmsg?
    if (allocated(ebc)) then
      call bc_fac%alloc_nxH_bc(hbc, stat, errmsg, omit_edge_list=ebc%index)
    else
      call bc_fac%alloc_nxH_bc(hbc, stat, errmsg)
    end if
    if (stat /= 0) return !TODO: augment errmsg?

    !TODO: ensure the BC cover the entire boundary (all boundary edges)
    !NB: nxH = 0 is the natural (i.e., do nothing) BC; could be the default?

    call params%get('steps-per-cycle', this%steps_per_cycle, stat, errmsg, default=20)
    if (stat /= 0) return
    if (this%steps_per_cycle < 1) then
      stat = 1
      errmsg = 'steps-per-cycle must be > 0'
      return
    end if

    call params%get('steady-state-tol', this%ss_tol, stat, errmsg, default=0.01_r8)
    if (stat /= 0) return
    if (this%ss_tol <= 0.0_r8) then
      stat = 1
      errmsg = 'steady-state-tol must be > 0.0'
      return
    end if

    call params%get('max-source-cycles', this%max_cycles, stat, errmsg, default=5)
    if (this%max_cycles < 1) then
      stat = 1
      errmsg = 'max-source-cycles must be > 0'
      return
    end if

    call params%get('graphics-output', this%graphics_output, stat, errmsg, default=.false.)
    if (stat /= 0) return
    if (this%graphics_output) then
      call params%get('graphics-file', filename, stat, errmsg)
      if (stat /= 0) return
      call graphics_output_init(this, filename, stat, errmsg)
      if (stat /= 0) return
      call export_scalar_cell_field(this, eps, 'permittivity')
      call export_scalar_cell_field(this, mu, 'permeability')
      call export_scalar_cell_field(this, sigma, 'conductivity')
    end if

    model_eps = eps_scf*eps
    if (params%is_parameter('c-ratio')) then ! apply a numerical regularization in void cells
      call params%get('c-ratio', c_ratio, stat, errmsg)
      if (stat /= 0) return
      if (c_ratio <= 0.0_r8 .or. c_ratio >= 1.0_r8) then
        stat = 1
        errmsg = 'c-ratio must be > 0 and < 1'
        return
      end if
      where (sigma == 0.0_r8) model_eps = model_eps / c_ratio**2
    end if

    model_sigma = sigma_scf*sigma

    !! Create and initialize the time-discretized model and solver.
    dt = 1.0_r8 / this%steps_per_cycle
    allocate(this%model)
    call this%model%init(this%mesh, model_eps, mu, model_sigma, dt, ebc, hbc, this%bscf)
    call this%solver%init(this%mesh, this%model, params, stat, errmsg)
    if (stat /= 0) return

  end subroutine init

  subroutine compute(this, q_avg, stat, errmsg)

    use parallel_communication, only: global_maxval, global_dot_product
    use truchas_logging_services
    use string_utilities, only: i_to_c

    class(tdme_joule_heat_sim), intent(inout) :: this
    real(r8), intent(out) :: q_avg(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    real(r8) :: t, q_avg_max, q_tot
    real(r8), allocatable :: efield(:), bfield(:), q(:), q_avg_last(:)
    character(len=256) :: string

    ASSERT(size(q_avg) == this%mesh%ncell_onP)

    allocate(efield(this%mesh%nedge), bfield(this%mesh%nface))
    allocate(q(this%mesh%ncell_onP))
    !TODO: make persistent?

    t = 0.0_r8
    efield = 0.0_r8
    bfield = 0.0_r8
    call this%solver%set_initial_state(t, efield, bfield)
    !NB: scaling client-provided fields would normally be necessary:
    !call this%solver%set_initial_state(t, efield/this%escf, bfield/this%bscf)
    q = 0.0_r8  ! assuming efield == 0

    if (this%graphics_output) call export_fields(this, t, this%escf*efield, this%bscf*bfield, q)

    do n = 1, this%max_cycles
      !! Time step through one source cycle while integrating the joule heat
      q_avg = 0.0_r8  ! use trapezoid rule for the time average
      do j = 1, this%steps_per_cycle
        q_avg = q_avg + 0.5_r8 * q
        call this%solver%step(t, efield, bfield, stat, errmsg)
        if (stat /= 0) then
          stat = -1
          errmsg = 'time step failure: ' // errmsg
          return
        end if
        call this%model%compute_joule_heat(efield, q)
        q = this%qscf * q
        q_avg = q_avg + 0.5_r8 * q
        if (this%graphics_output) call export_fields(this, t, this%escf*efield, this%bscf*bfield, q)
      end do

      !! Time-averaged Joule power density over the last source cycle.
      q_avg = q_avg / this%steps_per_cycle
      q_avg_max = global_maxval(q_avg)
      if (this%graphics_output) call export_scalar_cell_field(this, q_avg, 'Avg_Joule-'//i_to_c(n))

      !TODO: make output subject to verbosity level
      q_tot = global_dot_product(q_avg, abs(this%mesh%volume(:this%mesh%ncell_onP)))
      write(string,fmt='(t4,a,i4,2(a,es11.4))') &
          'Source cycle', n, ': |Q|_max=', q_avg_max, ', Q_total=', q_tot
      call TLS_info(trim(string))

      if (n > 1) then ! check for convergence to a steady-state periodic solution
        if (global_maxval(abs(q_avg-q_avg_last)) < this%ss_tol*q_avg_max) return
      end if
      q_avg_last = q_avg
    end do

    !! Not yet converged; return Q_AVG from last cycle
    stat = 1
    errmsg = 'not converged to steady-state'

  end subroutine compute

  !! Create the VTKHDF format graphics file with the given file name, write the
  !! mesh, and register the time-dependent fields that will be written.

  subroutine graphics_output_init(this, filename, stat, errmsg)

    use parallel_communication, only: is_IOP, broadcast

    class(tdme_joule_heat_sim), intent(inout) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: vec_mold(3,0), sca_mold(0)

    if (is_IOP) call this%viz_file%create(filename, stat, errmsg)
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast(errmsg)
      return
    end if

    call export_mesh(this)

    if (is_IOP) call this%viz_file%register_temporal_cell_dataset('E-field', vec_mold, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call this%viz_file%register_temporal_cell_dataset('B-field', vec_mold, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call this%viz_file%register_temporal_cell_dataset('Joule', sca_mold, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

  end subroutine graphics_output_init

  !! Write the mesh to the VTKHDF graphics file.

  subroutine export_mesh(this)

    use,intrinsic :: iso_fortran_env, only: int8
    use parallel_communication, only: is_IOP, broadcast

    class(tdme_joule_heat_sim), intent(inout) :: this

    integer, allocatable, target :: cnode(:,:)
    integer, allocatable :: xcnode(:)
    integer(int8), allocatable :: types(:)
    real(r8), allocatable :: x(:,:)
    integer, pointer :: connectivity(:)
    integer :: j, stat
    character(:), allocatable :: errmsg

    !! Collate the mesh data structure onto the IO process
    call this%mesh%get_global_cnode_array(cnode)
    call this%mesh%get_global_x_array(x)

    if (is_IOP) then
      xcnode = [(1+4*j, j=0, size(cnode,dim=2))]
      connectivity(1:size(cnode)) => cnode ! flattened view
      types = spread(VTK_TETRA, dim=1, ncopies=size(cnode,dim=2))
      call this%viz_file%write_mesh(x, connectivity, xcnode, types, stat, errmsg)
    end if
    call broadcast(stat)
    INSIST(stat == 0)

  end subroutine export_mesh

  !! Write the given cell-based scalar field to the VTKHDF file as a
  !! time-independent dataset with the given name.

  subroutine export_scalar_cell_field(this, field, name)
    use parallel_communication, only: is_IOP, broadcast, gather
    class(tdme_joule_heat_sim), intent(in) :: this
    real(r8), intent(in) :: field(:)
    character(*), intent(in) :: name
    integer :: stat
    character(:), allocatable :: errmsg
    real(r8), allocatable :: g_field(:)
    allocate(g_field(merge(this%mesh%cell_imap%global_size, 0, is_IOP)))
    call gather(field(:this%mesh%ncell_onP), g_field)
    if (is_IOP) call this%viz_file%write_cell_dataset(name, g_field, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)
  end subroutine

  !! Write the given time-dependent EM fields to the VTKHDF graphics file.

  subroutine export_fields(this, t, efield, bfield, qfield)

    use parallel_communication, only: is_IOP, broadcast, gather
    use mimetic_discretization, only: w1_vector_on_cells, w2_vector_on_cells

    class(tdme_joule_heat_sim), intent(inout) :: this
    real(r8), intent(in) :: t, efield(:), bfield(:), qfield(:)

    real(r8) :: v(3,this%mesh%ncell)
    real(r8), allocatable :: g_v(:,:), g_s(:)
    integer :: stat
    character(:), allocatable :: errmsg

    if (is_IOP) call this%viz_file%write_time_step(t)

    allocate(g_v(3,merge(this%mesh%cell_imap%global_size,0,is_iop)))
    allocate(g_s(merge(this%mesh%cell_imap%global_size,0,is_iop)))

    !! Interpolate cell average E-field from the primitive E-field edge circulations.
    v = w1_vector_on_cells(this%mesh, efield)
    call gather(v(:,:this%mesh%ncell_onP), g_v)
    if (is_IOP) call this%viz_file%write_temporal_cell_dataset('E-field', g_v, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    !! Interpolate cell average B-field from the primitive B-field face fluxes.
    v = w2_vector_on_cells(this%mesh, bfield)
    call gather(v(:,:this%mesh%ncell_onP), g_v)
    if (is_IOP) call this%viz_file%write_temporal_cell_dataset('B-field', g_v, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    call gather(qfield(:this%mesh%ncell_onP), g_s)
    if (is_IOP) call this%viz_file%write_temporal_cell_dataset('Joule', g_s, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

  end subroutine

end module tdme_joule_heat_sim_type
