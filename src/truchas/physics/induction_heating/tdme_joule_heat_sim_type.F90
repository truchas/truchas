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
  use ih_source_factory_type
  use scalar_func_class
  use EM_graphics_output  !TODO: modernize
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

  subroutine init(this, mesh, freq, eps, mu, sigma, bc_fac, params, stat, errmsg)

    use em_bc_factory_type
    use parameter_list_type
    use bndry_func1_class

    class(tdme_joule_heat_sim), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: freq, eps(:), mu(:), sigma(:)
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(bndry_func1), allocatable :: ebc, hbc
    real(r8) :: dt, eps_scf, sigma_scf

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

    call bc_fac%alloc_nxE_bc(ebc, stat, errmsg)
    if (stat /= 0) return !TODO: augment errmsg?
    if (allocated(ebc)) then
      call bc_fac%alloc_nxH_bc(hbc, stat, errmsg, omit_edge_list=ebc%index, scale_factor=(1.0_r8/freq))
    else
      call bc_fac%alloc_nxH_bc(hbc, stat, errmsg, scale_factor=(1.0_r8/freq))
    end if
    if (stat /= 0) return !TODO: augment errmsg?

    !TODO: ensure the BC cover the entire boundary (all boundary edges)
    !NB: nxH = 0 is the natural (i.e., do nothing) BC; could be the default?

    call params%get('steps-per-cycle', this%steps_per_cycle, stat=stat, errmsg=errmsg, default=20)
    if (stat /= 0) return
    if (this%steps_per_cycle < 1) then
      stat = 1
      errmsg = 'steps-per-cycle must be > 0'
      return
    end if

    call params%get('steady-state-tol', this%ss_tol, stat=stat, errmsg=errmsg, default=0.01_r8)
    if (stat /= 0) return
    if (this%ss_tol <= 0.0_r8) then
      stat = 1
      errmsg = 'steady-state-tol must be > 0.0'
      return
    end if

    call params%get('max-source-cycles', this%max_cycles, stat=stat, errmsg=errmsg, default=5)
    if (this%max_cycles < 1) then
      stat = 1
      errmsg = 'max-source-cycles must be > 0'
      return
    end if

    call params%get('graphics-output', this%graphics_output, stat=stat, errmsg=errmsg, default=.false.)
    if (stat /= 0) return
    if (this%graphics_output) call export_mesh(mesh, eps, mu, sigma)

    !! Create and initialize the time-discretized model and solver.
    dt = 1.0_r8 / this%steps_per_cycle
    allocate(this%model)
    call this%model%init(this%mesh, eps_scf*eps, mu, sigma_scf*sigma, dt, ebc, hbc)
    call this%solver%init(this%mesh, this%model, params, stat, errmsg)
    if (stat /= 0) return

  end subroutine init

  subroutine compute(this, q_avg, stat, errmsg)

    use parallel_communication, only: global_maxval, global_dot_product
    use truchas_logging_services

    class(tdme_joule_heat_sim), intent(inout) :: this
    real(r8), intent(out) :: q_avg(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    real(r8) :: t, q_avg_max, q_tot
    real(r8), allocatable :: efield(:), bfield(:), q(:), q_avg_last(:)
    character(len=256) :: string

    ASSERT(size(q_avg) == this%mesh%ncell)

    allocate(efield(this%mesh%nedge), bfield(this%mesh%nface), q(this%mesh%ncell))
    !TODO: make persistent?

    t = 0.0_r8
    efield = 0.0_r8
    bfield = 0.0_r8
    call this%solver%set_initial_state(t, efield, bfield)
    !NB: scaling client-provided fields would normally be necessary:
    !call this%solver%set_initial_state(t, efield/this%escf, bfield/this%bscf)
    q = 0.0_r8  ! assuming efield == 0

    if (this%graphics_output) then
      call initialize_field_output
      call export_fields(this%mesh, t, this%escf*efield, this%bscf*bfield, q)
    end if

    do n = 1, this%max_cycles
      !! Time step through one source cycle while integrating the joule heat
      if (this%graphics_output .and. n > 1) call initialize_field_output
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
        if (this%graphics_output) call export_fields(this%mesh, t, this%escf*efield, this%bscf*bfield, q)
      end do

      !! Time-averaged Joule power density over the last source cycle.
      q_avg = q_avg / this%steps_per_cycle
      q_avg_max = global_maxval(q_avg)
      if (this%graphics_output) call finalize_field_output(q_avg)

      !TODO: make output subject to verbosity level
      q_tot = global_dot_product(q_avg(:this%mesh%ncell_onP), abs(this%mesh%volume(:this%mesh%ncell_onP)))
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

end module tdme_joule_heat_sim_type
