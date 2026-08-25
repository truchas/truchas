!!
!! HT_2D_IC_SOLVER_TYPE
!!
!! This module defines the initial-condition solver used by the 2D thermal
!! transport solver. It takes cell temperature initial conditions and computes
!! the full thermal vector state and initial time derivative needed by the
!! implicit integrator.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Face temperatures are solved to satisfy their algebraic residual equations.
!! Cell and face temperature derivatives are computed by finite differencing
!! an advanced state with similarly solved face temperatures.
!!

#include "f90_assert.fpp"

module ht_2d_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use ht_2d_model_type
  use ht_2d_vector_type
  use parameter_list_type
  use parallel_communication, only: global_dot_product
  use simulation_environment_type
  use simulation_log_type, only: LOG_DETAIL
  implicit none
  private

  type, public :: ht_2d_ic_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    type(ht_2d_model), pointer :: model => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_udot
  end type

contains

  subroutine init(this, model, params)
    class(ht_2d_ic_solver), intent(out) :: this
    type(ht_2d_model), intent(in), target :: model
    type(parameter_list), intent(inout), target :: params
    this%model => model
    this%mesh => model%mesh
    this%params => params
  end subroutine


  subroutine compute(this, env, t, temp, u, udot, stat, errmsg)

    class(ht_2d_ic_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t, temp(:)
    type(ht_2d_vector), intent(inout) :: u, udot
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    ASSERT(size(temp) >= this%mesh%ncell_onP)

    call u%setval(0.0_r8)
    u%tc(:this%mesh%ncell_onP) = temp(:this%mesh%ncell_onP)
    call this%mesh%cell_imap%gather_offp(u%tc)

    !! Averaging provides a good initial guess for the algebraic face solve.
    call average_to_faces(this%mesh, u%tc, u%tf)
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      u%tf(this%model%bc_dir%index) = this%model%bc_dir%value
    end if
    if (allocated(this%model%bc_inflow)) then
      call this%model%bc_inflow%compute(t)
      u%tf(this%model%bc_inflow%index) = this%model%bc_inflow%value
    end if

    call this%model%H_of_T%compute_value(u%tc, u%hc(:this%mesh%ncell_onP))

    call compute_face_temp(this, env, t, u, stat, errmsg)
    if (stat /= 0) return

    call this%compute_udot(env, t, u, udot, stat, errmsg)

  end subroutine compute

  !! Compute the initial time derivative by advancing enthalpy one small step,
  !! solving the algebraic variables at that advanced state, and differencing.

  subroutine compute_udot(this, env, t, u, udot, stat, errmsg)

    class(ht_2d_ic_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t
    type(ht_2d_vector), intent(inout) :: u, udot
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(ht_2d_vector) :: f, advanced
    real(r8) :: dt, Tmin, Tmax
    integer :: j

    stat = 0

    call this%params%get('dt', dt, stat, errmsg)
    if (stat /= 0) return
    if (dt <= 0.0_r8) then
      stat = 1
      errmsg = '"dt" must be > 0.0'
      return
    end if

    call udot%setval(0.0_r8)
    call f%init(u)
    call f%setval(0.0_r8)
    call this%model%residual(t, u, udot, f)
    udot%hc(:this%mesh%ncell_onP) = -f%tc(:this%mesh%ncell_onP) / &
        this%mesh%volume(:this%mesh%ncell_onP)

    call advanced%init(u)
    call advanced%copy(u)
    advanced%hc(:this%mesh%ncell_onP) = u%hc(:this%mesh%ncell_onP) + &
        dt * udot%hc(:this%mesh%ncell_onP)

    do j = 1, this%mesh%ncell_onP
      Tmin = u%tc(j) - 1.0_r8
      Tmax = u%tc(j) + 1.0_r8
      call this%model%T_of_H%compute(j, advanced%hc(j), Tmin, Tmax, advanced%tc(j))
    end do
    call this%mesh%cell_imap%gather_offp(advanced%tc)

    call compute_face_temp(this, env, t+dt, advanced, stat, errmsg)
    if (stat /= 0) return

    udot%tc(:this%mesh%ncell_onP) = (advanced%tc(:this%mesh%ncell_onP) - &
        u%tc(:this%mesh%ncell_onP)) / dt
    udot%tf(:this%mesh%nface_onP) = (advanced%tf(:this%mesh%nface_onP) - &
        u%tf(:this%mesh%nface_onP)) / dt
    call udot%gather_offp

  end subroutine compute_udot

  !! Solve the face-temperature algebraic equations using the A22 diffusion
  !! block. This belongs to the initial-condition solver rather than the
  !! time-integration model because it constructs a consistent initial state.
  !! NB: The current 2D face-temperature system is linear for fixed cell
  !! temperatures. A future nonlinear face system will require a nonlinear
  !! solve here, rather than this single CG solve.

  subroutine compute_face_temp(this, env, t, u, stat, errmsg)

    use hypre_hybrid_type
    use mfd_2d_diff_matrix_type

    class(ht_2d_ic_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t
    type(ht_2d_vector), intent(inout) :: u
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_2d_diff_matrix), target :: dm
    type(parameter_list), target :: solver_params
    type(hypre_hybrid) :: solver
    type(ht_2d_vector) :: udot, f
    real(r8), allocatable :: coef(:), z(:)
    real(r8) :: norm, rel_tol
    integer :: max_itr, num_itr, num_dscg_itr, num_pcg_itr
    character(80) :: msg

    call this%params%get('rel-tol', rel_tol, stat, errmsg)
    if (stat /= 0) return
    if (rel_tol <= 0.0_r8) then
      stat = 1
      errmsg = '"rel-tol" must be > 0.0'
      return
    end if
    call this%params%get('max-iter', max_itr, stat, errmsg)
    if (stat /= 0) return
    if (max_itr <= 0) then
      stat = 1
      errmsg = '"max-iter" must be > 0'
      return
    end if

    call udot%init(u)
    call udot%setval(0.0_r8)
    call f%init(u)
    call f%setval(0.0_r8)
    call this%model%residual(t, u, udot, f)
    norm = sqrt(global_dot_product(f%tf(:this%mesh%nface_onP), f%tf(:this%mesh%nface_onP)))

    if (env%simlog%is_enabled(LOG_DETAIL)) then
      write (msg,'(a,es10.3)') 'ht_2d_ic_solver%compute_face_temp: initial ||rface||_2 = ', norm
      call env%simlog%info(trim(msg), level=LOG_DETAIL)
    end if
    if (norm == 0.0_r8) return

    allocate(coef(this%mesh%ncell))
    call this%model%conductivity%compute_value(u%tc, coef(:this%mesh%ncell_onP))
    call this%mesh%cell_imap%gather_offp(coef)
    call dm%init(this%model%disc)
    call dm%compute(coef)
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if
    if (allocated(this%model%bc_inflow)) then
      call this%model%bc_inflow%compute(t)
      call dm%set_dir_faces(this%model%bc_inflow%index)
    end if
    if (allocated(this%model%bc_htc)) then
      call this%model%bc_htc%compute(t, u%tf)
      call dm%incr_face_diag(this%model%bc_htc%index, this%model%bc_htc%deriv)
    end if
    if (allocated(this%model%bc_rad)) then
      call this%model%bc_rad%compute(t, u%tf)
      call dm%incr_face_diag(this%model%bc_rad%index, this%model%bc_rad%deriv)
    end if

    call solver_params%set('krylov-method', 'cg')
    call solver_params%set('max-ds-iter', max_itr)
    call solver_params%set('max-amg-iter', max_itr)
    call solver_params%set('rel-tol', rel_tol)
    if (env%simlog%is_enabled(LOG_DETAIL)) then
      call solver_params%set('print-level', 1)
      call solver_params%set('logging-level', 1)
    else
      call solver_params%set('print-level', 0)
      call solver_params%set('logging-level', 0)
    end if
    call solver%init(dm%a22, solver_params)
    call solver%setup()

    allocate(z(this%mesh%nface_onP))
    z = 0.0_r8
    call solver%solve(f%tf(:this%mesh%nface_onP), z, stat)
    u%tf(:this%mesh%nface_onP) = u%tf(:this%mesh%nface_onP) - z

    if (env%simlog%is_enabled(LOG_DETAIL)) then
      call solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(msg,'(3(a,i0),a,es9.2)') 'solve: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      call env%simlog%info(trim(msg), level=LOG_DETAIL)
    end if

    if (stat /= 0) then
      call solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(msg,'(3(a,i0),a,es9.2)') 'failed to converge: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      errmsg = trim(msg)
    end if

  end subroutine compute_face_temp


  subroutine average_to_faces(mesh, ucell, uface)

    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: ucell(:)
    real(r8), intent(out) :: uface(:)

    integer :: j
    integer :: scale(size(uface))

    ASSERT(size(ucell) == mesh%ncell)
    ASSERT(size(uface) == mesh%nface)

    uface = 0.0_r8
    scale = 0
    do j = 1, mesh%ncell
      associate (cface => mesh%cface(mesh%cstart(j):mesh%cstart(j+1)-1))
        uface(cface) = uface(cface) + ucell(j)
        scale(cface) = scale(cface) + 1
      end associate
    end do
    call mesh%face_imap%gather_offp(uface)
    call mesh%face_imap%gather_offp(scale)

    where (scale == 0)
      uface = 0.0_r8
    elsewhere
      uface = uface / scale
    end where

  end subroutine average_to_faces

end module ht_2d_ic_solver_type
