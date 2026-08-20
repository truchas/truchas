!!
!! NS_2D_SOLVER_TYPE
!!
!! This module defines NS_2D_SOLVER, the isothermal incompressible
!! Navier--Stokes step algorithm.  It uses the common flow model, state, and
!! pressure projection, with an explicit first-order donor-cell treatment of
!! momentum transport in the predictor RHS.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ns_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parallel_communication, only: global_minval
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_momentum_solver_type
  use flow_2d_projection_solver_type
  use flow_2d_projection_update_type
  use flow_2d_ic_solver_type
  implicit none
  private

  type, public :: ns_2d_solver
    private
    type(flow_2d_model), pointer :: model => null()  ! unowned reference
    type(flow_2d_state), pointer :: state => null()  ! unowned reference
    type(flow_2d_momentum_solver) :: momentum_solver
    type(flow_2d_projection_solver), pointer :: projection_solver => null()
    type(flow_2d_projection_update) :: projection_update
    type(flow_2d_ic_solver), pointer :: ic_solver => null()
    real(r8), allocatable :: rhs(:,:), grad_p(:,:)
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: step
    procedure :: courant_time_step
    final :: delete
  end type

contains

  subroutine init(this, model, state, momentum_params, projection_params)
    class(ns_2d_solver), intent(out) :: this
    type(flow_2d_model), target, intent(in) :: model
    type(flow_2d_state), target, intent(inout) :: state
    type(parameter_list), target, intent(in) :: momentum_params, projection_params

    this%model => model
    this%state => state
    allocate(this%rhs(2, model%mesh%ncell_onP), this%grad_p(2, model%mesh%ncell), this%projection_solver, &
        this%ic_solver)
    call this%momentum_solver%init(model%momentum, momentum_params)
    call this%projection_solver%init(model%projection, projection_params)
    call this%projection_update%init(model%mesh, model%operators, model%projection, this%projection_solver, &
        model%body_acceleration)
    call this%ic_solver%init(model, momentum_params, projection_params)
  end subroutine


  subroutine delete(this)
    type(ns_2d_solver), intent(inout) :: this

    if (associated(this%projection_solver)) deallocate(this%projection_solver)
    if (associated(this%ic_solver)) deallocate(this%ic_solver)
  end subroutine


  !! Set STATE from an input velocity. The common initial-condition solver
  !! projects the velocity and computes an initial pressure with its temporary
  !! Stokes step, as mainline does when it omits initial momentum transport.
  subroutine set_initial_state(this, time, dt, velocity, stat)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    call this%ic_solver%solve(time, dt, velocity, this%state, stat)
  end subroutine


  !! Advance STATE from TIME to TIME + DT. Momentum transport is evaluated
  !! from the old cell and face velocity fields and added explicitly to the
  !! otherwise implicit mass-plus-viscous predictor.
  subroutine step(this, time, dt, stat, errmsg)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg

    integer :: c
    character(:), allocatable :: bc_errmsg

    ASSERT(dt > 0.0_r8)
    call this%model%compute_bc(time, dt, stat, bc_errmsg)
    if (stat /= 0) then
      if (present(errmsg)) errmsg = bc_errmsg
      return
    end if
    call this%model%pressure_gradient(this%state%p_cc, this%grad_p)
    call this%model%momentum%assemble(dt, this%model%density_c, this%model%viscosity_f, &
        this%model%bc, this%rhs)
    call this%model%momentum%add_advective_rhs(dt, this%model%density_c, this%state%vel_cc, &
        this%state%vel_fn, this%model%bc, this%rhs)
    do c = 1, size(this%rhs,2)
      this%rhs(:,c) = this%rhs(:,c) + this%model%density_c(c)*this%model%mesh%volume(c)*this%state%vel_cc(:,c) - &
          dt*this%model%mesh%volume(c)*this%grad_p(:,c)
    end do
    call this%momentum_solver%setup()
    call this%momentum_solver%solve(this%rhs, this%state%vel_cc(:,1:size(this%rhs,2)), stat)
    if (stat /= 0) return
    call this%model%mesh%cell_imap%gather_offp(this%state%vel_cc)
    call this%projection_update%correct(dt, this%model%inv_density_c, this%model%inv_density_f, &
        this%model%bc, this%state, stat)
  end subroutine


  !! Return the maximum step size satisfying the specified convective Courant
  !! number for the old face-normal velocity.
  function courant_time_step(this, courant_number) result(dt)
    class(ns_2d_solver), intent(in) :: this
    real(r8), intent(in) :: courant_number
    real(r8) :: dt

    integer :: f
    real(r8) :: dt_local

    ASSERT(courant_number > 0.0_r8 .and. courant_number <= 1.0_r8)
    dt_local = huge(1.0_r8)
    do f = 1, this%model%mesh%nface_onP
      if (this%state%vel_fn(f) == 0.0_r8) cycle
      dt_local = min(dt_local, this%model%operators%normal_distance(f)/abs(this%state%vel_fn(f)))
    end do
    dt = courant_number*global_minval(dt_local)
  end function

end module ns_2d_solver_type
