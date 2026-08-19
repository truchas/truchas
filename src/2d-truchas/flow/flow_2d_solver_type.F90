!!
!! FLOW_2D_SOLVER_TYPE
!!
!! This module defines FLOW_2D_SOLVER, the first isothermal incompressible
!! flow-step algorithm.  It solves the implicit momentum predictor and then
!! applies an incremental pressure correction.  The mesh-associated model and
!! solution state remain owned by the caller.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_momentum_solver_type
  use flow_2d_projection_solver_type
  use flow_2d_projection_update_type
  use flow_2d_ic_solver_type
  implicit none
  private

  type, public :: flow_2d_solver
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
    final :: delete
  end type

contains

  subroutine init(this, model, state, momentum_params, projection_params)
    class(flow_2d_solver), intent(out) :: this
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


  subroutine set_initial_state(this, time, dt, velocity, stat)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    call this%ic_solver%solve(time, dt, velocity, this%state, stat)
  end subroutine


  subroutine delete(this)
    type(flow_2d_solver), intent(inout) :: this

    if (associated(this%projection_solver)) deallocate(this%projection_solver)
    if (associated(this%ic_solver)) deallocate(this%ic_solver)
  end subroutine


  !! Advance STATE from TIME to TIME + DT.  The first implementation has an
  !! implicit Stokes predictor with body acceleration but no advective momentum
  !! term.
  subroutine step(this, time, dt, stat)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt
    integer, intent(out) :: stat

    integer :: c

    ASSERT(dt > 0.0_r8)
    call this%model%compute_bc(time, dt)
    call this%model%pressure_gradient(this%state%p_cc, this%grad_p)
    call this%model%momentum%assemble(dt, this%model%density_c, this%model%viscosity_f, &
        this%model%bc, this%rhs)
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

end module flow_2d_solver_type
