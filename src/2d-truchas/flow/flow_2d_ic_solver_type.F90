!!
!! FLOW_2D_IC_SOLVER_TYPE
!!
!! This module defines FLOW_2D_IC_SOLVER, which repairs a supplied initial
!! cell velocity to satisfy velocity boundary conditions and discrete
!! continuity.  It then computes an initial pressure by an artificial Stokes
!! predictor/projection step while retaining the repaired velocity state.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_momentum_solver_type
  use flow_2d_projection_solver_type
  use flow_2d_projection_update_type
  use parallel_communication, only: global_maxval
  implicit none
  private

  type, public :: flow_2d_ic_solver
    private
    type(flow_2d_model), pointer :: model => null()  ! unowned reference
    type(flow_2d_momentum_solver) :: momentum_solver
    type(flow_2d_projection_solver), pointer :: projection_solver => null()
    type(flow_2d_projection_update) :: projection_update
    real(r8), allocatable :: rhs(:,:), velocity_cc(:,:), velocity_fn(:)
  contains
    procedure :: init
    procedure :: solve
    final :: delete
  end type

contains

  subroutine init(this, model, momentum_params, projection_params)
    class(flow_2d_ic_solver), intent(out) :: this
    type(flow_2d_model), target, intent(in) :: model
    type(parameter_list), target, intent(in) :: momentum_params, projection_params

    this%model => model
    allocate(this%projection_solver, this%rhs(2,model%mesh%ncell_onP), &
        this%velocity_cc(2,model%mesh%ncell), this%velocity_fn(model%mesh%nface))
    call this%momentum_solver%init(model%momentum, momentum_params)
    call this%projection_solver%init(model%projection, projection_params)
    call this%projection_update%init(model%mesh, model%operators, model%projection, this%projection_solver)
  end subroutine


  !! Set STATE from a supplied cell-centered velocity.  The accepted velocity
  !! is its discrete projection onto the velocity BCs and continuity
  !! constraint.  The temporary artificial step supplies initial pressure but
  !! its velocity update is rejected.
  subroutine solve(this, time, dt, velocity, state, stat)
    class(flow_2d_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt, velocity(:,:)
    type(flow_2d_state), intent(inout) :: state
    integer, intent(out) :: stat

    integer :: c

    ASSERT(dt > 0.0_r8)
    ASSERT(size(velocity,1) == 2)
    ASSERT(size(velocity,2) >= this%model%mesh%ncell_onP)
    stat = 0
    call state%set_zero()
    state%vel_cc(:,1:this%model%mesh%ncell_onP) = velocity(:,1:this%model%mesh%ncell_onP)
    call this%model%bc%compute_initial(time)
    call this%projection_update%project_velocity(dt, this%model%inv_density_c, this%model%inv_density_f, &
        this%model%bc, state, stat)
    if (stat /= 0) return
    this%velocity_cc = state%vel_cc
    this%velocity_fn = state%vel_fn

    !! Compute pressure from a temporary Stokes step.  The repaired velocity
    !! remains the initial datum after this procedure returns.
    state%p_cc = 0.0_r8
    call this%model%momentum%assemble(dt, this%model%density_c, this%model%viscosity_f, &
        this%model%bc, this%rhs)
    do c = 1, this%model%mesh%ncell_onP
      this%rhs(:,c) = this%rhs(:,c) + this%model%density_c(c)*this%model%mesh%volume(c)*state%vel_cc(:,c)
    end do
    state%vel_cc = 0.0_r8
    if (global_maxval(maxval(abs(this%rhs))) > 0.0_r8) then
      call this%momentum_solver%setup()
      call this%momentum_solver%solve(this%rhs, state%vel_cc(:,1:this%model%mesh%ncell_onP), stat)
      if (stat /= 0) return
    end if
    call this%model%mesh%cell_imap%gather_offp(state%vel_cc)
    call this%projection_update%correct(dt, this%model%inv_density_c, this%model%inv_density_f, &
        this%model%bc, state, stat)
    if (stat /= 0) return
    state%vel_cc = this%velocity_cc
    state%vel_fn = this%velocity_fn
  end subroutine


  subroutine delete(this)
    type(flow_2d_ic_solver), intent(inout) :: this

    if (associated(this%projection_solver)) deallocate(this%projection_solver)
  end subroutine

end module flow_2d_ic_solver_type
