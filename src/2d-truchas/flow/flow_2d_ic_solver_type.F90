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
  use simulation_environment_type
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
    real(r8), allocatable :: rhs(:,:), grad_p(:,:), velocity_cc(:,:), velocity_fn(:)
  contains
    procedure :: init
    procedure :: solve
    final :: delete
  end type

contains

  subroutine init(this, model, momentum_params, projection_params, stat, errmsg)
    class(flow_2d_ic_solver), intent(out) :: this
    type(flow_2d_model), target, intent(in) :: model
    type(parameter_list), target, intent(in), optional :: momentum_params
    type(parameter_list), target, intent(in) :: projection_params
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    this%model => model
    allocate(this%projection_solver, this%rhs(2,model%mesh%ncell_onP), this%grad_p(2,model%mesh%ncell), &
        this%velocity_cc(2,model%mesh%ncell), this%velocity_fn(model%mesh%nface))
    if (present(momentum_params)) then
      call this%momentum_solver%init(model%momentum, momentum_params, stat, errmsg)
      if (present(stat)) then
        if (stat /= 0) return
      end if
    end if
    call this%projection_solver%init(model%projection, projection_params, stat, errmsg)
    if (present(stat)) then
      if (stat /= 0) return
    end if
    call this%projection_update%init(model%mesh, model%operators, model%projection, this%projection_solver, &
        model%body_acceleration)
  end subroutine


  !! Set STATE from a supplied cell-centered velocity.  The accepted velocity
  !! is its discrete projection onto the velocity BCs and continuity
  !! constraint.  The temporary artificial step supplies initial pressure but
  !! its velocity update is rejected.
  subroutine solve(this, env, time, dt, velocity, state, stat)
    class(flow_2d_ic_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: time, dt, velocity(:,:)
    type(flow_2d_state), intent(inout) :: state
    integer, intent(out) :: stat

    integer :: c, num_itr, num_dscg_itr, num_pcg_itr
    real(r8) :: rel_res_norm
    character(:), allocatable :: errmsg
    logical :: projection_solved

    ASSERT(dt > 0.0_r8)
    ASSERT(size(velocity,1) == 2)
    ASSERT(size(velocity,2) >= this%model%mesh%ncell_onP)
    stat = 0
    call state%set_zero()
    state%vel_cc(:,1:this%model%mesh%ncell_onP) = velocity(:,1:this%model%mesh%ncell_onP)
    call this%model%bc%compute_initial(time)
    call this%model%bc%check_velocity_flux(stat, errmsg)
    if (stat /= 0) return
    !! Start the temporary solve from zero pressure.  As in mainline, the
    !! predictor/projection sequence computes the initial pressure using the
    !! same pressure-gradient and body-force operators used during stepping.
    state%p_cc = 0.0_r8
    call this%projection_update%project_velocity(dt, this%model%matl_props%inv_density_c, &
        this%model%matl_props%inv_density_f, &
        this%model%bc, state, stat, projection_solved)
    if (projection_solved) then
      call this%projection_solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
      call write_solver_metrics(env, 'flow.initial.velocity-projection', num_itr, num_dscg_itr, num_pcg_itr, &
          rel_res_norm, stat)
    else
      call env%simlog%info('flow.initial.velocity-projection method=none reason=zero-rhs status=skipped')
    end if
    if (stat /= 0) return
    this%velocity_cc = state%vel_cc
    this%velocity_fn = state%vel_fn

    !! Compute pressure from a temporary Stokes step.  The repaired velocity
    !! remains the initial datum after this procedure returns.
    !! The mainline initial predictor starts with no committed pressure
    !! gradient.  Body force enters through the projection below, using the
    !! same gravity-head operator as an ordinary step.
    this%grad_p = 0.0_r8
    call this%model%assemble_momentum(dt, this%rhs)
    do c = 1, this%model%mesh%ncell_onP
      this%rhs(:,c) = this%rhs(:,c) + this%model%matl_props%density_c(c)*this%model%mesh%volume(c)*state%vel_cc(:,c)
      this%rhs(:,c) = this%rhs(:,c) - dt*this%model%mesh%volume(c)*this%grad_p(:,c)
    end do
    state%vel_cc = 0.0_r8
    if (this%model%inviscid) then
      call this%model%momentum%solve_inviscid(this%model%matl_props%density_c, this%rhs, &
          state%vel_cc(:,1:this%model%mesh%ncell_onP))
      call env%simlog%info('flow.initial.momentum method=inviscid-direct status=ok')
    else if (global_maxval(maxval(abs(this%rhs))) > 0.0_r8) then
      call this%momentum_solver%setup()
      call this%momentum_solver%solve(this%rhs, state%vel_cc(:,1:this%model%mesh%ncell_onP), stat)
      call this%momentum_solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
      call write_solver_metrics(env, 'flow.initial.momentum', num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm, stat)
      if (stat /= 0) return
    else
      call env%simlog%info('flow.initial.momentum method=none reason=zero-rhs status=skipped')
    end if
    call this%model%mesh%cell_imap%gather_offp(state%vel_cc)
    call this%projection_update%correct(dt, this%model%matl_props%inv_density_c, &
        this%model%matl_props%inv_density_f, this%model%matl_props%density_delta_c, this%model%bc, state, stat, &
        initial=.true., solved=projection_solved)
    if (projection_solved) then
      call this%projection_solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
      call write_solver_metrics(env, 'flow.initial.pressure-projection', num_itr, num_dscg_itr, num_pcg_itr, &
          rel_res_norm, stat)
    else
      call env%simlog%info('flow.initial.pressure-projection method=none reason=zero-rhs status=skipped')
    end if
    if (stat /= 0) return
    state%vel_cc = this%velocity_cc
    state%vel_fn = this%velocity_fn
  end subroutine


  subroutine write_solver_metrics(env, name, num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm, stat)
    type(simulation_environment), intent(in) :: env
    character(*), intent(in) :: name
    integer, intent(in) :: num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8), intent(in) :: rel_res_norm

    character(256) :: line
    character(7) :: status

    if (stat == 0) then
      status = 'ok'
    else
      status = 'failed'
    end if
    write(line,'(a,i0,a,i0,a,i0,a,es0.5,a,a)') trim(name)//' iterations=', num_itr, &
        ' dscg=', num_dscg_itr, ' amg=', num_pcg_itr, ' rel_residual=', rel_res_norm, ' status=', trim(status)
    call env%simlog%info(trim(line))
  end subroutine


  subroutine delete(this)
    type(flow_2d_ic_solver), intent(inout) :: this

    if (associated(this%projection_solver)) deallocate(this%projection_solver)
  end subroutine

end module flow_2d_ic_solver_type
