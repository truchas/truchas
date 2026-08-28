!!
!! FLOW_2D_SOLVER_TYPE
!!
!! This module defines FLOW_2D_SOLVER, the first isothermal incompressible
!! flow-step algorithm.  It solves the implicit momentum predictor and then
!! applies an incremental pressure correction.  It owns the evolving flow
!! state and provides views of that state to its caller.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
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
    type(flow_2d_state) :: state
    type(flow_2d_momentum_solver) :: momentum_solver
    type(flow_2d_projection_solver), pointer :: projection_solver => null()
    type(flow_2d_projection_update) :: projection_update
    type(flow_2d_ic_solver), pointer :: ic_solver => null()
    real(r8), allocatable :: rhs(:,:), grad_p(:,:)
    integer(int64) :: nstep = 0_int64
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    final :: delete
  end type

contains

  subroutine init(this, env, model, momentum_params, projection_params)
    class(flow_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(in) :: model
    type(parameter_list), target, intent(in), optional :: momentum_params
    type(parameter_list), target, intent(in) :: projection_params

    this%model => model
    call this%state%init(model%mesh)
    allocate(this%rhs(2, model%mesh%ncell_onP), this%grad_p(2, model%mesh%ncell), this%projection_solver, &
        this%ic_solver)
    if (present(momentum_params)) call this%momentum_solver%init(model%momentum, momentum_params)
    call this%projection_solver%init(model%projection, projection_params)
    call this%projection_update%init(model%mesh, model%operators, model%projection, this%projection_solver, &
        model%body_acceleration)
    if (present(momentum_params)) then
      call this%ic_solver%init(model, momentum_params, projection_params)
    else
      call this%ic_solver%init(model, projection_params=projection_params)
    end if
  end subroutine


  subroutine set_initial_state(this, time, dt, velocity, stat)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    call this%ic_solver%solve(time, dt, velocity, this%state, stat)
    if (stat == 0) this%nstep = 0_int64
  end subroutine


  !! Return no-copy views of the cell-centered pressure and velocity.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(flow_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    pressure => this%state%p_cc
    velocity => this%state%vel_cc
  end subroutine


  !! Return a no-copy view of the face-normal velocity.
  subroutine get_face_velocity(this, velocity)
    class(flow_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    velocity => this%state%vel_fn
  end subroutine


  subroutine delete(this)
    type(flow_2d_solver), intent(inout) :: this

    if (associated(this%projection_solver)) deallocate(this%projection_solver)
    if (associated(this%ic_solver)) deallocate(this%ic_solver)
  end subroutine


  !! Advance STATE from T_N to T_NP1.  The time step is derived from the two
  !! endpoint times so callers retain exact target times. The first
  !! implementation has an implicit Stokes predictor with body acceleration
  !! but no advective momentum term.
  subroutine step(this, t_n, t_np1, stat, errmsg)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg

    integer :: c
    real(r8) :: dt
    character(:), allocatable :: bc_errmsg

    dt = t_np1 - t_n
    ASSERT(dt > 0.0_r8)
    call this%model%compute_bc(t_n, dt, stat, bc_errmsg)
    if (stat /= 0) then
      if (present(errmsg)) errmsg = bc_errmsg
      return
    end if
    call this%model%pressure_gradient(this%state%p_cc, this%grad_p)
    call this%model%assemble_momentum(dt, this%rhs)
    do c = 1, size(this%rhs,2)
      this%rhs(:,c) = this%rhs(:,c) + this%model%matl_props%density_c(c)*this%model%mesh%volume(c)*this%state%vel_cc(:,c) - &
          dt*this%model%mesh%volume(c)*this%grad_p(:,c)
    end do
    if (this%model%inviscid) then
      call this%model%momentum%solve_inviscid(this%model%matl_props%density_c, this%rhs, &
          this%state%vel_cc(:,1:size(this%rhs,2)))
    else
      call this%momentum_solver%setup()
      call this%momentum_solver%solve(this%rhs, this%state%vel_cc(:,1:size(this%rhs,2)), stat)
      if (stat /= 0) return
    end if
    call this%model%mesh%cell_imap%gather_offp(this%state%vel_cc)
    call this%projection_update%correct(dt, this%model%matl_props%inv_density_c, &
        this%model%matl_props%inv_density_f, this%model%matl_props%density_delta_c, this%model%bc, this%state, stat)
    if (stat == 0) this%nstep = this%nstep + 1_int64
  end subroutine


  !! Declare the temporal scalar fields published by this solver.
  !! These fields are updated at each requested solution output and written
  !! by the simulation's output writer.
  subroutine init_temporal_output(this, data)
    class(flow_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine


  !! Set the current values of the temporal scalar fields published by this
  !! solver.
  subroutine set_temporal_output(this, data)
    class(flow_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data

    call data%set('NStep', this%nstep)
  end subroutine

end module flow_2d_solver_type
