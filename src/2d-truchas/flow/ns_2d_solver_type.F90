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
  use simulation_environment_type
  use parameter_list_type
  use parallel_communication, only: global_minval
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_material_transport_type
  use flow_2d_momentum_solver_type
  use flow_2d_projection_solver_type
  use flow_2d_projection_update_type
  use flow_2d_ic_solver_type
  implicit none
  private

  type, public :: ns_2d_solver
    private
    type(flow_2d_model), pointer :: model => null()  ! unowned reference
    type(flow_2d_state) :: state
    type(flow_2d_momentum_solver) :: momentum_solver
    type(flow_2d_projection_solver), pointer :: projection_solver => null()
    type(flow_2d_projection_update) :: projection_update
    type(flow_2d_ic_solver), pointer :: ic_solver => null()
    type(flow_2d_material_transport) :: material_transport
    real(r8), allocatable :: rhs(:,:), grad_p(:,:)
  contains
    procedure :: init
    procedure :: set_buoyancy_temperature
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: step
    procedure :: advance_momentum
    procedure :: courant_time_step
    final :: delete
  end type

contains

  subroutine init(this, env, model, momentum_params, projection_params)
    class(ns_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(in) :: model
    type(parameter_list), target, intent(in), optional :: momentum_params
    type(parameter_list), target, intent(in) :: projection_params

    this%model => model
    call this%state%init(model%mesh)
    allocate(this%rhs(2, model%mesh%ncell_onP), this%grad_p(2, model%mesh%ncell), &
        this%projection_solver, this%ic_solver)
    call this%material_transport%init(env, model%mesh, size(model%density))
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


  subroutine set_buoyancy_temperature(this, temperature)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    call this%model%set_buoyancy_temperature(temperature)
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


  !! Return no-copy views of the cell-centered pressure and velocity.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(ns_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    pressure => this%state%p_cc
    velocity => this%state%vel_cc
  end subroutine


  !! Return a no-copy view of the face-normal velocity.
  subroutine get_face_velocity(this, velocity)
    class(ns_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    velocity => this%state%vel_fn
  end subroutine


  !! Advance STATE from T_N to T_NP1. The time step is derived from the two
  !! endpoint times so callers retain exact target times. This is the
  !! isothermal wrapper: it first obtains material transport from the old face
  !! velocity and then advances momentum and pressure.
  subroutine step(this, t_n, t_np1, stat, errmsg)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg

    call this%material_transport%advance(t_n, t_np1, this%state%vel_fn)
    call this%advance_momentum(t_n, t_np1, this%material_transport%flux_volumes, stat, errmsg)
  end subroutine


  !! Advance momentum and pressure from T_N to T_NP1 using material-resolved
  !! flux volumes already constructed for the pending step. This separate
  !! operation lets a coupled solver advect material and thermal enthalpy
  !! before it updates the flow state.
  subroutine advance_momentum(this, t_n, t_np1, flux_volumes, stat, errmsg)
    class(ns_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t_n, t_np1
    real(r8), intent(in) :: flux_volumes(:,:)
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg

    integer :: c
    real(r8) :: dt
    character(:), allocatable :: bc_errmsg

    dt = t_np1 - t_n
    ASSERT(dt > 0.0_r8)
    ASSERT(size(flux_volumes,1) == size(this%model%density))
    ASSERT(size(flux_volumes,2) == size(this%model%mesh%cface))
    call this%model%compute_bc(t_n, dt, stat, bc_errmsg)
    if (stat /= 0) then
      if (present(errmsg)) errmsg = bc_errmsg
      return
    end if
    call this%model%pressure_gradient(this%state%p_cc, this%grad_p)
    call this%model%assemble_momentum(dt, this%rhs)
    call this%model%momentum%add_advective_rhs(this%model%density, this%state%vel_cc, &
        flux_volumes, &
        this%model%bc, this%rhs)
    do c = 1, size(this%rhs,2)
      this%rhs(:,c) = this%rhs(:,c) + this%model%density_c(c)*this%model%mesh%volume(c)*this%state%vel_cc(:,c) - &
          dt*this%model%mesh%volume(c)*this%grad_p(:,c)
    end do
    if (this%model%inviscid) then
      call this%model%momentum%solve_inviscid(this%model%density_c, this%rhs, &
          this%state%vel_cc(:,1:size(this%rhs,2)))
    else
      call this%momentum_solver%setup()
      call this%momentum_solver%solve(this%rhs, this%state%vel_cc(:,1:size(this%rhs,2)), stat)
      if (stat /= 0) return
    end if
    call this%model%mesh%cell_imap%gather_offp(this%state%vel_cc)
    call this%projection_update%correct(dt, this%model%inv_density_c, this%model%inv_density_f, &
        this%model%density_delta_c, this%model%bc, this%state, stat)
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
