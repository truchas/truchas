!!
!! FLOW_2D_SOLVER_TYPE
!!
!! This module defines FLOW_2D_SOLVER, the common two-dimensional flow
!! mechanics solver.  It advances the implicit momentum predictor and
!! incremental pressure correction, while allowing a caller to interpose
!! material or other coupled-physics transport before the momentum update.
!! It owns the committed and pending flow states and provides views of them
!! to its caller.  Material transport and time-step orchestration belong to
!! higher-level simulation solvers.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use parallel_communication, only: global_minval
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_momentum_solver_type
  use flow_2d_projection_solver_type
  use flow_2d_projection_update_type
  use flow_2d_ic_solver_type
  use flow_domain_types
  implicit none
  private

  type, public :: flow_2d_solver
    private
    type(flow_2d_model), pointer :: model => null()  ! unowned reference
    type(flow_2d_state) :: state
    type(flow_2d_state) :: pending_state
    type(flow_2d_momentum_solver) :: momentum_solver
    type(flow_2d_projection_solver), pointer :: projection_solver => null()
    type(flow_2d_projection_update) :: projection_update
    type(flow_2d_ic_solver), pointer :: ic_solver => null()
    real(r8), allocatable :: rhs(:,:), grad_p(:,:)
    real(r8) :: courant_number = 0.5_r8
    integer(int64) :: nstep = 0_int64
    logical :: step_is_pending = .false.
  contains
    procedure :: init
    procedure :: set_volume_fractions
    procedure :: set_initial_material_state
    procedure :: set_buoyancy_temperature
    procedure :: set_initial_state
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: step
    procedure :: advance_momentum
    procedure :: commit_step
    procedure :: reject_step
    procedure :: courant_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    final :: delete
  end type

contains

  subroutine init(this, env, model, momentum_params, projection_params, courant_number, stat, errmsg)
    class(flow_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(in) :: model
    type(parameter_list), target, intent(in), optional :: momentum_params
    type(parameter_list), target, intent(in) :: projection_params
    real(r8), intent(in), optional :: courant_number
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    this%model => model
    if (present(courant_number)) then
      ASSERT(courant_number > 0.0_r8 .and. courant_number <= 1.0_r8)
      this%courant_number = courant_number
    end if
    call this%state%init(model%mesh)
    call this%pending_state%init(model%mesh)
    allocate(this%rhs(2, model%mesh%ncell_onP), this%grad_p(2, model%mesh%ncell), this%projection_solver, &
        this%ic_solver)
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
    if (present(momentum_params)) then
      call this%ic_solver%init(model, momentum_params, projection_params, stat, errmsg)
    else
      call this%ic_solver%init(model, projection_params=projection_params, stat=stat, errmsg=errmsg)
    end if
    if (present(stat)) then
      if (stat /= 0) return
    end if
  end subroutine


  subroutine set_volume_fractions(this, vfrac)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    call this%model%set_volume_fractions(vfrac)
  end subroutine


  subroutine set_initial_material_state(this, vfrac, temperature)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%model%set_initial_material_state(vfrac, temperature)
  end subroutine


  subroutine set_buoyancy_temperature(this, temperature)
    class(flow_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    call this%model%set_buoyancy_temperature(temperature)
  end subroutine


  subroutine set_initial_state(this, env, time, dt, velocity, stat)
    class(flow_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: time, dt, velocity(:,:)
    integer, intent(out) :: stat

    call this%ic_solver%solve(env, time, dt, velocity, this%state, stat)
    if (stat /= 0) return
    this%pending_state%vel_cc = this%state%vel_cc
    this%pending_state%vel_fn = this%state%vel_fn
    this%pending_state%p_cc = this%state%p_cc
    this%step_is_pending = .false.
    this%nstep = 0_int64
  end subroutine


  !! Return no-copy views of the cell-centered pressure and velocity.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(flow_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    if (this%step_is_pending) then
      pressure => this%pending_state%p_cc
      velocity => this%pending_state%vel_cc
    else
      pressure => this%state%p_cc
      velocity => this%state%vel_cc
    end if
  end subroutine


  !! Return a no-copy view of the face-normal velocity.
  subroutine get_face_velocity(this, velocity)
    class(flow_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    if (this%step_is_pending) then
      velocity => this%pending_state%vel_fn
    else
      velocity => this%state%vel_fn
    end if
  end subroutine


  subroutine delete(this)
    type(flow_2d_solver), intent(inout) :: this

    if (associated(this%projection_solver)) deallocate(this%projection_solver)
    if (associated(this%ic_solver)) deallocate(this%ic_solver)
  end subroutine


  !! Advance the flow mechanics from T_N to T_NP1.  The time step is derived
  !! from the two endpoint times so callers retain exact target times.  If
  !! FLUX_VOLUMES is present, its material-resolved values provide the
  !! explicit momentum-advection contribution.  The result remains pending
  !! until COMMIT_STEP is called.
  subroutine step(this, env, t_n, t_np1, stat, errmsg)
    class(flow_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg
    character(256) :: line

    write(line,'(a,i0,a,es0.5,a,es0.5,a)') 'step=', this%nstep + 1_int64, &
        ' attempt=1 t0=', t_n, ' dt=', t_np1 - t_n, ' cause=explicit'
    call env%simlog%begin_section(trim(line))
    call this%advance_momentum(env, t_n, t_np1, stat, errmsg)
    if (stat == 0) then
      call this%commit_step()
      this%nstep = this%nstep + 1_int64
      call env%simlog%end_section('step-end status=accepted')
    else
      call env%simlog%end_section('step-end status=failed')
    end if
  end subroutine


  subroutine advance_momentum(this, env, t_n, t_np1, stat, errmsg, flux_volumes)
    class(flow_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg
    real(r8), intent(in), optional :: flux_volumes(:,:)

    integer :: c, num_itr, num_dscg_itr, num_pcg_itr
    real(r8) :: dt
    real(r8) :: rel_res_norm
    character(:), allocatable :: bc_errmsg
    logical :: projection_solved

    ASSERT(.not.this%step_is_pending)
    dt = t_np1 - t_n
    ASSERT(dt > 0.0_r8)
    call env%timer%start('flow/momentum')
    this%pending_state%vel_cc = this%state%vel_cc
    this%pending_state%vel_fn = this%state%vel_fn
    this%pending_state%p_cc = this%state%p_cc
    call this%model%compute_bc(t_n, dt, stat, bc_errmsg)
    if (stat /= 0) then
      call env%timer%stop('flow/momentum')
      if (present(errmsg)) errmsg = bc_errmsg
      return
    end if
    call this%model%pressure_gradient(this%state%p_cc, this%grad_p)
    call this%model%assemble_momentum(dt, this%rhs)
    if (present(flux_volumes)) then
      ASSERT(size(flux_volumes,1) == size(this%model%matl_props%density))
      ASSERT(size(flux_volumes,2) == size(this%model%mesh%cface))
      call this%model%momentum%add_advective_rhs(this%model%matl_props%density, this%state%vel_cc, &
          flux_volumes, this%model%matl_props%cell_t, this%model%matl_props%face_t, this%model%bc, this%rhs)
    end if
    do c = 1, size(this%rhs,2)
      if (this%model%matl_props%cell_t(c) > regular_t) then
        this%rhs(:,c) = 0.0_r8
        cycle
      end if
      this%rhs(:,c) = this%rhs(:,c) + this%model%matl_props%density_c_old(c)*this%model%mesh%volume(c)* &
          this%state%vel_cc(:,c) - dt*this%model%mesh%volume(c)*this%grad_p(:,c)
    end do
    if (this%model%inviscid) then
      call this%model%momentum%solve_inviscid(this%model%matl_props%density_c, &
          this%model%matl_props%cell_t, this%rhs, &
          this%pending_state%vel_cc(:,1:size(this%rhs,2)))
      call env%timer%stop('flow/momentum')
      call env%simlog%info('flow.momentum method=inviscid-direct status=ok')
    else
      call this%momentum_solver%setup()
      call this%momentum_solver%solve(this%rhs, this%pending_state%vel_cc(:,1:size(this%rhs,2)), stat)
      call env%timer%stop('flow/momentum')
      call this%momentum_solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
      call write_solver_metrics(env, 'flow.momentum', num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm, stat)
      if (stat /= 0) return
    end if
    call this%model%mesh%cell_imap%gather_offp(this%pending_state%vel_cc)
    call env%timer%start('flow/projection')
    call this%projection_update%correct(dt, this%model%matl_props%inv_density_c, &
        this%model%matl_props%inv_density_f, this%model%matl_props%density_delta_c, &
        this%model%matl_props%cell_t, this%model%matl_props%face_t, this%model%bc, &
        this%pending_state, stat, solved=projection_solved)
    call env%timer%stop('flow/projection')
    if (projection_solved) then
      call this%projection_solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
      call write_solver_metrics(env, 'flow.projection', num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm, stat)
    else
      call env%simlog%info('flow.projection method=none reason=zero-rhs status=skipped')
    end if
    if (stat /= 0) return
    this%step_is_pending = .true.
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
    write(line,'(a,i0,a,i0,a,i0,a,es0.5,a,a)') trim(name) // &
        ' iterations=', num_itr, ' dscg=', num_dscg_itr, ' amg=', num_pcg_itr, &
        ' rel_residual=', rel_res_norm, ' status=', trim(status)
    call env%simlog%info(trim(line))
  end subroutine


  subroutine commit_step(this)
    class(flow_2d_solver), intent(inout) :: this

    if (this%step_is_pending) then
      this%state%vel_cc = this%pending_state%vel_cc
      this%state%vel_fn = this%pending_state%vel_fn
      this%state%p_cc = this%pending_state%p_cc
      call this%model%accept_material_state()
      this%step_is_pending = .false.
    end if
  end subroutine


  subroutine reject_step(this)
    class(flow_2d_solver), intent(inout) :: this

    this%step_is_pending = .false.
  end subroutine


  !! Return the maximum step size satisfying the configured convective Courant
  !! number for the current face-normal velocity.
  function courant_time_step(this) result(dt)
    class(flow_2d_solver), intent(in) :: this
    real(r8) :: dt

    integer :: f
    real(r8) :: dt_local
    real(r8), pointer :: velocity(:)

    ASSERT(this%courant_number > 0.0_r8 .and. this%courant_number <= 1.0_r8)
    dt_local = huge(1.0_r8)
    call this%get_face_velocity(velocity)
    do f = 1, this%model%mesh%nface_onP
      if (velocity(f) == 0.0_r8) cycle
      dt_local = min(dt_local, this%model%operators%normal_distance(f)/abs(velocity(f)))
    end do
    dt = this%courant_number*global_minval(dt_local)
  end function


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
