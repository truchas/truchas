!!
!! FLOW_2D_SIM_TYPE
!!
!! This module defines FLOW_2D_SIM, which owns the mesh, single-fluid flow
!! model, state, solver, time-integration control, and output writer for a
!! two-dimensional isothermal incompressible-flow simulation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module flow_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use vector_func_class
  use vector_func_factories, only: alloc_vector_func
  use vector_func_projection
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_solver_type
  use flow_2d_vtkhdf_output
  use time_step_sync_type
  use simulation_environment_type
  implicit none
  private

  type, public :: flow_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(flow_2d_model), pointer :: model => null()
    type(flow_2d_state), pointer :: state => null()
    type(flow_2d_solver), pointer :: solver => null()
    type(flow_2d_vtkhdf_writer) :: output
    real(r8) :: t_init, tlast, hlast
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    real(r8), allocatable :: tout(:)
    type(time_step_sync) :: ts_sync
  contains
    final :: delete
    procedure :: init
    procedure :: run
  end type

contains

  subroutine delete(this)
    type(flow_2d_sim), intent(inout) :: this

    call this%output%close()
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%state)) deallocate(this%state)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine


  subroutine init(this, params, stat, errmsg, env)

    use signal_handler, only: init_signal_handler, SIGURG
    class(flow_2d_sim), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simulation_environment), intent(in) :: env

    type(parameter_list), pointer :: mesh_params, model_params, bc_params, solver_params
    type(parameter_list), pointer :: momentum_params, projection_params, control_params
    class(vector_func), allocatable :: initial_velocity_func
    real(r8), allocatable :: initial_velocity(:,:)
    real(r8) :: density, viscosity
    real(r8), allocatable :: body_acceleration(:)

    stat = 0
    call init_signal_handler(SIGURG)
    if (.not.params%is_sublist('mesh')) then
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
      return
    end if
    mesh_params => params%sublist('mesh')
    this%mesh => new_unstr_2d_mesh(mesh_params, stat, errmsg, env%simlog)
    if (stat /= 0) then
      errmsg = 'processing ' // mesh_params%path() // ': ' // errmsg
      return
    end if

    if (.not.params%is_sublist('flow-model')) then
      stat = 1
      errmsg = 'missing "flow-model" sublist parameter'
      return
    end if
    model_params => params%sublist('flow-model')
    call model_params%get('density', density, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    call model_params%get('viscosity', viscosity, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    call model_params%get('body-acceleration', body_acceleration, stat=stat, errmsg=errmsg, &
        default=[0.0_r8, 0.0_r8])
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    if (density <= 0.0_r8 .or. viscosity <= 0.0_r8) then
      stat = 1
      errmsg = 'processing ' // model_params%path() // ': require density > 0 and viscosity > 0'
      return
    end if
    if (.not.model_params%is_sublist('bc')) then
      stat = 1
      errmsg = 'missing "bc" sublist parameter in ' // model_params%path()
      return
    end if
    bc_params => model_params%sublist('bc')
    allocate(this%model)
    call this%model%init(this%mesh, bc_params, density, viscosity, stat, errmsg, body_acceleration, env)
    if (stat /= 0) then
      errmsg = 'processing ' // bc_params%path() // ': ' // errmsg
      return
    end if
    allocate(this%state)
    call this%state%init(this%mesh)

    if (.not.params%is_sublist('flow-solver')) then
      stat = 1
      errmsg = 'missing "flow-solver" sublist parameter'
      return
    end if
    solver_params => params%sublist('flow-solver')
    if (.not.solver_params%is_sublist('momentum-solver') .or. &
        .not.solver_params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'flow-solver requires "momentum-solver" and "projection-solver" sublists'
      return
    end if
    momentum_params => solver_params%sublist('momentum-solver')
    projection_params => solver_params%sublist('projection-solver')
    allocate(this%solver)
    call this%solver%init(this%model, this%state, momentum_params, projection_params)

    if (.not.params%is_sublist('sim-control')) then
      stat = 1
      errmsg = 'missing "sim-control" sublist parameter'
      return
    end if
    control_params => params%sublist('sim-control')
    call control_params%get('initial-time', this%t_init, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call control_params%get('initial-time-step', this%dt_init, stat, errmsg)
    if (stat /= 0) return
    call control_params%get('min-time-step', this%dt_min, stat, errmsg)
    if (stat /= 0) return
    call control_params%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call control_params%get('time-step-growth', this%dt_grow, default=1.05_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call control_params%get('output-times', this%tout, stat, errmsg)
    if (stat /= 0) return
    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. size(this%tout) == 0 .or. &
        any(this%tout <= this%t_init) .or. any(this%tout(2:) <= this%tout(:size(this%tout)-1))) then
      stat = 1
      errmsg = 'processing ' // control_params%path() // ': require 0 < min-time-step <= initial-time-step <= ' // &
          'max-time-step, time-step-growth >= 1, and strictly increasing output-times after initial-time'
      return
    end if
    this%ts_sync = time_step_sync(4)

    allocate(initial_velocity(2,this%mesh%ncell_onP), source=0.0_r8)
    if (params%is_parameter('initial-velocity')) then
      call alloc_vector_func(params, 'initial-velocity', initial_velocity_func, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing initial-velocity: ' // errmsg
        return
      end if
      if (initial_velocity_func%dim /= 2) then
        stat = 1
        errmsg = 'processing initial-velocity: require a two-component vector function'
        return
      end if
      call project_vector_func_to_cell_centers(this%mesh, initial_velocity_func, initial_velocity)
    end if
    call this%solver%set_initial_state(this%t_init, this%dt_init, initial_velocity, stat)
    if (stat /= 0) then
      errmsg = 'initializing flow state failed'
      return
    end if

    call this%output%open(this%mesh, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'opening VTKHDF output: ' // errmsg
      return
    end if

  end subroutine


  subroutine run(this, stat, errmsg)
    class(flow_2d_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: time, hnext, t_write

    stat = 0
    time = this%t_init
    call this%output%write_solution(time, this%state%p_cc, this%state%vel_cc)
    t_write = time
    hnext = this%dt_init
    this%tlast = time
    this%hlast = hnext
    do n = 1, size(this%tout)
      call integrate(this, this%tout(n), hnext, time, stat, errmsg)
      if (stat < 0 .and. time == t_write) exit
      call this%output%write_solution(time, this%state%p_cc, this%state%vel_cc)
      t_write = time
      if (stat /= 0) exit
    end do
    if (stat > 0) then
      stat = 0
      deallocate(errmsg)
    end if
    call this%output%close()
  end subroutine


  subroutine integrate(this, tout, hnext, time, stat, errmsg)

    use signal_handler, only: read_signal, SIGURG

    class(flow_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: tout
    real(r8), intent(inout) :: hnext
    real(r8), intent(out) :: time
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical :: sig_rcvd
    real(r8) :: dt

    do
      time = this%ts_sync%next_time(tout, this%tlast, this%hlast, hnext)
      dt = time - this%tlast
      if (dt < this%dt_min) then
        stat = -1
        errmsg = 'next time step is too small'
        time = this%tlast
        return
      end if
      call this%solver%step(this%tlast, dt, stat, errmsg)
      if (stat /= 0) then
        if (.not.allocated(errmsg)) errmsg = 'flow solver step failed'
        time = this%tlast
        return
      end if
      this%hlast = dt
      this%tlast = time
      hnext = min(this%dt_grow*dt, this%dt_max)

      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if
      if (time == tout) return
    end do
  end subroutine

end module flow_2d_sim_type
