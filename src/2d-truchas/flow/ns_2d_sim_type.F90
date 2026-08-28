!!
!! NS_2D_SIM_TYPE
!!
!! This module defines NS_2D_SIM, which owns the mesh, single-fluid
!! Navier--Stokes model, solver, output schedule, and output writer for a
!! two-dimensional isothermal incompressible-flow simulation.  The solver
!! owns the flow state and time-step policy; the simulation supplies target
!! output times.  Each step is restricted by a convective Courant limit.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ns_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use vector_func_class
  use vector_func_factories, only: alloc_vector_func
  use vector_func_projection
  use flow_2d_model_type
  use ns_2d_solver_type
  use flow_2d_vtkhdf_writer_type
  use simulation_environment_type
  implicit none
  private

  type, public :: ns_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(flow_2d_model), pointer :: model => null()
    type(ns_2d_solver), pointer :: solver => null()
    type(flow_2d_vtkhdf_writer) :: output
    real(r8) :: t_init
    real(r8), allocatable :: tout(:)
  contains
    final :: delete
    procedure :: init
    procedure :: run
  end type

contains

  subroutine delete(this)
    type(ns_2d_sim), intent(inout) :: this

    call this%output%close()
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine


  subroutine init(this, env, params, stat, errmsg)

    use signal_handler, only: init_signal_handler, SIGURG
    class(ns_2d_sim), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: mesh_params, model_params, bc_params, solver_params
    type(parameter_list), pointer :: momentum_params, projection_params, control_params
    class(vector_func), allocatable :: initial_velocity_func
    real(r8), allocatable :: initial_velocity(:,:)
    real(r8) :: density, viscosity
    real(r8), allocatable :: body_acceleration(:)
    logical :: inviscid

    stat = 0
    call init_signal_handler(SIGURG)
    if (.not.params%is_sublist('mesh')) then
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
      return
    end if
    mesh_params => params%sublist('mesh')
    this%mesh => new_unstr_2d_mesh(env, mesh_params, stat, errmsg)
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
    call model_params%get('inviscid', inviscid, default=.false., stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    call model_params%get('density', density, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    if (.not.inviscid) then
      call model_params%get('viscosity', viscosity, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // model_params%path() // ': ' // errmsg
        return
      end if
    end if
    call model_params%get('body-acceleration', body_acceleration, stat=stat, errmsg=errmsg, &
        default=[0.0_r8, 0.0_r8])
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    if (inviscid) then
      if (density <= 0.0_r8) then
        stat = 1
        errmsg = 'processing ' // model_params%path() // ': require density > 0'
        return
      end if
    else if (density <= 0.0_r8 .or. viscosity <= 0.0_r8) then
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
    if (inviscid) then
      call this%model%init(env, this%mesh, bc_params, density=[density], stat=stat, errmsg=errmsg, &
          body_acceleration=body_acceleration, inviscid=.true.)
    else
      call this%model%init(env, this%mesh, bc_params, density=[density], viscosity=viscosity, stat=stat, &
          errmsg=errmsg, body_acceleration=body_acceleration)
    end if
    if (stat /= 0) then
      errmsg = 'processing ' // bc_params%path() // ': ' // errmsg
      return
    end if
    if (.not.params%is_sublist('flow-solver')) then
      stat = 1
      errmsg = 'missing "flow-solver" sublist parameter'
      return
    end if
    solver_params => params%sublist('flow-solver')
    if (.not.solver_params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'flow-solver requires a "projection-solver" sublist'
      return
    end if
    projection_params => solver_params%sublist('projection-solver')
    allocate(this%solver)
    if (inviscid) then
      call this%solver%init(env, this%model, projection_params=projection_params)
    else
      if (.not.solver_params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a "momentum-solver" sublist'
        return
      end if
      momentum_params => solver_params%sublist('momentum-solver')
      call this%solver%init(env, this%model, momentum_params, projection_params)
    end if

    if (.not.params%is_sublist('sim-control')) then
      stat = 1
      errmsg = 'missing "sim-control" sublist parameter'
      return
    end if
    control_params => params%sublist('sim-control')
    call control_params%get('initial-time', this%t_init, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call control_params%get('output-times', this%tout, stat, errmsg)
    if (stat /= 0) return
    if (size(this%tout) == 0) then
      stat = 1
      errmsg = 'processing ' // control_params%path() // ': require nonempty output-times'
      return
    end if
    if (any(this%tout <= this%t_init) .or. any(this%tout(2:) <= this%tout(:size(this%tout)-1))) then
      stat = 1
      errmsg = 'processing ' // control_params%path() // ': require strictly increasing output-times after initial-time'
      return
    end if
    call this%solver%init_time_stepper(control_params, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // control_params%path() // ': ' // errmsg
      return
    end if

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
    call this%solver%set_initial_state(this%t_init, this%solver%initial_time_step(), initial_velocity, stat)
    if (stat /= 0) then
      errmsg = 'initializing flow state failed'
      return
    end if

    call this%output%open(env, this%mesh, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'opening VTKHDF output: ' // errmsg
      return
    end if
  end subroutine


  subroutine run(this, stat, errmsg)
    class(ns_2d_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: time, t_write
    real(r8), pointer :: pressure(:), velocity(:,:)

    stat = 0
    time = this%solver%last_time()
    call this%solver%get_cell_flow_soln(pressure, velocity)
    call this%output%write_solution(time, pressure, velocity)
    t_write = time
    do n = 1, size(this%tout)
      call this%solver%integrate(this%tout(n), stat, errmsg)
      time = this%solver%last_time()
      if (stat < 0 .and. time == t_write) exit
      call this%output%write_solution(time, pressure, velocity)
      t_write = time
      if (stat /= 0) exit
    end do
    if (stat > 0) then
      stat = 0
      deallocate(errmsg)
    end if
    call this%output%close()
  end subroutine


end module ns_2d_sim_type
