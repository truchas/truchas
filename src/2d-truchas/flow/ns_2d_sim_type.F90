!!
!! NS_2D_SIM_TYPE
!!
!! This module defines NS_2D_SIM, which owns the mesh, single-fluid
!! Navier--Stokes model, state, solver, output, and fixed maximum time-step
!! control for a two-dimensional isothermal incompressible-flow simulation.
!! The actual step is further restricted by a convective Courant limit.
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
  use flow_2d_state_type
  use ns_2d_solver_type
  use flow_2d_vtkhdf_output
  implicit none
  private

  type, public :: ns_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(flow_2d_model), pointer :: model => null()
    type(flow_2d_state), pointer :: state => null()
    type(ns_2d_solver), pointer :: solver => null()
    type(flow_2d_vtkhdf_writer) :: output
    real(r8) :: initial_time, time_step, final_time, courant_number
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
    if (associated(this%state)) deallocate(this%state)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine


  subroutine init(this, params, stat, errmsg)
    class(ns_2d_sim), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: mesh_params, model_params, bc_params, solver_params
    type(parameter_list), pointer :: momentum_params, projection_params, control_params
    class(vector_func), allocatable :: initial_velocity_func
    real(r8), allocatable :: initial_velocity(:,:)
    real(r8) :: density, viscosity
    real(r8), allocatable :: body_acceleration(:)

    stat = 0
    if (.not.params%is_sublist('mesh')) then
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
      return
    end if
    mesh_params => params%sublist('mesh')
    this%mesh => new_unstr_2d_mesh(mesh_params, stat, errmsg)
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
    call this%model%init(this%mesh, bc_params, density, viscosity, stat, errmsg, body_acceleration)
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
    call control_params%get('initial-time', this%initial_time, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call control_params%get('time-step', this%time_step, stat, errmsg)
    if (stat /= 0) return
    call control_params%get('final-time', this%final_time, stat, errmsg)
    if (stat /= 0) return
    call control_params%get('courant-number', this%courant_number, default=0.5_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%time_step <= 0.0_r8 .or. this%final_time < this%initial_time .or. &
        this%courant_number <= 0.0_r8 .or. this%courant_number > 1.0_r8) then
      stat = 1
      errmsg = 'processing ' // control_params%path() // &
          ': require time-step > 0, final-time >= initial-time, and courant-number in (0,1]'
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
    call this%solver%set_initial_state(this%initial_time, this%time_step, initial_velocity, stat)
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
    class(ns_2d_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: time, dt

    stat = 0
    time = this%initial_time
    call this%output%write_solution(time, this%state%p_cc, this%state%vel_cc)
    do while (time < this%final_time)
      dt = min(this%time_step, this%solver%courant_time_step(this%courant_number), this%final_time - time)
      call this%solver%step(time, dt, stat, errmsg)
      if (stat /= 0) then
        if (.not.allocated(errmsg)) errmsg = 'Navier--Stokes solver step failed'
        return
      end if
      time = time + dt
      call this%output%write_solution(time, this%state%p_cc, this%state%vel_cc)
    end do
    call this%output%close()
  end subroutine

end module ns_2d_sim_type
