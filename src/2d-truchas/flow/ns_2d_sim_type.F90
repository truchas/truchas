!!
!! NS_2D_SIM_TYPE
!!
!! This module defines NS_2D_SIM, which owns the mesh, material composition,
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
  use material_database_type
  use material_model_type
  use material_composition_type
  use flow_2d_material_layout_type
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
    type(material_database) :: matl_db
    type(material_model) :: matl_model
    type(material_composition), pointer :: composition => null()
    type(flow_2d_model), pointer :: model => null()
    type(ns_2d_solver), pointer :: solver => null()
    type(flow_2d_vtkhdf_writer) :: output
    type(parameter_list) :: temporal_output
    real(r8) :: t_init
    real(r8), allocatable :: tout(:)
  contains
    final :: delete
    procedure :: init
    procedure :: run
    procedure :: write_solution
  end type

contains

  subroutine delete(this)
    type(ns_2d_sim), intent(inout) :: this

    call this%output%close()
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%composition)) deallocate(this%composition)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine


  subroutine init(this, env, params, stat, errmsg)

    use material_factory, only: load_material_database
    use signal_handler, only: init_signal_handler, SIGURG
    class(ns_2d_sim), intent(out) :: this
    type(simulation_environment), intent(inout) :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: mesh_params, model_params, bc_params, solver_params, tracking_params, materials_params
    type(parameter_list), pointer :: momentum_params, projection_params, control_params
    type(parameter_list_iterator) :: piter
    type(flow_2d_material_layout) :: material_layout
    class(vector_func), allocatable :: initial_velocity_func
    real(r8), allocatable :: initial_velocity(:,:)
    real(r8), allocatable :: body_acceleration(:)
    real(r8), allocatable :: vfrac(:,:), temperature(:)
    character(:), allocatable :: matl_name(:), region_name(:), region_matl_name(:)
    integer, allocatable :: fluid_material_ids(:)
    integer :: i, rlev
    logical :: inviscid

    stat = 0
    call init_signal_handler(SIGURG)
    if (.not.params%is_sublist('mesh')) then
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
      return
    end if
    mesh_params => params%sublist('mesh')
    call env%simlog%begin_section('Constructing mesh.')
    this%mesh => new_unstr_2d_mesh(env, mesh_params, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Mesh construction failed.')
      errmsg = 'processing ' // mesh_params%path() // ': ' // errmsg
      return
    end if
    call env%simlog%end_section('Mesh construction complete.')

    if (.not.params%is_sublist('materials')) then
      stat = 1
      errmsg = 'missing "materials" sublist parameter'
      return
    end if
    materials_params => params%sublist('materials')
    call env%simlog%info('Loading material database.')
    call load_material_database(this%matl_db, materials_params, stat, errmsg)
    if (stat /= 0) return
    allocate(this%composition)
    if (params%is_sublist('material-regions')) then
      call env%simlog%begin_section('Reading material regions.')
      call get_material_region_names(params%sublist('material-regions'), matl_name, stat, errmsg, &
          region_name, region_matl_name)
      if (stat /= 0) then
        call env%simlog%end_section('Material-region processing failed.')
        return
      end if
      do i = 1, size(region_name)
        call env%simlog%info('Region "' // trim(region_name(i)) // '": material="' // &
            trim(region_matl_name(i)) // '".')
      end do
      call env%simlog%end_section('Material regions read.')
      call params%get('material-region-refinement-level', rlev, default=6, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (rlev < 0) then
        stat = 1
        errmsg = '"material-region-refinement-level" must be >= 0'
        return
      end if
    else
      piter = parameter_list_iterator(materials_params, sublists_only=.true.)
      if (piter%count() /= 1) then
        stat = 1
        errmsg = 'multiple materials require a "material-regions" sublist'
        return
      end if
      matl_name = [piter%name()]
    end if
    call env%simlog%begin_section('Constructing material model.')
    call this%matl_model%init(matl_name, this%matl_db, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Material model construction failed.')
      return
    end if
    do i = 1, size(matl_name)
      call env%simlog%info('Using material "' // trim(matl_name(i)) // '".')
    end do
    call env%simlog%end_section('Material model complete.')
    call env%simlog%begin_section('Constructing material distribution.')
    if (params%is_sublist('material-regions')) then
      call this%composition%init(env, this%mesh, this%matl_model, params%sublist('material-regions'), rlev, stat, errmsg)
    else
      call env%simlog%info('Assigning material "' // trim(matl_name(1)) // '" uniformly.')
      call this%composition%init_uniform(this%mesh, this%matl_model, 1, stat, errmsg)
    end if
    if (stat /= 0) then
      call env%simlog%end_section('Material distribution construction failed.')
      return
    end if
    call env%simlog%end_section('Material distribution complete.')

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
    call model_params%get('body-acceleration', body_acceleration, stat=stat, errmsg=errmsg, &
        default=[0.0_r8, 0.0_r8])
    if (stat /= 0) then
      errmsg = 'processing ' // model_params%path() // ': ' // errmsg
      return
    end if
    if (.not.model_params%is_sublist('bc')) then
      stat = 1
      errmsg = 'missing "bc" sublist parameter in ' // model_params%path()
      return
    end if
    bc_params => model_params%sublist('bc')
    if (.not.params%is_sublist('flow-solver')) then
      stat = 1
      errmsg = 'missing "flow-solver" sublist parameter'
      return
    end if
    solver_params => params%sublist('flow-solver')
    if (.not.solver_params%is_sublist('volume-tracking')) then
      stat = 1
      errmsg = 'flow-solver requires a volume-tracking sublist for material flow'
      return
    end if
    tracking_params => solver_params%sublist('volume-tracking')
    call material_layout%init(this%matl_model, tracking_params, stat, errmsg)
    if (stat /= 0) return
    allocate(fluid_material_ids(material_layout%num_real_fluid()))
    call material_layout%get_real_fluid_material_ids(fluid_material_ids)
    allocate(this%model)
    call this%model%init_material(env, this%mesh, bc_params, this%matl_model, fluid_material_ids, stat, errmsg, &
        body_acceleration=body_acceleration, inviscid=inviscid)
    if (stat /= 0) then
      errmsg = 'processing ' // bc_params%path() // ': ' // errmsg
      return
    end if
    if (.not.solver_params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'flow-solver requires a "projection-solver" sublist'
      return
    end if
    projection_params => solver_params%sublist('projection-solver')
    allocate(this%solver)
    if (inviscid) then
      call this%solver%init(env, this%model, projection_params=projection_params, material_layout=material_layout, &
          tracking_params=tracking_params)
    else
      if (.not.solver_params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a "momentum-solver" sublist'
        return
      end if
      momentum_params => solver_params%sublist('momentum-solver')
      call this%solver%init(env, this%model, momentum_params, projection_params, material_layout, tracking_params)
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
    allocate(vfrac(material_layout%num_material(), this%mesh%ncell), temperature(this%mesh%ncell_onP), source=0.0_r8)
    call material_layout%get_reduced_volume_fractions(this%composition, vfrac)
    call this%mesh%cell_imap%gather_offp(vfrac)
    call this%solver%set_initial_material_state(vfrac, temperature)
    call this%solver%set_initial_state(this%t_init, this%solver%initial_time_step(), initial_velocity, stat)
    if (stat /= 0) then
      errmsg = 'initializing flow state failed'
      return
    end if

    call this%solver%init_temporal_output(this%temporal_output)
    call this%output%open(env, this%mesh, this%temporal_output, stat, errmsg, this%matl_model)
    if (stat /= 0) then
      errmsg = 'opening VTKHDF output: ' // errmsg
      return
    end if
  end subroutine


  subroutine run(this, env, stat, errmsg)
    class(ns_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: time, t_write

    stat = 0
    if (associated(env%timer)) call env%timer%start('integration')
    time = this%solver%last_time()
    if (associated(env%timer)) call env%timer%start('output')
    call this%write_solution(time)
    if (associated(env%timer)) call env%timer%stop('output')
    t_write = time
    do n = 1, size(this%tout)
      call this%solver%integrate(env, this%tout(n), stat, errmsg)
      time = this%solver%last_time()
      if (stat < 0 .and. time == t_write) exit
      if (associated(env%timer)) call env%timer%start('output')
      call this%write_solution(time)
      if (associated(env%timer)) call env%timer%stop('output')
      t_write = time
      if (stat /= 0) exit
    end do
    if (stat > 0) then
      stat = 0
      deallocate(errmsg)
    end if
    call this%output%close()
    if (associated(env%timer)) call env%timer%stop('integration')
  end subroutine


  subroutine write_solution(this, time)
    class(ns_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: time

    real(r8), pointer :: pressure(:), velocity(:,:), vfrac(:,:)

    call this%solver%get_cell_flow_soln(pressure, velocity)
    call this%solver%set_temporal_output(this%temporal_output)
    call this%solver%get_volume_fractions(vfrac)
    call this%output%write_solution(time, pressure, velocity, this%temporal_output, vfrac)
  end subroutine


end module ns_2d_sim_type
