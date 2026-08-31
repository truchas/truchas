!!
!! NS_HT_2D_SIM_TYPE
!!
!! This module defines NS_HT_2D_SIM, which owns the setup and output for a
!! single-liquid two-dimensional non-isothermal Navier--Stokes simulation.
!! Time-step policy is owned by its coupled solver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ns_ht_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use material_database_type
  use material_model_type
  use material_composition_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use scalar_func_projection
  use vector_func_class
  use vector_func_factories, only: alloc_vector_func
  use vector_func_projection
  use flow_2d_model_type
  use flow_2d_material_layout_type
  use ht_2d_model_type
  use ns_ht_2d_solver_type
  use ns_ht_2d_vtkhdf_writer_type
  use simulation_environment_type
  implicit none
  private

  type, public :: ns_ht_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(material_database) :: matl_db
    type(material_model) :: matl_model
    type(material_composition), pointer :: composition => null()
    type(flow_2d_model), pointer :: flow_model => null()
    type(ht_2d_model), pointer :: ht_model => null()
    type(ns_ht_2d_solver), pointer :: solver => null()
    type(ns_ht_2d_vtkhdf_writer) :: output
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
    type(ns_ht_2d_sim), intent(inout) :: this

    ! Finalization may occur at an unexpected time and is not collective.
    ! Leave collective HDF5 cleanup to the explicit simulation shutdown.
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%ht_model)) deallocate(this%ht_model)
    if (associated(this%flow_model)) deallocate(this%flow_model)
    if (associated(this%composition)) deallocate(this%composition)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine

  subroutine init(this, env, params, stat, errmsg)
    use material_factory, only: load_material_database
    use material_utilities, only: add_enthalpy_prop
    use signal_handler, only: init_signal_handler, SIGURG

    class(ns_ht_2d_sim), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist, flow_bc, solver_params, flow_solver_params, tracking_params
    type(parameter_list_iterator) :: piter
    type(flow_2d_material_layout) :: material_layout
    class(scalar_func), allocatable :: initial_temp
    class(vector_func), allocatable :: initial_velocity_func
    character(:), allocatable :: matl_name(:)
    real(r8), allocatable :: body_acceleration(:), temp(:), velocity(:,:)
    integer, allocatable :: fluid_material_ids(:)
    integer :: rlev
    logical :: inviscid

    stat = 0
    call env%simlog%info('Initializing the non-isothermal flow simulation')
    call init_signal_handler(SIGURG)
    if (.not.params%is_sublist('mesh')) then
      stat = 1; errmsg = 'missing "mesh" sublist parameter'; return
    end if
    plist => params%sublist('mesh')
    this%mesh => new_unstr_2d_mesh(env, plist, stat, errmsg)
    if (stat /= 0) return

    if (.not.params%is_sublist('materials')) then
      stat = 1; errmsg = 'missing "materials" sublist parameter'; return
    end if
    plist => params%sublist('materials')
    call load_material_database(this%matl_db, plist, stat, errmsg)
    if (stat /= 0) return
    allocate(this%composition)
    if (params%is_sublist('material-regions')) then
      plist => params%sublist('material-regions')
      call get_material_region_names(plist, matl_name, stat, errmsg)
      if (stat /= 0) return
      call params%get('material-region-refinement-level', rlev, stat, errmsg, default=6)
      if (stat /= 0) return
      if (rlev < 0) then
        stat = 1
        errmsg = '"material-region-refinement-level" must be >= 0'
        return
      end if
    else
      piter = parameter_list_iterator(plist, sublists_only=.true.)
      if (piter%count() /= 1) then
        stat = 1; errmsg = 'multiple materials require a "material-regions" sublist'; return
      end if
      matl_name = [piter%name()]
    end if
    call this%matl_model%init(matl_name, this%matl_db, stat, errmsg)
    if (stat /= 0) return
    if (params%is_sublist('material-regions')) then
      plist => params%sublist('material-regions')
      call this%composition%init(env, this%mesh, this%matl_model, plist, rlev, stat, errmsg)
    else
      call this%composition%init_uniform(this%mesh, this%matl_model, 1, stat, errmsg)
    end if
    if (stat /= 0) return
    call add_enthalpy_prop(this%matl_model, stat, errmsg)
    if (stat /= 0) return

    if (.not.params%is_sublist('flow-model')) then
      stat = 1; errmsg = 'missing "flow-model" sublist parameter'; return
    end if
    plist => params%sublist('flow-model')
    call plist%get('inviscid', inviscid, default=.false., stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call plist%get('body-acceleration', body_acceleration, stat=stat, errmsg=errmsg, default=[0.0_r8,0.0_r8])
    if (stat /= 0) return
    if (.not.plist%is_sublist('bc')) then
      stat = 1; errmsg = 'flow-model requires a "bc" sublist'; return
    end if
    flow_bc => plist%sublist('bc')
    if (.not.params%is_sublist('solver')) then
      stat = 1; errmsg = 'missing "solver" sublist'; return
    end if
    solver_params => params%sublist('solver')
    if (.not.solver_params%is_sublist('flow')) then
      stat = 1; errmsg = 'solver requires a "flow" sublist'; return
    end if
    flow_solver_params => solver_params%sublist('flow')
    tracking_params => flow_solver_params%sublist('volume-tracking')
    call material_layout%init(this%matl_model, tracking_params, stat, errmsg)
    if (stat /= 0) return
    allocate(fluid_material_ids(material_layout%num_real_fluid()))
    call material_layout%get_real_fluid_material_ids(fluid_material_ids)
    allocate(this%flow_model)
    call this%flow_model%init_material(env, this%mesh, flow_bc, this%matl_model, fluid_material_ids, stat, errmsg, &
        body_acceleration=body_acceleration, inviscid=inviscid, boussinesq=.true.)
    if (stat /= 0) return

    if (.not.params%is_sublist('ht-model')) then
      stat = 1; errmsg = 'missing "ht-model" sublist parameter'; return
    end if
    plist => params%sublist('ht-model')
    allocate(this%ht_model)
    call this%ht_model%init(env, this%mesh, this%matl_model, this%composition, plist, stat, errmsg, advection=.true.)
    if (stat /= 0) return

    allocate(this%solver)
    call this%solver%init(env, this%flow_model, this%ht_model, this%matl_model, this%composition, &
        solver_params, stat, errmsg)
    if (stat /= 0) return

    if (.not.params%is_sublist('sim-control')) then
      stat = 1; errmsg = 'missing "sim-control" sublist parameter'; return
    end if
    plist => params%sublist('sim-control')
    call plist%get('initial-time', this%t_init, default=0.0_r8, stat=stat, errmsg=errmsg); if (stat /= 0) return
    call plist%get('output-times', this%tout, stat, errmsg); if (stat /= 0) return
    if (size(this%tout) == 0 .or. any(this%tout <= this%t_init) .or. any(this%tout(2:) <= this%tout(:size(this%tout)-1))) then
      stat = 1; errmsg = 'output-times must be strictly increasing and after initial-time'; return
    end if

    allocate(velocity(2,this%mesh%ncell_onP), source=0.0_r8)
    allocate(temp(this%mesh%ncell_onP))
    if (params%is_parameter('initial-velocity')) then
      call alloc_vector_func(params, 'initial-velocity', initial_velocity_func, stat, errmsg)
      if (stat /= 0) return
      if (initial_velocity_func%dim /= 2) then
        stat = 1; errmsg = 'initial-velocity must be two-dimensional'; return
      end if
      call project_vector_func_to_cell_centers(this%mesh, initial_velocity_func, velocity)
    end if
    if (.not.params%is_parameter('initial-temperature')) then
      stat = 1; errmsg = 'missing "initial-temperature" parameter'; return
    end if
    call alloc_scalar_func(params, 'initial-temperature', initial_temp, stat, errmsg)
    if (stat /= 0) return
    call project_scalar_func_to_cell_centers(this%mesh, initial_temp, temp)
    call this%solver%set_initial_state(env, this%t_init, velocity, temp, stat, errmsg)
    if (stat /= 0) return
    call this%solver%init_temporal_output(this%temporal_output)
    call this%output%open(env, this%mesh, this%matl_model, this%temporal_output, stat, errmsg)
  end subroutine

  subroutine run(this, env, stat, errmsg)
    class(ns_ht_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: time, time_written
    character(80) :: string(2)

    call env%timer%start('integration')

    time = this%solver%last_time()
    call this%write_solution(env, time)
    time_written = time
    call env%simlog%info('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', time
    call env%simlog%info(string(1))

    do n = 1, size(this%tout)
      call this%solver%integrate(env, this%tout(n), stat, errmsg)
      time = this%solver%last_time()
      if (stat < 0 .and. time == time_written) exit
      call this%write_solution(env, time)
      time_written = time
      call this%solver%write_metrics(string)
      call env%simlog%info('')
      call env%simlog%info(string(1))
      call env%simlog%info(string(2))
      if (stat /= 0) exit
    end do

    if (stat > 0) then
      call env%simlog%info('')
      call env%simlog%info(errmsg // ': current solution written, and now terminating ...')
      stat = 0
      deallocate(errmsg)
    else if (stat < 0) then
      call env%simlog%info('')
      errmsg = 'unrecoverable integration failure: ' // errmsg
      call env%simlog%info(errmsg)
    end if

    call env%simlog%info('')
    write(string(1),'(a,es12.5)') 'Completed integration to T = ', time
    call env%simlog%info(string(1))
    call this%output%close()
    call env%timer%stop('integration')
  end subroutine


  subroutine write_solution(this, env, time)
    class(ns_ht_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: time

    real(r8), pointer :: p(:), velocity(:,:)
    real(r8), allocatable :: H(:), T(:)

    call env%timer%start('output')
    allocate(H(this%mesh%ncell_onP), T(this%mesh%ncell_onP))
    call this%solver%get_cell_flow_soln(p, velocity)
    call this%solver%get_cell_heat_soln(H)
    call this%solver%get_cell_temp_soln(T)
    call this%solver%set_temporal_output(this%temporal_output)
    call this%output%write_solution(time, p, velocity, H, T, this%composition%vfrac, this%temporal_output)
    call env%timer%stop('output')
  end subroutine

end module ns_ht_2d_sim_type
