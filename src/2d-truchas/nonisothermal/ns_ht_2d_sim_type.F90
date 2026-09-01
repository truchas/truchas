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
  use material_distribution_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use scalar_func_projection
  use vector_func_class
  use vector_func_factories, only: alloc_vector_func
  use vector_func_projection
  use flow_2d_model_type
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
    type(material_distribution), pointer :: matl_dist => null()
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
    if (associated(this%matl_dist)) deallocate(this%matl_dist)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine

  subroutine init(this, env, params, stat, errmsg)
    use material_factory, only: load_material_database
    use material_utilities, only: add_enthalpy_prop
    use signal_handler, only: init_signal_handler, SIGURG

    class(ns_ht_2d_sim), intent(out) :: this
    type(simulation_environment), intent(inout) :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: mesh_params, model_params, solver_params, materials_params
    type(parameter_list), pointer :: control_params
    character(:), allocatable :: matl_name(:)
    real(r8), allocatable :: temp(:), velocity(:,:)
    integer :: i, rlev
    logical :: inviscid

    stat = 0

    call init_signal_handler(SIGURG)

    !! Construct the mesh.
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

    !! Load the material database.
    if (.not.params%is_sublist('materials')) then
      stat = 1
      errmsg = 'missing "materials" sublist parameter'
      return
    end if
    materials_params => params%sublist('materials')
    call env%simlog%info('Loading material database.')
    call load_material_database(this%matl_db, materials_params, stat, errmsg)
    if (stat /= 0) return

    !! Identify the simulation materials.
    call env%simlog%begin_section('Reading material regions.')
    call read_material_regions(stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Material-region processing failed.')
      return
    end if
    call env%simlog%end_section('Material regions read.')

    !! Construct the material model.
    call env%simlog%begin_section('Constructing material model.')
    call this%matl_model%init(matl_name, this%matl_db, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Material model construction failed.')
      return
    end if
    do i = 1, size(matl_name)
      call env%simlog%info('Using material "' // trim(matl_name(i)) // '".')
    end do
    call env%simlog%info('Adding enthalpy property.')
    call add_enthalpy_prop(this%matl_model, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Material model construction failed.')
      return
    end if
    call env%simlog%end_section('Material model complete.')

    !! Construct the material distribution.
    call env%simlog%begin_section('Constructing material distribution.')
    allocate(this%matl_dist)
    if (params%is_sublist('material-regions')) then
      call this%matl_dist%init(env, this%mesh, this%matl_model, params%sublist('material-regions'), rlev, stat, errmsg)
    else
      call env%simlog%info('Assigning material "' // trim(matl_name(1)) // '" uniformly.')
      call this%matl_dist%init_uniform(this%mesh, this%matl_model, 1, stat, errmsg)
    end if
    if (stat /= 0) then
      call env%simlog%end_section('Material distribution construction failed.')
      return
    end if
    call env%simlog%end_section('Material distribution complete.')

    !! Construct the Navier-Stokes flow model.
    if (.not.params%is_sublist('flow-model')) then
      stat = 1
      errmsg = 'missing "flow-model" sublist parameter'
      return
    end if
    model_params => params%sublist('flow-model')
    call env%simlog%begin_section('Constructing flow model.')
    call construct_flow_model(model_params, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Flow model construction failed.')
      return
    end if
    call env%simlog%end_section('Flow model complete.')

    !! Construct the thermal model.
    if (.not.params%is_sublist('ht-model')) then
      stat = 1
      errmsg = 'missing "ht-model" sublist parameter'
      return
    end if
    call env%simlog%begin_section('Constructing thermal model.')
    call construct_thermal_model(params%sublist('ht-model'), stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Thermal model construction failed.')
      return
    end if
    call env%simlog%end_section('Thermal model complete.')

    !! Construct the coupled solver.
    if (.not.params%is_sublist('solver')) then
      stat = 1
      errmsg = 'missing "solver" sublist'
      return
    end if
    solver_params => params%sublist('solver')
    call env%simlog%begin_section('Constructing coupled solver.')
    allocate(this%solver)
    call this%solver%init(env, this%flow_model, this%ht_model, this%matl_model, this%matl_dist, &
        solver_params, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Coupled solver construction failed.')
      return
    end if
    call env%simlog%end_section('Coupled solver complete.')

    !! Configure integration.
    if (.not.params%is_sublist('sim-control')) then
      stat = 1
      errmsg = 'missing "sim-control" sublist parameter'
      return
    end if
    control_params => params%sublist('sim-control')
    call env%simlog%begin_section('Configuring integration.')
    call configure_integration(control_params, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Integration configuration failed.')
      return
    end if
    call env%simlog%end_section('Integration configuration complete.')

    !! Initialize the flow state.
    call env%simlog%begin_section('Initializing flow state.')
    call initialize_flow_state(velocity, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Flow-state initialization failed.')
      return
    end if
    call env%simlog%end_section('Flow-state initialization complete.')

    !! Initialize the thermal state.
    call env%simlog%begin_section('Initializing thermal state.')
    call initialize_thermal_state(temp, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Thermal-state initialization failed.')
      return
    end if
    call env%simlog%end_section('Thermal-state initialization complete.')

    !! Set the initial solver state.
    call env%simlog%begin_section('Setting initial solver state.')
    call this%solver%set_initial_state(env, this%t_init, velocity, temp, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Initial solver-state setup failed.')
      return
    end if
    call env%simlog%end_section('Initial solver state ready.')

    !! Create the output file and write the mesh.
    call this%solver%init_temporal_output(this%temporal_output)
    call env%simlog%info('Opening VTKHDF output.')
    call this%output%open(env, this%mesh, this%matl_model, this%temporal_output, stat, errmsg)
    if (stat /= 0) errmsg = 'opening VTKHDF output: ' // errmsg

  contains

    subroutine read_material_regions(stat, errmsg)

      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      type(parameter_list_iterator) :: piter
      character(:), allocatable :: region_name(:), region_matl_name(:)
      integer :: i

      stat = 0
      if (params%is_sublist('material-regions')) then
        call get_material_region_names(params%sublist('material-regions'), matl_name, stat, errmsg, &
            region_name, region_matl_name)
        if (stat /= 0) return
        do i = 1, size(region_name)
          call env%simlog%info('Region "' // trim(region_name(i)) // '": material="' // &
              trim(region_matl_name(i)) // '".')
        end do
        call params%get('material-region-refinement-level', rlev, default=6, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        if (rlev < 0) then
          stat = 1
          errmsg = '"material-region-refinement-level" must be >= 0'
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

    end subroutine read_material_regions


    subroutine construct_flow_model(params, stat, errmsg)

      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      type(parameter_list), pointer :: bc_params
      real(r8), allocatable :: body_acceleration(:)
      character(96) :: message

      stat = 0
      call params%get('inviscid', inviscid, default=.false., stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // params%path() // ': ' // errmsg
        return
      end if
      call params%get('body-acceleration', body_acceleration, stat=stat, errmsg=errmsg, &
          default=[0.0_r8, 0.0_r8])
      if (stat /= 0) then
        errmsg = 'processing ' // params%path() // ': ' // errmsg
        return
      end if
      if (.not.params%is_sublist('bc')) then
        stat = 1
        errmsg = 'missing "bc" sublist parameter in ' // params%path()
        return
      end if
      bc_params => params%sublist('bc')

      allocate(this%flow_model)
      if (inviscid) then
        call env%simlog%info('Using inviscid flow.')
      else
        call env%simlog%info('Using viscous flow.')
      end if
      if (any(body_acceleration /= 0.0_r8)) then
        write(message, '(a,es11.4,a,es11.4,a)') 'Using body acceleration [', body_acceleration(1), ', ', &
            body_acceleration(2), '].'
        call env%simlog%info(trim(message))
      end if

      call this%flow_model%init_core(env, this%mesh, bc_params, stat, errmsg, &
          body_acceleration=body_acceleration, inviscid=inviscid)
      if (stat /= 0) errmsg = 'processing ' // bc_params%path() // ': ' // errmsg

    end subroutine construct_flow_model


    subroutine construct_thermal_model(params, stat, errmsg)

      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      allocate(this%ht_model)
      call this%ht_model%init(env, this%mesh, this%matl_model, this%matl_dist, params, stat, errmsg, advection=.true.)

    end subroutine construct_thermal_model


    subroutine configure_integration(params, stat, errmsg)

      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      character(96) :: message

      stat = 0
      call params%get('initial-time', this%t_init, default=0.0_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // params%path() // ': ' // errmsg
        return
      end if
      call params%get('output-times', this%tout, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // params%path() // ': ' // errmsg
        return
      end if
      if (size(this%tout) == 0) then
        stat = 1
        errmsg = 'processing ' // params%path() // ': require nonempty output-times'
        return
      end if
      if (any(this%tout <= this%t_init) .or. any(this%tout(2:) <= this%tout(:size(this%tout)-1))) then
        stat = 1
        errmsg = 'processing ' // params%path() // ': require strictly increasing output-times after initial-time'
        return
      end if
      write(message,'(a,es11.4)') 'Initial time: ', this%t_init
      call env%simlog%info(trim(message))
      write(message,'(a,es11.4)') 'Initial time step: ', this%solver%initial_time_step()
      call env%simlog%info(trim(message))
      if (size(this%tout) == 1) then
        write(message,'(a,es11.4)') 'Output schedule: 1 time at ', this%tout(1)
      else
        write(message,'(a,i0,a,es11.4,a,es11.4)') 'Output schedule: ', size(this%tout), &
            ' times from ', this%tout(1), ' through ', this%tout(size(this%tout))
      end if
      call env%simlog%info(trim(message))

    end subroutine configure_integration


    subroutine initialize_flow_state(velocity, stat, errmsg)

      real(r8), allocatable, intent(out) :: velocity(:,:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      class(vector_func), allocatable :: initial_velocity_func

      stat = 0
      allocate(velocity(2,this%mesh%ncell_onP), source=0.0_r8)
      if (params%is_parameter('initial-velocity')) then
        call env%simlog%info('Projecting user-supplied initial velocity.')
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
        call project_vector_func_to_cell_centers(this%mesh, initial_velocity_func, velocity)
      else
        call env%simlog%info('Using zero initial velocity.')
      end if

    end subroutine initialize_flow_state


    subroutine initialize_thermal_state(temp, stat, errmsg)

      real(r8), allocatable, intent(out) :: temp(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      class(scalar_func), allocatable :: initial_temp

      stat = 0
      if (.not.params%is_parameter('initial-temperature')) then
        stat = 1
        errmsg = 'missing "initial-temperature" parameter'
        return
      end if
      call env%simlog%info('Projecting user-supplied initial temperature.')
      allocate(temp(this%mesh%ncell_onP), source=0.0_r8)
      call alloc_scalar_func(params, 'initial-temperature', initial_temp, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing initial-temperature: ' // errmsg
        return
      end if
      call project_scalar_func_to_cell_centers(this%mesh, initial_temp, temp)

    end subroutine initialize_thermal_state

  end subroutine

  subroutine run(this, env, stat, errmsg)
    class(ns_ht_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    character(256) :: line
    real(r8) :: time, t_write

    stat = 0
    time = this%solver%last_time()
    write(line,'(a,es0.5)') 'integration-begin t0=', time
    call env%simlog%info('')
    call env%simlog%begin_section(trim(line))
    call env%timer%start('integration')
    call env%timer%start('output')
    call this%write_solution(time)
    call env%timer%stop('output')
    write(line,'(a,es0.5)') 'output t=', time
    call env%simlog%info(trim(line))
    t_write = time
    do n = 1, size(this%tout)
      call this%solver%integrate(env, this%tout(n), stat, errmsg)
      time = this%solver%last_time()
      if (stat < 0 .and. time == t_write) exit
      call env%timer%start('output')
      call this%write_solution(time)
      call env%timer%stop('output')
      write(line,'(a,es0.5)') 'output t=', time
      call env%simlog%info(trim(line))
      t_write = time
      if (stat /= 0) exit
    end do
    if (stat > 0) then
      stat = 0
      deallocate(errmsg)
      write(line,'(a,es0.5,a,i0,a)') 'integration-end t=', time, ' nstep=', this%solver%num_steps(), &
          ' status=interrupted reason=signal'
    else if (stat == 0) then
      write(line,'(a,es0.5,a,i0,a)') 'integration-end t=', time, ' nstep=', this%solver%num_steps(), &
          ' status=complete'
    else
      write(line,'(a,es0.5,a,i0,a)') 'integration-end t=', time, ' nstep=', this%solver%num_steps(), &
          ' status=failed'
    end if
    call this%output%close()
    call env%simlog%end_section(trim(line))
    call env%timer%stop('integration')
  end subroutine


  subroutine write_solution(this, time)
    class(ns_ht_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: time

    real(r8), pointer :: p(:), velocity(:,:)
    real(r8), allocatable :: H(:), T(:)

    allocate(H(this%mesh%ncell_onP), T(this%mesh%ncell_onP))
    call this%solver%get_cell_flow_soln(p, velocity)
    call this%solver%get_cell_heat_soln(H)
    call this%solver%get_cell_temp_soln(T)
    call this%solver%set_temporal_output(this%temporal_output)
    call this%output%write_solution(time, p, velocity, H, T, this%matl_dist%vfrac, this%temporal_output)
  end subroutine

end module ns_ht_2d_sim_type
