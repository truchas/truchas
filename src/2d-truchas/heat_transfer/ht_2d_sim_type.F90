!!
!! HT_2D_SIM_TYPE
!!
!! This module defines HT_2D_SIM, which owns the mesh, material state,
!! thermal model, solver, time-integration control, and output writer for a
!! two-dimensional thermal transport simulation.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, August 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ht_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use unstr_2d_mesh_type
  use material_database_type
  use material_model_type
  use material_distribution_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use scalar_func_projection
  use ht_2d_model_type
  use ht_2d_solver_type
  use ht_2d_vtkhdf_writer_type
  use simulation_environment_type
  use simulation_type
  implicit none
  private

  type, extends(simulation), public :: ht_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(material_database) :: matl_db
    type(material_model) :: matl_model
    type(material_distribution), pointer :: matl_dist => null()
    type(ht_2d_model), pointer :: model => null()
    type(ht_2d_solver), pointer :: solver => null()
    type(ht_2d_vtkhdf_writer) :: output
    type(parameter_list) :: temporal_output
    !! Integration control
    real(r8) :: t_init
    real(r8), allocatable :: tout(:)
  contains
    final :: ht_2d_sim_delete
    procedure :: init
    procedure :: run
    procedure :: write_solution
  end type

contains

  subroutine ht_2d_sim_delete(this)
    type(ht_2d_sim), intent(inout) :: this
    call this%output%close()
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%matl_dist)) deallocate(this%matl_dist)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine ht_2d_sim_delete


  subroutine init(this, env, params, stat, errmsg)

    use unstr_2d_mesh_factory
    use signal_handler, only: init_signal_handler, SIGURG
    use material_factory, only: load_material_database

    class(ht_2d_sim), intent(out) :: this
    type(simulation_environment), intent(inout) :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: mesh_params, materials_params
    type(parameter_list), pointer :: model_params, solver_params, control_params
    character(:), allocatable :: matl_name(:)
    real(r8), allocatable :: temp(:)
    integer :: i, rlev

    stat = 0
    !! Catch SIGURG signals.
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
    call env%simlog%end_section('Material model complete.')

    !! Construct the material distribution.
    call env%simlog%begin_section('Constructing material distribution.')
    allocate(this%matl_dist)
    if (params%is_sublist('material-regions')) then
      call this%matl_dist%init(env, this%mesh, this%matl_model, params%sublist('material-regions'), rlev, stat, errmsg)
    else
      call this%matl_dist%init_uniform(this%mesh, this%matl_model, 1, stat, errmsg)
    end if
    if (stat /= 0) then
      call env%simlog%end_section('Material distribution construction failed.')
      return
    end if
    call env%simlog%end_section('Material distribution complete.')

    !! Construct the thermal model.
    if (.not.params%is_sublist('ht-model')) then
      stat = 1
      errmsg = 'missing "ht-model" sublist parameter'
      return
    end if
    model_params => params%sublist('ht-model')
    call env%simlog%begin_section('Constructing thermal model.')
    call env%timer%start('ht-model')
    call construct_thermal_model(model_params, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Thermal model construction failed.')
      return
    end if
    call env%simlog%end_section('Thermal model complete.')
    call env%timer%stop('ht-model')

    !! Construct the thermal solver.
    if (.not.params%is_sublist('ht-solver')) then
      stat = 1
      errmsg = 'missing "ht-solver" sublist parameter'
      return
    end if
    solver_params => params%sublist('ht-solver')
    call env%simlog%begin_section('Constructing thermal solver.')
    call env%timer%start('ht-solver')
    allocate(this%solver)
    call this%solver%init(env, this%model, solver_params, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Thermal solver construction failed.')
      return
    end if
    call env%simlog%end_section('Thermal solver complete.')
    call env%timer%stop('ht-solver')

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
    call this%solver%set_initial_state(env, this%t_init, temp, stat, errmsg)
    if (stat /= 0) then
      call env%simlog%end_section('Initial solver-state setup failed.')
      return
    end if
    call env%simlog%end_section('Initial solver state ready.')

    !! Create the output file and write the mesh.
    call this%solver%init_temporal_output(this%temporal_output)
    call env%simlog%info('Opening VTKHDF output.')
    call this%output%open(env, this%mesh, this%matl_model, this%temporal_output, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing VTKHDF output: ' // errmsg
      return
    end if

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


    subroutine construct_thermal_model(model_params, stat, errmsg)

      type(parameter_list), intent(inout) :: model_params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      type(parameter_list), pointer :: constants
      real(r8) :: value

      stat = 0
      if (params%is_sublist('physical-constants')) then
        constants => params%sublist('physical-constants')
        if (constants%is_parameter('stefan-boltzmann')) then
          call constants%get('stefan-boltzmann', value, stat, errmsg)
          if (stat /= 0) return
          call model_params%set('stefan-boltzmann', value)
        end if
        if (constants%is_parameter('absolute-zero')) then
          call constants%get('absolute-zero', value, stat, errmsg)
          if (stat /= 0) return
          call model_params%set('absolute-zero', value)
        end if
      end if
      allocate(this%model)
      call this%model%init(env, this%mesh, this%matl_model, this%matl_dist, model_params, stat, errmsg)

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
      call this%solver%init_time_stepper(params, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // params%path() // ': ' // errmsg
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
      call alloc_scalar_func(params, 'initial-temperature', initial_temp, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing initial-temperature: ' // errmsg
        return
      end if
      allocate(temp(this%mesh%ncell_onP), source=0.0_r8)
      call project_scalar_func_to_cell_centers(this%mesh, initial_temp, temp)

    end subroutine initialize_thermal_state

  end subroutine init


  subroutine run(this, env, stat, errmsg)

    class(ht_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: time, t_write
    character(256) :: line

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

  end subroutine run

  !! Write the current cell solution to the VTKHDF output.

  subroutine write_solution(this, t)

    class(ht_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: t

    real(r8), allocatable :: Hcell(:), Tcell(:)

    allocate(Hcell(this%mesh%ncell_onP), Tcell(this%mesh%ncell_onP))

    call this%solver%get_cell_heat_soln(Hcell)
    call this%solver%get_cell_temp_soln(Tcell)

    call this%solver%set_temporal_output(this%temporal_output)
    call this%output%write_solution(t, Hcell, Tcell, this%matl_model, this%matl_dist%vfrac, &
        this%temporal_output)

  end subroutine write_solution


end module ht_2d_sim_type
