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
  use unstr_2d_mesh_type
  use material_database_type
  use material_model_type
  use material_composition_type
  use scalar_func_factories
  use scalar_func_projection
  use ht_2d_model_type
  use ht_2d_solver_type
  use ht_2d_vtkhdf_writer_type
  use time_step_sync_type
  use parallel_communication
  use simulation_environment_type
  use simulation_log_type, only: LOG_DETAIL
  implicit none
  private

  type, public :: ht_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(material_database) :: matl_db
    type(material_model) :: matl_model
    type(material_composition), pointer :: composition => null()
    type(ht_2d_model), pointer :: model => null()
    type(ht_2d_solver), pointer :: solver => null()
    integer :: integrator_log_unit = 0
    type(ht_2d_vtkhdf_writer) :: output
    !! Integration control
    real(r8) :: t_init
    real(r8) :: tlast, hlast
    real(r8) :: dt_init, dt_min, dt_max
    integer  :: max_try
    real(r8), allocatable :: tout(:)
    type(time_step_sync) :: ts_sync
  contains
    final :: ht_2d_sim_delete
    procedure :: init
    procedure :: run
    procedure :: write_solution
  end type

contains

  subroutine ht_2d_sim_delete(this)
    type(ht_2d_sim), intent(inout) :: this
    if (this%integrator_log_unit /= 0) close(this%integrator_log_unit)
    call this%output%close()
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%composition)) deallocate(this%composition)
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine ht_2d_sim_delete


  subroutine init(this, env, params, stat, errmsg)

    use parameter_list_type
    use unstr_2d_mesh_factory
    use signal_handler, only: init_signal_handler, SIGURG
    use material_factory, only: load_material_database
    use material_utilities, only: add_enthalpy_prop

    class(ht_2d_sim), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist, integrator_params
    class(scalar_func), allocatable :: f
    character(:), allocatable :: context, matl_names(:)
    character(:), allocatable :: integrator_log_file, integrator_log_path
    character(256) :: iomsg
    real(r8), allocatable :: temp(:)
    integer :: rlev
    type(parameter_list_iterator) :: piter

    stat = 0
    call env%timer%start('initialization')
    call env%simlog%info('Initializing the simulation', LOG_DETAIL)

    !! Catch SIGURG signals.
    call init_signal_handler(SIGURG)

    !! Create the mesh.
    call env%timer%start('mesh')
    if (params%is_sublist('mesh')) then
      plist => params%sublist('mesh')
      context = 'processing ' // plist%path() // ': '
      this%mesh => new_unstr_2d_mesh(env, plist, stat, errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
    else
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
      return
    end if
    call env%timer%stop('mesh')

    !! Load the material database and initialize the material model
    if (params%is_sublist('materials')) then
      plist => params%sublist('materials')
      context = 'processing ' // plist%path() // ': '
      call load_material_database(this%matl_db, plist, stat, errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
    else
      stat = 1
      errmsg = 'missing "materials" sublist parameter'
      return
    end if
    allocate(this%composition)
    if (params%is_sublist('material-regions')) then
      plist => params%sublist('material-regions')
      context = 'processing ' // plist%path() // ': '
      call get_material_region_names(plist, matl_names, stat, errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
      !! Recursive geometric refinement level. The finest unresolved
      !! subtriangle has a linear size approximately 2**(-rlev) times that
      !! of its initial cell triangle.
      call params%get('material-region-refinement-level', rlev, stat, errmsg, default=6)
      if (stat /= 0) then
        errmsg = 'processing material-region-refinement-level: ' // errmsg
        return
      else if (rlev < 0) then
        stat = 1
        errmsg = '"material-region-refinement-level" must be >= 0'
        return
      end if
    else
      plist => params%sublist('materials')
      piter = parameter_list_iterator(plist, sublists_only=.true.)
      if (piter%count() /= 1) then
        stat = 1
        errmsg = 'multiple materials require a "material-regions" sublist'
        return
      end if
      matl_names = [piter%name()]
    end if
    call this%matl_model%init(matl_names, this%matl_db, stat, errmsg)
    if (stat /= 0) return
    if (this%matl_model%have_void) then
      stat = 1
      errmsg = '2D heat transport does not yet support VOID material regions'
      return
    end if
    if (params%is_sublist('material-regions')) then
      plist => params%sublist('material-regions')
      call this%composition%init(this%mesh, this%matl_model, plist, rlev, stat, errmsg)
    else
      call this%composition%init_uniform(this%mesh, this%matl_model, 1, stat, errmsg)
    end if
    if (stat /= 0) then
      errmsg = 'initializing material composition: ' // errmsg
      return
    end if

    !! Initialize enthalpy
    call add_enthalpy_prop(this%matl_model, stat, errmsg)
    if (stat /= 0) return

    !! Create the heat conduction model.
    call env%timer%start('ht-model')
    if (params%is_sublist('ht-model')) then
      plist => params%sublist('ht-model')
      context = 'processing ' // plist%path() // ': '
      allocate(this%model)
      call this%model%init(env, this%mesh, this%matl_model, this%composition, plist, stat, errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
    else
      stat = 1
      errmsg = 'missing "ht-model" sublist parameter'
      return
    end if
    call env%timer%stop('ht-model')

    !! Create the heat conduction solver.
    call env%timer%start('ht-solver')
    if (params%is_sublist('ht-solver')) then
      plist => params%sublist('ht-solver')
      allocate(this%solver)
      context = 'processing ' // plist%path() // ': '
      call this%solver%init(env, this%model, plist, stat, errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
      if (plist%is_sublist('integrator')) then
        integrator_params => plist%sublist('integrator')
        if (integrator_params%is_parameter('integrator-log-file')) then
          context = 'processing ' // integrator_params%path() // ': '
          call integrator_params%get('integrator-log-file', integrator_log_file, stat, errmsg)
          if (stat /= 0) then
            errmsg = context // errmsg
            return
          end if
          if (len_trim(integrator_log_file) == 0) then
            stat = 1
            errmsg = context//'"integrator-log-file" must not be empty'
            return
          end if
          if (integrator_log_file(1:1) == '/') then
            integrator_log_path = trim(integrator_log_file)
          else
            integrator_log_path = trim(env%output_dir)//'/'//trim(integrator_log_file)
          end if
          if (is_IOP) then
            open(newunit=this%integrator_log_unit, file=integrator_log_path, status='replace', action='write', &
                iostat=stat, iomsg=iomsg)
            if (stat /= 0) errmsg = context // trim(iomsg)
          end if
          call broadcast_iop_status(stat, errmsg)
          if (stat /= 0) return
          if (is_IOP) then
            call this%solver%set_integrator_log(this%integrator_log_unit)
          end if
        end if
      end if
    else
      stat = 1
      errmsg = 'missing "ht-solver" sublist parameter'
      return
    end if
    call env%timer%stop('ht-solver')

    !! Create output file.
    call this%output%open(env, this%mesh, this%matl_model, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing VTKHDF output: ' // errmsg
      return
    end if

    !! Simulation control parameters
    if (params%is_sublist('sim-control')) then
      plist => params%sublist('sim-control')
      context = 'processing ' // plist%path() // ': '
      call plist%get('initial-time', this%t_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
      call plist%get('initial-time-step', this%dt_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      else if (this%dt_init <= 0.0_r8) then
        stat = 1
        errmsg = context//'"initial-time-step" must be > 0.0'
        return
      end if
      call plist%get('min-time-step', this%dt_min, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
      call plist%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
      call plist%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      else if (this%dt_min > this%dt_init) then
        stat = 1
        errmsg = context//'require "min-time-step" <= "initial-time-step"'
        return
      end if
      call plist%get('output-times', this%tout, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
      !TODO: check for strictly increasing values in TOUT, TOUT > t_init, or sort
      !and cull those < t_init.
    else
      stat = 1
      errmsg = 'missing "sim-control" sublist parameter'
      return
    end if

    this%ts_sync = time_step_sync(4)

    !! Generate the initial temperature field
    call env%timer%start('initial-state')
    if (params%is_parameter('initial-temperature')) then
      context = 'processing initial-temperature: '
      call alloc_scalar_func(params, 'initial-temperature', f, stat, errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
    else
      stat = 1
      errmsg = 'missing "initial-temperature" parameter'
      return
    end if
    allocate(temp(this%mesh%ncell_onP))
    call project_scalar_func_to_cell_centers(this%mesh, f, temp)

    !! Define the initial heat conduction state
    call this%solver%set_initial_state(env, this%t_init, this%dt_init, temp, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'initializing thermal state: ' // errmsg
      return
    end if
    call env%timer%stop('initial-state')

    call env%timer%stop('initialization')

  end subroutine init


  subroutine broadcast_iop_status(stat, errmsg)
    integer, intent(inout) :: stat
    character(:), allocatable, intent(inout) :: errmsg
    integer :: n

    call broadcast(stat)
    if (stat /= 0) then
      if (is_IOP) n = len(errmsg)
      call broadcast(n)
      if (.not.is_IOP) allocate(character(len=n) :: errmsg)
      call broadcast(errmsg)
    end if
  end subroutine broadcast_iop_status


  subroutine run(this, env, stat, errmsg)

    class(ht_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: t, hnext, t_write
    character(80) :: string(2)

    call env%timer%start('integration')

    !! Write the initial solution
    t = this%solver%last_time()
    call this%write_solution(env, t)
    t_write = t ! keep track of the last write time

    call env%simlog%info('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    call env%simlog%info(string(1))

    hnext = this%dt_init; this%tlast = t; this%hlast = hnext
    do n = 1, size(this%tout)
      call integrate(this, this%tout(n), hnext, t, stat, errmsg)
      if (stat < 0 .and. t == t_write) exit
      call this%write_solution(env, t)
      t_write = t ! keep track of the last write time
      call this%solver%write_metrics(string)
      call env%simlog%info('')
      call env%simlog%info(string(1))
      call env%simlog%info(string(2))
      if (stat /= 0) exit
    end do

    if (stat > 0) then  ! caught a signal
      call env%simlog%info('')
      call env%simlog%info(errmsg // ': current solution written, and now terminating ...')
      stat = 0  ! this is a successful return
      deallocate(errmsg)
    else if (stat < 0) then
      call env%simlog%info('')
      errmsg = 'unrecoverable integration failure: ' // errmsg
      call env%simlog%info(errmsg)
    end if

    call env%simlog%info('')
    write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', t
    call env%simlog%info(string(1))

    call this%output%close()

    call env%timer%stop('integration')

  end subroutine run

  !! Integrate the system toward TOUT. The final time reached is returned in
  !! T, and the solver retains the corresponding committed thermal state.
  !! Nominally T equals TOUT, but it is earlier after a time-stepping failure
  !! or user interrupt. The input value of HNEXT is the initial step size;
  !! its return value is the suggested next step size.
  !! STAT returns a negative value if a time stepping failure occurs, with an
  !! explanatory message in ERRMSG.  If the process recieves the SIGURG
  !! signal, the procedure returns at the end of the next time step with a
  !! positive value of STAT, with T and the solver's committed state set
  !! accordingly.

  subroutine integrate(this, tout, hnext, t, stat, errmsg)

    use signal_handler, only: read_signal, SIGURG

    class(ht_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: tout
    real(r8), intent(inout) :: hnext
    real(r8), intent(out) :: t
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical :: sig_rcvd

    do
      !! Time for next step; nominally TLAST+HNEXT but possibly adjusted
      t = this%ts_sync%next_time(tout, this%tlast, this%hlast, hnext)
      call step(this, t, hnext, stat, errmsg)
      if (stat /= 0) then
        !! Return the last good solution
        t = this%tlast
        return
      end if
      this%hlast = t - this%tlast
      this%tlast = t

      hnext = min(hnext, this%dt_max)

      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if

      if (t == tout) return
    end do

  end subroutine integrate

  !! Take a resilient step.  Nominally this takes a single step from the
  !! current time to time T. However if the step attempt fails, the procedure
  !! will re-attempt with successively
  !! smaller step sizes, until the step is successful or the number of
  !! attempts exceeds a maximum or the step size gets too small.  Thus the
  !! return value of T may differ from its input value.  HNEXT returns the
  !! suggested next time step.  STAT returns a negative value if the step
  !! was ultimately unsuccessful, with an explanatory message in ERRMSG.

  subroutine step(this, t,  hnext, stat, errmsg)

    class(ht_2d_sim), intent(inout) :: this
    real(r8), intent(inout) :: t
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: tlast

    tlast = this%solver%last_time()

    do n = 1, this%max_try
      call this%solver%step(t, hnext, stat)
      if (stat == 0) then ! success
        call this%solver%commit_step
        return
      end if
      t = tlast + hnext
      if (t - tlast < this%dt_min) then
        stat = -1
        errmsg = 'next time step is too small'
        return
      end if
    end do

    stat = -2
    errmsg = 'unable to take a time step'

  end subroutine step

  !! Write the current cell solution to the VTKHDF output.

  subroutine write_solution(this, env, t)

    class(ht_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t

    real(r8), allocatable :: Hcell(:), Tcell(:)

    call env%timer%start('output')

    allocate(Hcell(this%mesh%ncell_onP), Tcell(this%mesh%ncell_onP))

    call this%solver%get_cell_heat_soln(Hcell)
    call this%solver%get_cell_temp_soln(Tcell)

    call this%output%write_solution(t, Hcell, Tcell, this%composition%vfrac)

    call env%timer%stop('output')

  end subroutine write_solution


end module ht_2d_sim_type
