!TODO: finish documentation
!!
!! HT_2D_SIM_TYPE
!!
!! This module defines a derived type that encapsulates a heat transfer
!! simulation.  This drives the time integration and generates the output.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ht_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use material_database_type
  use material_model_type
  use material_composition_type
  use scalar_func_factories
  use scalar_func_projection
  use mfd_2d_disc_type
  use ht_2d_model_type
  use ht_2d_solver_type
  use ht_2d_vtkhdf_output
  use time_step_sync_type
  use parallel_communication
  use truchas_logging_services
  use timer_tree_type
  implicit none
  private

  type, public:: ht_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(mfd_2d_disc), pointer :: disc => null()
    type(material_database) :: matl_db
    type(material_model) :: matl_model
    type(material_composition), pointer :: composition => null()
    type(ht_2d_model), pointer :: model => null()
    type(ht_2d_solver), pointer :: solver => null()
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
  end type ht_2d_sim

contains

  subroutine ht_2d_sim_delete(this)
    type(ht_2d_sim), intent(inout) :: this
    if (associated(this%mesh)) deallocate(this%mesh)
    if (associated(this%disc)) deallocate(this%disc)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%solver)) deallocate(this%solver)
    if (associated(this%composition)) deallocate(this%composition)
  end subroutine ht_2d_sim_delete


  subroutine init(this, params)

    use parameter_list_type
    use unstr_2d_mesh_factory
    use signal_handler, only: init_signal_handler, SIGURG
    use material_factory, only: load_material_database
    use material_utilities, only: add_enthalpy_prop

    class(ht_2d_sim), intent(out) :: this
    type(parameter_list) :: params

    type(parameter_list), pointer :: plist
    class(scalar_func), allocatable :: f
    character(:), allocatable :: errmsg, context, matl_names(:)
    real(r8), allocatable :: temp(:)
    real(r8) :: rel_tol
    integer :: stat, max_itr, rlev
    type(parameter_list_iterator) :: piter

    call start_timer('initialization')
    call TLS_info('Initializing the simulation', TLS_VERB_NOISY)

    !! Catch SIGURG signals.
    call init_signal_handler(SIGURG)

    !! Create the mesh.
    call start_timer('mesh')
    if (params%is_sublist('mesh')) then
      plist => params%sublist('mesh')
      context = 'processing ' // plist%path() // ': '
      this%mesh => new_unstr_2d_mesh(plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
    else
      call TLS_fatal('missing "mesh" sublist parameter')
    end if
    call stop_timer('mesh')

    !! Create the discretization object.
    call start_timer('mfd-discretization')
    allocate(this%disc)
    call this%disc%init(this%mesh)
    call stop_timer('mfd-discretization')

    !! Load the material database and initialize the material model
    !TODO: input name instead of hardwiring it
    if (params%is_sublist('materials')) then
      plist => params%sublist('materials')
      context = 'processing ' // plist%path() // ': '
      call load_material_database(this%matl_db, plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
    else
      call TLS_fatal('missing "materials" sublist parameter')
    end if
    allocate(this%composition)
    if (params%is_sublist('material-regions')) then
      plist => params%sublist('material-regions')
      context = 'processing ' // plist%path() // ': '
      call get_material_region_names(plist, matl_names, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      !! Recursive geometric refinement level. The finest unresolved
      !! subtriangle has a linear size approximately 2**(-rlev) times that
      !! of its initial cell triangle.
      call params%get('material-region-refinement-level', rlev, stat, errmsg, default=6)
      if (stat /= 0) call TLS_fatal('processing material-region-refinement-level: ' // errmsg)
      if (rlev < 0) call TLS_fatal('"material-region-refinement-level" must be >= 0')
    else
      plist => params%sublist('materials')
      piter = parameter_list_iterator(plist, sublists_only=.true.)
      if (piter%count() /= 1) call TLS_fatal('multiple materials require a "material-regions" sublist')
      matl_names = [piter%name()]
    end if
    call this%matl_model%init(matl_names, this%matl_db, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    if (this%matl_model%have_void) call TLS_fatal('2D heat transport does not yet support VOID material regions')
    if (params%is_sublist('material-regions')) then
      plist => params%sublist('material-regions')
      call this%composition%init(this%mesh, this%matl_model, plist, rlev, stat, errmsg)
    else
      call this%composition%init_uniform(this%mesh, this%matl_model, 1, stat, errmsg)
    end if
    if (stat /= 0) call TLS_fatal('initializing material composition: ' // errmsg)

    !! Initialize enthalpy
    call add_enthalpy_prop(this%matl_model, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

    !! Create the heat conduction model.
    call start_timer('ht-model')
    if (params%is_sublist('ht-model')) then
      plist => params%sublist('ht-model')
      context = 'processing ' // plist%path() // ': '
      allocate(this%model)
      call this%model%init(this%disc, this%matl_model, this%composition, plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
    else
      call TLS_fatal('missing "ht-model" sublist parameter')
    end if
    call stop_timer('ht-model')

    !! Create the heat conduction solver.
    call start_timer('ht-solver')
    if (params%is_sublist('ht-solver')) then
      plist => params%sublist('ht-solver')
      allocate(this%solver)
      call this%solver%init(this%model, plist)
    else
      call TLS_fatal('missing "ht-solver" sublist parameter')
    end if
    call stop_timer('ht-solver')

    !! Create output file.
    call this%output%open(this%mesh, stat, errmsg)
    if (stat /= 0) call TLS_fatal('processing VTKHDF output: ' // errmsg)

    !! Simulation control parameters
    if (params%is_sublist('sim-control')) then
      plist => params%sublist('sim-control')
      context = 'processing ' // plist%path() // ': '
      call plist%get('initial-time', this%t_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('initial-time-step', this%dt_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      if (this%dt_init <= 0.0_r8) call TLS_fatal(context//'"initial-time-step" must be > 0.0')
      call plist%get('min-time-step', this%dt_min, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      if (this%dt_min > this%dt_init) call TLS_fatal(context//'require "min-time-step" <= "initial-time-step"')
      call plist%get('output-times', this%tout, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      !TODO: check for strictly increasing values in TOUT, TOUT > t_init, or sort
      !and cull those < t_init.
    else
      call TLS_fatal('missing "sim-control" sublist parameter')
    end if

    this%ts_sync = time_step_sync(4)

    !! Generate the initial temperature field
    call start_timer('initial-state')
    if (params%is_parameter('initial-temperature')) then
      context = 'processing initial-temperature: '
      call alloc_scalar_func(params, 'initial-temperature', f, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
    else
      call TLS_fatal('missing "initial-temperature" parameter')
    end if
    allocate(temp(this%mesh%ncell_onP))
    call project_scalar_func_to_cell_centers(this%mesh, f, this%t_init, temp)

    !! Define the initial heat conduction state
    call this%solver%set_initial_state(this%t_init, this%dt_init, temp)
    call stop_timer('initial-state')

    call stop_timer('initialization')

  end subroutine init


  subroutine run(this, stat, errmsg)

    class(ht_2d_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: t, hnext, t_write
    character(80) :: string(2)

    call start_timer('integration')

    !! Write the initial solution
    t = this%solver%time()
    call this%write_solution(t)
    t_write = t ! keep track of the last write time

    call TLS_info('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    call TLS_info(string(1))

    hnext = this%dt_init; this%tlast = t; this%hlast = hnext
    do n = 1, size(this%tout)
      call integrate(this, this%tout(n), hnext, t, stat, errmsg)
      if (stat < 0 .and. t == t_write) exit
      call this%write_solution(t)
      t_write = t ! keep track of the last write time
      call this%solver%write_metrics(string)
      call TLS_info('')
      call TLS_info(string(1))
      call TLS_info(string(2))
      if (stat /= 0) exit
    end do

    if (stat > 0) then  ! caught a signal
      call TLS_info('')
      call TLS_info(errmsg // ': current solution written, and now terminating ...')
      stat = 0  ! this is a successful return
      deallocate(errmsg)
    else if (stat < 0) then
      call TLS_info('')
      errmsg = 'unrecoverable integration failure: ' // errmsg
      call TLS_info(errmsg)
    end if

    call TLS_info('')
    write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', t
    call TLS_info(string(1))

    call this%output%close()

    call stop_timer('integration')

  end subroutine run

  !! This integrates the system to the target time TOUT.  The final solution
  !! achieved is returned in T and THIS%U.  Nominally this will be at time
  !! TOUT, however it will be at some earlier time when there is a failure in
  !! the time stepping (or a user interrupt).  The input value of HNEXT is the
  !! initial step size to take, and its return value is the suggested next
  !! step size (nominally the value passed to the next call to INTEGRATE).
  !! STAT returns a negative value if a time stepping failure occurs, with an
  !! explanatory message in ERRMSG.  If the process recieves the SIGURG
  !! signal, the procedure returns at the end of the next time step with a
  !! positive value of STAT, with T and THIS%U set accordingly.

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
  !! current time to time T and returns the solution in THIS%U.  However if
  !! the step attempt fails, the procedure will re-attempt with successively
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

    tlast = this%solver%time()

    do n = 1, this%max_try
      call this%solver%step(t, hnext, stat)
      if (stat == 0) then ! success
        call this%solver%commit_pending_state
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

  !! This auxiliary subroutine writes the solution to a GMV format viz file.

  subroutine write_solution(this, t)

    class(ht_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: t

    real(r8), allocatable :: Hcell(:), Tcell(:)

    call start_timer('output')

    allocate(Hcell(this%mesh%ncell_onP), Tcell(this%mesh%ncell_onP))

    call this%solver%get_cell_heat_soln(Hcell)
    call this%solver%get_cell_temp_soln(Tcell)

    call this%output%write_solution(t, Hcell, Tcell)

    call stop_timer('output')

  end subroutine write_solution


end module ht_2d_sim_type
