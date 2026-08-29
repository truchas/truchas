program test_ns_ht_2d_solver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use simulation_environment_type
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use material_database_type
  use material_model_type
  use material_composition_type
  use material_factory, only: load_material_database
  use material_utilities, only: add_enthalpy_prop
  use flow_2d_model_type
  use ht_2d_model_type
  use ht_2d_solver_type
  use ns_ht_2d_solver_type
  use time_step_sync_type
  implicit none

  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  call init_parallel_communication
  call fhypre_initialize
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_SILENT)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  allocate(env%timer)
  call env%simlog%init(env%comm, 'test_ns_ht_2d_solver.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call fail('initializing simulation log: ' // errmsg)

  status = 0
  call test_zero_velocity_thermal_equivalence

  call halt_parallel_communication
  stop status

contains

  subroutine test_zero_velocity_thermal_equivalence
    type(unstr_2d_mesh), pointer :: mesh
    type(material_database) :: database
    type(material_model) :: matl_model
    type(material_composition), target :: composition
    type(flow_2d_model), target :: flow_model
    type(ht_2d_model), target :: ht_model
    type(ht_2d_model), target :: standalone_ht_model
    type(ns_ht_2d_solver), target :: solver
    type(ht_2d_solver), target :: thermal_solver
    type(parameter_list), pointer :: matl_params, flow_bc_params, ht_params, standalone_ht_params
    type(parameter_list), target :: solver_params
    type(parameter_list), pointer :: plist, thermal_params
    real(r8), allocatable :: velocity(:,:), temp(:), temp_result(:), heat_result(:), thermal_temp(:), thermal_heat(:)
    integer :: n, ntry
    real(r8), parameter :: equivalence_tol = 1.0e-12_r8
    real(r8) :: hlast, hnext, thermal_hnext, t, t_np1
    type(time_step_sync) :: ts_sync

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call parameter_list_from_json_string( &
      '{"liquid":{"properties":{"fluid":true,"conductivity":1.0,"density":1.0,"specific-heat":1.0}}}', &
      matl_params, errmsg)
    if (.not.associated(matl_params)) call fail('parsing material input: ' // errmsg)
    call load_material_database(database, matl_params, stat, errmsg)
    if (stat /= 0) call fail('loading material database: ' // errmsg)
    call matl_model%init(['liquid'], database, stat, errmsg)
    if (stat /= 0) call fail('initializing material model: ' // errmsg)
    call composition%init_uniform(mesh, matl_model, 1, stat, errmsg)
    if (stat /= 0) call fail('initializing material composition: ' // errmsg)
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call fail('adding enthalpy property: ' // errmsg)

    call parameter_list_from_json_string( &
      '{"wall":{"type":"no-slip","face-set-ids":[1,2,3,4]}}', flow_bc_params, errmsg)
    if (.not.associated(flow_bc_params)) call fail('parsing flow BC input: ' // errmsg)
    call flow_model%init(env, mesh, flow_bc_params, [1.0_r8], 1.0_r8, stat, errmsg)
    if (stat /= 0) call fail('initializing flow model: ' // errmsg)

    call parameter_list_from_json_string( &
      '{"bc":{"all":{"type":"temperature","face-set-ids":[1,2,3,4],"temp":2.0}},"source":{}}', &
      ht_params, errmsg)
    if (.not.associated(ht_params)) call fail('parsing thermal-model input: ' // errmsg)
    call ht_model%init(env, mesh, matl_model, composition, ht_params, stat, errmsg, advection=.true.)
    if (stat /= 0) call fail('initializing thermal model: ' // errmsg)
    call parameter_list_from_json_string( &
      '{"bc":{"all":{"type":"temperature","face-set-ids":[1,2,3,4],"temp":2.0}},"source":{}}', &
      standalone_ht_params, errmsg)
    if (.not.associated(standalone_ht_params)) call fail('parsing standalone thermal-model input: ' // errmsg)
    call standalone_ht_model%init(env, mesh, matl_model, composition, standalone_ht_params, stat, errmsg)
    if (stat /= 0) call fail('initializing standalone thermal model: ' // errmsg)

    plist => solver_params%sublist('flow')
    plist => plist%sublist('momentum-solver')
    call set_flow_solver_params(plist)
    plist => solver_params%sublist('flow')
    plist => plist%sublist('projection-solver')
    call set_flow_solver_params(plist)
    thermal_params => solver_params%sublist('thermal')
    call set_thermal_solver_params(thermal_params)
    call solver_params%set('initial-time-step', 1.0e-3_r8)
    call solver_params%set('min-time-step', 1.0e-8_r8)
    call solver_params%set('max-time-step', 1.0e-3_r8)
    call solver_params%set('time-step-growth', 1.0_r8)
    call solver_params%set('courant-number', 0.5_r8)
    call solver_params%set('max-try-at-step', 4)
    call thermal_solver%init(env, standalone_ht_model, thermal_params, stat, errmsg)
    if (stat /= 0) call fail('initializing standalone thermal solver: ' // errmsg)
    call solver%init(env, flow_model, ht_model, matl_model, composition, &
        solver_params, stat, errmsg)
    if (stat /= 0) call fail('initializing coupled solver: ' // errmsg)

    call mesh%init_cell_centroid
    allocate(velocity(2,mesh%ncell_onP), temp(mesh%ncell_onP), temp_result(mesh%ncell_onP), &
        heat_result(mesh%ncell_onP), thermal_temp(mesh%ncell_onP), thermal_heat(mesh%ncell_onP))
    velocity = 0.0_r8
    temp = 2.0_r8 + sin(acos(-1.0_r8)*mesh%cell_centroid(1,:mesh%ncell_onP)) * &
        sin(acos(-1.0_r8)*mesh%cell_centroid(2,:mesh%ncell_onP))
    call thermal_solver%set_initial_state(env, 0.0_r8, 1.0e-3_r8, temp, stat, errmsg)
    if (stat /= 0) call fail('initializing standalone thermal state: ' // errmsg)
    call solver%set_initial_state(env, 0.0_r8, velocity, temp, stat, errmsg)
    if (stat /= 0) call fail('initializing coupled state: ' // errmsg)
    call require_zero_face_velocity(solver, mesh%nface, 'zero initial velocity was not preserved')
    call thermal_solver%get_cell_temp_soln(thermal_temp)
    call thermal_solver%get_cell_heat_soln(thermal_heat)
    call thermal_solver%step(1.0e-3_r8, thermal_hnext, stat)
    call require(stat == 0, 'thermal step for reject test failed')
    if (stat /= 0) return
    call thermal_solver%reject_step()
    call thermal_solver%get_cell_temp_soln(temp_result)
    call thermal_solver%get_cell_heat_soln(heat_result)
    call require(maxval(abs(temp_result - thermal_temp)) == 0.0_r8, &
        'thermal temperature was not restored by reject_step')
    call require(maxval(abs(heat_result - thermal_heat)) == 0.0_r8, &
        'thermal enthalpy was not restored by reject_step')
    ts_sync = time_step_sync(4)
    t = 0.0_r8
    hlast = 1.0e-3_r8
    hnext = hlast
    do n = 1, 10
      do while (t < real(n,r8)*1.0e-3_r8)
        t_np1 = ts_sync%next_time(real(n,r8)*1.0e-3_r8, t, hlast, hnext)
        do ntry = 1, 4
          call thermal_solver%step(t_np1, thermal_hnext, stat)
          if (stat == 0) exit
          t_np1 = t + thermal_hnext
          if (t_np1 - t < 1.0e-8_r8) then
            call fail('standalone thermal next step is too small')
          end if
        end do
        call require(stat == 0, 'standalone thermal step failed')
        if (stat /= 0) return
        call thermal_solver%commit_step
        hlast = t_np1 - t
        hnext = min(thermal_hnext, hlast, 1.0e-3_r8)
        t = t_np1
      end do
      call solver%integrate(env, real(n,r8)*1.0e-3_r8, stat, errmsg)
      call require(stat == 0, 'zero-velocity coupled integration failed')
      if (stat /= 0) return
      call require_zero_face_velocity(solver, mesh%nface, 'zero velocity did not remain zero')
      call thermal_solver%get_cell_temp_soln(thermal_temp)
      call thermal_solver%get_cell_heat_soln(thermal_heat)
      call solver%get_cell_temp_soln(temp_result)
      call solver%get_cell_heat_soln(heat_result)
      call require(maxval(abs(temp_result - thermal_temp)) < equivalence_tol, &
          'zero-velocity coupled and standalone thermal temperatures differ')
      call require(maxval(abs(heat_result - thermal_heat)) < equivalence_tol, &
          'zero-velocity coupled and standalone thermal enthalpies differ')
    end do
    call require(t == 1.0e-2_r8 .and. solver%last_time() == t, &
        'coupled step did not reach its requested endpoint')
  end subroutine


  subroutine require_zero_face_velocity(solver, nface, message)
    type(ns_ht_2d_solver), target, intent(in) :: solver
    integer, intent(in) :: nface
    character(*), intent(in) :: message
    real(r8), pointer :: velocity(:)

    call solver%get_face_velocity(velocity)
    call require(all(velocity == 0.0_r8), message)
  end subroutine


  subroutine set_flow_solver_params(params)
    type(parameter_list), intent(inout) :: params

    call params%set('rel-tol', 1.0e-10_r8)
    call params%set('max-ds-iter', 100)
    call params%set('max-amg-iter', 100)
  end subroutine


  subroutine set_thermal_solver_params(params)
    type(parameter_list), intent(inout) :: params
    type(parameter_list), pointer :: plist

    plist => params%sublist('preconditioner')
    call plist%set('method', 'BoomerAMG')
    plist => plist%sublist('params')
    call plist%set('num-cycles', 1)
    plist => params%sublist('error-norm')
    call plist%set('rel-t-tol', 1.0e-6_r8)
    call plist%set('rel-h-tol', 1.0e-6_r8)
    plist => params%sublist('integrator')
    call plist%set('nlk-max-iter', 5)
    call plist%set('nlk-tol', 0.01_r8)
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine


  subroutine fail(message)
    character(*), intent(in) :: message

    if (is_IOP) print '("ERROR: ",a)', message
    error stop 1
  end subroutine

end program test_ns_ht_2d_solver
