program test_ns_2d_solver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_model_type
  use ns_2d_solver_type
  implicit none

  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  call init_parallel_communication
  call fhypre_initialize
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_ns_2d_solver.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_step

  call halt_parallel_communication
  stop status

contains

  subroutine test_step
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_model), target :: model
    type(ns_2d_solver), target :: solver
    type(parameter_list), target :: bc_params, momentum_params, projection_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: velocity(:,:), flux(:), flux_volumes(:,:), pressure_save(:), velocity_save(:,:), &
        face_velocity_save(:), pressure_trial(:), velocity_trial(:,:), face_velocity_trial(:)
    real(r8), pointer :: pressure(:), velocity_state(:,:), velocity_face(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1,2,3,4])
    call model%init(env, mesh, bc_params, [1.0_r8], 0.1_r8, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes model initialization failed')
    if (stat /= 0) return
    call set_solver_params(momentum_params)
    call set_solver_params(projection_params)
    call solver%init(env, model, momentum_params, projection_params)
    allocate(velocity(2,mesh%ncell_onP), flux(mesh%ncell_onP), flux_volumes(1,size(mesh%cface)))
    flux_volumes = 0.0_r8
    velocity = spread([1.0_r8, -0.5_r8], dim=2, ncopies=mesh%ncell_onP)
    call solver%set_initial_state(env, 0.0_r8, 0.01_r8, velocity, stat)
    call require(stat == 0, 'Navier--Stokes initial-condition solve did not converge')
    if (stat /= 0) return
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    allocate(pressure_save(size(pressure)), velocity_save(size(velocity_state,1),size(velocity_state,2)), &
        face_velocity_save(size(velocity_face)), pressure_trial(size(pressure)), &
        velocity_trial(size(velocity_state,1),size(velocity_state,2)), face_velocity_trial(size(velocity_face)))
    pressure_save = pressure
    velocity_save = velocity_state
    face_velocity_save = velocity_face
    call solver%advance_momentum(env, 0.0_r8, 0.01_r8, flux_volumes, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes momentum update did not converge')
    if (stat /= 0) return
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    pressure_trial = pressure
    velocity_trial = velocity_state
    face_velocity_trial = velocity_face
    call solver%reject_step()
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    call require(maxval(abs(pressure - pressure_save)) == 0.0_r8, &
        'flow pressure was not restored by reject_step')
    call require(maxval(abs(velocity_state - velocity_save)) == 0.0_r8, &
        'cell velocity was not restored by reject_step')
    call require(maxval(abs(velocity_face - face_velocity_save)) == 0.0_r8, &
        'face velocity was not restored by reject_step')
    call solver%advance_momentum(env, 0.0_r8, 0.01_r8, flux_volumes, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes second momentum update did not converge')
    if (stat /= 0) return
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    pressure_trial = pressure
    velocity_trial = velocity_state
    face_velocity_trial = velocity_face
    call solver%commit_step()
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    call require(maxval(abs(pressure - pressure_trial)) == 0.0_r8, &
        'flow pressure changed when committing the pending state')
    call require(maxval(abs(velocity_state - velocity_trial)) == 0.0_r8, &
        'cell velocity changed when committing the pending state')
    call require(maxval(abs(velocity_face - face_velocity_trial)) == 0.0_r8, &
        'face velocity changed when committing the pending state')
    call solver%step(env, 0.0_r8, 0.01_r8, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes solver step did not converge')
    if (stat /= 0) return
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    call model%operators%divergence(velocity_face, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, 'Navier--Stokes step did not make face velocity solenoidal')
    call require(solver%courant_time_step() > 0.0_r8, 'Navier--Stokes Courant time step is not positive')
  end subroutine


  subroutine set_solver_params(params)
    type(parameter_list), intent(inout) :: params

    call params%set('rel-tol', 1.0e-10_r8)
    call params%set('max-ds-iter', 100)
    call params%set('max-amg-iter', 100)
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_ns_2d_solver
