program test_flow_2d_solver

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
  use flow_2d_solver_type
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
  call env%simlog%init(env%comm, 'test_flow_2d_solver.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_step
  call test_pressure_drive
  call test_incompatible_flux

  call halt_parallel_communication
  stop status

contains

  subroutine test_step
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_model), target :: model
    type(flow_2d_solver), target :: solver
    type(parameter_list), target :: bc_params, momentum_params, projection_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: flux(:)
    real(r8), pointer :: pressure(:), velocity_state(:,:), velocity_face(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1,2,3,4])
    call model%init(env, mesh, bc_params, [1.0_r8], 1.0_r8, stat, errmsg)
    call require(stat == 0, 'flow model initialization failed')
    if (stat /= 0) return
    call momentum_params%set('rel-tol', 1.0e-10_r8)
    call momentum_params%set('max-ds-iter', 100)
    call momentum_params%set('max-amg-iter', 100)
    call projection_params%set('rel-tol', 1.0e-10_r8)
    call projection_params%set('max-ds-iter', 100)
    call projection_params%set('max-amg-iter', 100)
    call solver%init(model, momentum_params, projection_params)

    call solver%get_cell_flow_soln(pressure, velocity_state)
    velocity_state = spread([1.0_r8, -0.5_r8], dim=2, ncopies=mesh%ncell)
    call solver%step(0.0_r8, 1.0_r8, stat)
    call require(stat == 0, 'flow solver step did not converge')
    allocate(flux(mesh%ncell_onP))
    call solver%get_face_velocity(velocity_face)
    call model%operators%divergence(velocity_face, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, 'flow solver step did not make face velocity solenoidal')
  end subroutine


  subroutine test_incompatible_flux
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_model), target :: model
    type(flow_2d_solver), target :: solver
    type(parameter_list), target :: bc_params, momentum_params, projection_params
    type(parameter_list), pointer :: plist
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('inlet')
    call plist%set('type', 'velocity')
    call plist%set('face-set-ids', [1])
    call plist%set('velocity', [1.0_r8, 0.0_r8])
    plist => bc_params%sublist('walls')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [2,3,4])
    call model%init(env, mesh, bc_params, [1.0_r8], 1.0_r8, stat, errmsg)
    call require(stat == 0, 'incompatible-flux model initialization failed')
    if (stat /= 0) return
    call momentum_params%set('rel-tol', 1.0e-10_r8)
    call momentum_params%set('max-ds-iter', 100)
    call momentum_params%set('max-amg-iter', 100)
    call projection_params%set('rel-tol', 1.0e-10_r8)
    call projection_params%set('max-ds-iter', 100)
    call projection_params%set('max-amg-iter', 100)
    call solver%init(model, momentum_params, projection_params)
    call solver%step(0.0_r8, 1.0_r8, stat, errmsg)
    call require(stat /= 0, 'incompatible prescribed flux was not rejected')
  end subroutine


  !! A pressure boundary leaves the normal velocity unconstrained.  This is
  !! the minimal pressure-driven channel case: pressure decreases from the
  !! left face set to the right one, while the horizontal walls are no-slip.
  subroutine test_pressure_drive
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_model), target :: model
    type(flow_2d_solver), target :: solver
    type(parameter_list), target :: bc_params, momentum_params, projection_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: flux(:), expected_velocity(:)
    real(r8), pointer :: pressure(:), velocity_state(:,:), velocity_face(:)
    real(r8) :: time, mean_velocity
    character(:), allocatable :: errmsg
    integer :: stat, n

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [16, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('inlet')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    call plist%set('pressure', 1.0_r8)
    plist => bc_params%sublist('outlet')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [2])
    call plist%set('pressure', 0.0_r8)
    plist => bc_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [3,4])
    call model%init(env, mesh, bc_params, [1.0_r8], 1.0_r8, stat, errmsg)
    call require(stat == 0, 'pressure-driven flow model initialization failed')
    if (stat /= 0) return
    call momentum_params%set('rel-tol', 1.0e-10_r8)
    call momentum_params%set('max-ds-iter', 100)
    call momentum_params%set('max-amg-iter', 100)
    call projection_params%set('rel-tol', 1.0e-10_r8)
    call projection_params%set('max-ds-iter', 100)
    call projection_params%set('max-amg-iter', 100)
    call solver%init(model, momentum_params, projection_params)

    time = 0.0_r8
    do n = 1, 50
      call solver%step(time, real(n, r8), stat)
      call require(stat == 0, 'pressure-driven flow step did not converge')
      if (stat /= 0) return
      time = real(n, r8)
    end do
    allocate(flux(mesh%ncell_onP))
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    call model%operators%divergence(velocity_face, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, 'pressure-driven flow is not solenoidal')
    mean_velocity = global_sum(sum(velocity_state(1,1:mesh%ncell_onP))) / &
        global_sum(real(mesh%ncell_onP, r8))
    call require(mean_velocity > 1.0e-3_r8, 'pressure gradient did not drive flow from inlet to outlet')
    call mesh%init_cell_centroid
    allocate(expected_velocity(mesh%ncell_onP))
    expected_velocity = 0.5_r8 * mesh%cell_centroid(2,1:mesh%ncell_onP) * &
        (1.0_r8 - mesh%cell_centroid(2,1:mesh%ncell_onP))
    call require(maxval(abs(velocity_state(1,1:mesh%ncell_onP) - expected_velocity)) < 5.0e-3_r8, &
        'long-time pressure-driven flow does not match the Poiseuille profile')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_solver
