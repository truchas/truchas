program test_flow_2d_momentum_solver

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
  use flow_2d_operators_type
  use flow_2d_bc_type
  use flow_2d_momentum_type
  use flow_2d_momentum_solver_type
  use flow_domain_types
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
  call env%simlog%init(env%comm, 'test_flow_2d_momentum_solver.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_solver

  call halt_parallel_communication
  stop status

contains

  subroutine test_solver
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_operators), target :: operators
    type(flow_2d_momentum), target :: momentum
    type(flow_2d_momentum_solver) :: solver
    type(flow_2d_bc) :: bc
    type(parameter_list), target :: bc_params, solver_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: density(:), viscosity(:), rhs(:,:), velocity(:,:)
    integer, allocatable :: cell_t(:), face_t(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call operators%init(mesh)
    call momentum%init(mesh, operators)
    allocate(density(mesh%ncell), viscosity(mesh%nface), rhs(2,mesh%ncell_onP), velocity(2,mesh%ncell), &
        cell_t(mesh%ncell), face_t(mesh%nface))
    density = 0.0_r8
    viscosity = 1.0_r8
    cell_t = regular_t
    face_t = regular_t

    plist => bc_params%sublist('wall')
    call plist%set('type', 'velocity')
    call plist%set('face-set-ids', [1,2,3,4])
    call plist%set('velocity', [1.5_r8, -0.75_r8])
    call bc%init(env, mesh, bc_params, stat, errmsg)
    call require(stat == 0, 'velocity boundary condition initialization failed')
    call bc%compute(0.0_r8)
    call momentum%assemble(1.0_r8, density, viscosity, cell_t, face_t, bc, rhs)

    call solver_params%set('rel-tol', 1.0e-10_r8)
    call solver_params%set('max-ds-iter', 100)
    call solver_params%set('max-amg-iter', 100)
    call solver%init(momentum, solver_params)
    call solver%setup()
    velocity = 0.0_r8
    call solver%solve(rhs, velocity(:,1:mesh%ncell_onP), stat)
    call require(stat == 0, 'momentum solve did not converge')
    call require(maxval(abs(velocity(:,1:mesh%ncell_onP) - spread([1.5_r8, -0.75_r8], &
        dim=2, ncopies=mesh%ncell_onP))) < 1.0e-8_r8, 'momentum solution is incorrect')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_momentum_solver
