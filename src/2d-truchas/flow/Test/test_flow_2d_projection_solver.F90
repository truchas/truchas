program test_flow_2d_projection_solver

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
  use flow_2d_projection_type
  use flow_2d_projection_solver_type
  use flow_2d_bc_type
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
  call env%simlog%init(env%comm, 'test_flow_2d_projection_solver.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_solver

  call halt_parallel_communication
  stop status

contains

  subroutine test_solver
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_operators), target :: operators
    type(flow_2d_projection), target :: projection
    type(flow_2d_projection_solver) :: solver
    type(flow_2d_bc) :: bc
    type(parameter_list), target :: bc_params, solver_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: inv_density_f(:), rhs(:), pressure(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call operators%init(mesh)
    call projection%init(mesh, operators)
    allocate(inv_density_f(mesh%nface), rhs(mesh%ncell_onP), pressure(mesh%ncell))
    inv_density_f = 1.0_r8

    plist => bc_params%sublist('outlet')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    call plist%set('pressure', 3.0_r8)
    call bc%init(env, mesh, bc_params, stat, errmsg)
    call require(stat == 0, 'pressure boundary condition initialization failed')
    call bc%compute(0.0_r8)
    call projection%assemble(inv_density_f, bc, rhs)

    call solver_params%set('rel-tol', 1.0e-10_r8)
    call solver_params%set('max-ds-iter', 100)
    call solver_params%set('max-amg-iter', 100)
    call solver%init(projection, solver_params)
    call solver%setup()
    pressure = 0.0_r8
    call solver%solve(rhs, pressure, stat)
    call require(stat == 0, 'projection solve did not converge')
    call require(maxval(abs(pressure(1:mesh%ncell_onP) - 3.0_r8)) < 1.0e-8_r8, &
        'projection solution is incorrect')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_projection_solver
