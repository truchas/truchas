program test_flow_2d_projection

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use parameter_list_type
  use simulation_environment_type
  use flow_2d_operators_type
  use flow_2d_projection_type
  use flow_2d_bc_type
  use pcsr_matrix_type
  implicit none

  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_flow_2d_projection.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_projection

  call halt_parallel_communication
  stop status

contains

  subroutine test_projection
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_operators), target :: operators
    type(flow_2d_projection), target :: projection
    type(pcsr_matrix), pointer :: matrix
    type(parameter_list), target :: neumann_params, dirichlet_params
    type(parameter_list), pointer :: plist
    type(flow_2d_bc) :: bc
    real(r8), allocatable :: inv_density_f(:), p(:), one(:), result(:), rhs(:)
    real(r8), parameter :: gradient(2) = [1.0_r8, 2.0_r8]
    character(:), allocatable :: errmsg
    integer :: c, stat
    logical :: linear_harmonic, pinned, dirichlet_solution

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call rotate_mesh(mesh, 30.0_r8)
    call operators%init(mesh)
    call projection%init(mesh, operators)
    allocate(inv_density_f(mesh%nface), p(mesh%ncell), one(mesh%ncell), result(mesh%ncell_onP), &
        rhs(mesh%ncell_onP))
    inv_density_f = 1.0_r8
    one = 1.0_r8
    do c = 1, mesh%ncell
      p(c) = dot_product(gradient, mesh%cell_centroid(:,c))
    end do

    plist => neumann_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1])
    call bc%init(env, mesh, neumann_params, stat, errmsg)
    call require(stat == 0, 'Neumann boundary condition initialization failed')
    call bc%compute(0.0_r8)
    call projection%assemble(inv_density_f, bc, rhs)
    matrix => projection%matrix()
    call matrix%matvec(one, result)
    pinned = maxval(abs(result)) > 1.0e-12_r8
    call require(global_any(pinned), 'all-Neumann pressure conditions did not remove the constant nullspace')
    call require(maxval(abs(rhs)) < 1.0e-12_r8, 'homogeneous Neumann condition changed projection RHS')

    call matrix%matvec(p, result)
    linear_harmonic = .true.
    do c = 1, mesh%ncell_onP
      if (all(mesh%cnhbr(mesh%cstart(c):mesh%cstart(c+1)-1) > 0)) then
        linear_harmonic = linear_harmonic .and. abs(result(c)) < 1.0e-12_r8
      end if
    end do
    call require(linear_harmonic, 'projection operator is not rotation invariant on an orthogonal mesh')

    plist => dirichlet_params%sublist('outlet')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    call plist%set('pressure', 3.0_r8)
    call bc%init(env, mesh, dirichlet_params, stat, errmsg)
    call require(stat == 0, 'Dirichlet boundary condition initialization failed')
    call bc%compute(0.0_r8)
    call projection%assemble(inv_density_f, bc, rhs)
    p = 3.0_r8
    call matrix%matvec(p, result)
    dirichlet_solution = maxval(abs(result - rhs)) < 1.0e-12_r8
    call require(dirichlet_solution, 'constant pressure Dirichlet solution does not satisfy projection system')
    call require(global_any(abs(rhs) > 1.0e-12_r8), 'pressure Dirichlet condition did not contribute to projection RHS')
  end subroutine


  subroutine rotate_mesh(mesh, angle)
    type(unstr_2d_mesh), intent(inout) :: mesh
    real(r8), intent(in) :: angle

    real(r8) :: theta, rotation(2,2)

    theta = angle*acos(-1.0_r8)/180.0_r8
    rotation = reshape([cos(theta), sin(theta), -sin(theta), cos(theta)], [2,2])
    mesh%x = matmul(rotation, mesh%x)
    call mesh%compute_geometry()
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_projection
