program test_flow_2d_momentum

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_operators_type
  use flow_2d_bc_type
  use flow_2d_momentum_type
  use pbsr_matrix_type
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
  call env%simlog%init(env%comm, 'test_flow_2d_momentum.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_momentum
  call test_advective_transport

  call halt_parallel_communication
  stop status

contains

  subroutine test_momentum
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_operators), target :: operators
    type(flow_2d_momentum), target :: momentum
    type(flow_2d_bc) :: bc
    type(parameter_list), target :: velocity_params, slip_params
    type(parameter_list), pointer :: plist
    type(pbsr_matrix), pointer :: matrix
    real(r8), allocatable :: density(:), viscosity(:), velocity(:,:), rhs(:,:), result(:,:)
    character(:), allocatable :: errmsg
    integer :: stat, f, c, entry
    logical :: dirichlet_solution, coupled_slip

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call rotate_mesh(mesh, 30.0_r8)
    call operators%init(mesh)
    call momentum%init(mesh, operators)
    allocate(density(mesh%ncell), viscosity(mesh%nface), velocity(2,mesh%ncell), &
        rhs(2,mesh%ncell_onP), result(2,mesh%ncell_onP))
    density = 0.0_r8
    viscosity = 1.0_r8
    velocity = spread([1.5_r8, -0.75_r8], dim=2, ncopies=mesh%ncell)

    plist => velocity_params%sublist('wall')
    call plist%set('type', 'velocity')
    call plist%set('face-set-ids', [1,2,3,4])
    call plist%set('velocity', [1.5_r8, -0.75_r8])
    call bc%init(mesh, velocity_params, stat, errmsg)
    call require(stat == 0, 'velocity boundary condition initialization failed')
    call bc%compute(0.0_r8)
    call momentum%assemble(1.0_r8, density, viscosity, bc, rhs)
    matrix => momentum%matrix()
    call matrix%matvec(velocity, result)
    dirichlet_solution = maxval(abs(result - rhs)) < 1.0e-12_r8
    call require(dirichlet_solution, 'constant velocity Dirichlet solution does not satisfy momentum system')

    plist => slip_params%sublist('wall')
    call plist%set('type', 'free-slip')
    call plist%set('face-set-ids', [1])
    call bc%init(mesh, slip_params, stat, errmsg)
    call require(stat == 0, 'free-slip boundary condition initialization failed')
    call bc%compute(0.0_r8)
    call momentum%assemble(1.0_r8, density, viscosity, bc, rhs)
    matrix => momentum%matrix()
    coupled_slip = .false.
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle
      if (abs(mesh%unit_normal(1,f)*mesh%unit_normal(2,f)) < 1.0e-12_r8) cycle
      c = mesh%fcell(1,f)
      entry = matrix%graph%index(c, c)
      coupled_slip = abs(matrix%values(1,2,entry)) > 1.0e-12_r8 .and. &
          abs(matrix%values(1,2,entry) - matrix%values(2,1,entry)) < 1.0e-12_r8
      if (coupled_slip) exit
    end do
    call require(global_any(coupled_slip), 'free-slip condition did not add a symmetric velocity block coupling')
  end subroutine


  subroutine test_advective_transport
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_operators), target :: operators
    type(flow_2d_momentum), target :: momentum
    type(flow_2d_bc) :: bc
    type(parameter_list), target :: velocity_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: density(:), velocity_cc(:,:), velocity_fn(:), flux_volumes(:,:), rhs(:,:)
    real(r8) :: transported_velocity(2), material_fraction(2), result(2), expected(2)
    character(:), allocatable :: errmsg
    integer :: stat, f, n

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call operators%init(mesh)
    call momentum%init(mesh, operators)
    allocate(density(2), velocity_cc(2,mesh%ncell), velocity_fn(mesh%nface), &
        flux_volumes(2,size(mesh%cface)), rhs(2,mesh%ncell_onP))
    density = [1.0_r8, 3.0_r8]
    material_fraction = [0.25_r8, 0.75_r8]
    transported_velocity = [1.5_r8, -0.75_r8]
    velocity_cc = spread(transported_velocity, dim=2, ncopies=mesh%ncell)
    velocity_fn = 0.0_r8
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) > 0) velocity_fn(f) = dot_product(transported_velocity, mesh%unit_normal(:,f))
    end do
    call mesh%face_imap%gather_offp(velocity_fn)
    call compute_flux_volumes(mesh, 0.25_r8, velocity_fn, material_fraction, flux_volumes)
    call bc%init(mesh, velocity_params, stat, errmsg)
    call require(stat == 0, 'advective-transport boundary condition initialization failed')
    if (stat /= 0) return
    call bc%compute(0.0_r8)
    rhs = 0.0_r8
    call momentum%add_advective_rhs(density, velocity_cc, flux_volumes, bc, rhs)
    result = [global_sum(sum(rhs(1,:))), global_sum(sum(rhs(2,:)))]
    call require(maxval(abs(result)) < 1.0e-12_r8, 'interior advective transport is not conservative')

    plist => velocity_params%sublist('inlet')
    call plist%set('type', 'velocity')
    call plist%set('face-set-ids', [1])
    call plist%set('velocity', transported_velocity)
    call bc%init(mesh, velocity_params, stat, errmsg)
    call require(stat == 0, 'advective-inflow boundary condition initialization failed')
    if (stat /= 0) return
    call bc%compute(0.0_r8)
    velocity_fn = 0.0_r8
    expected = 0.0_r8
    do n = 1, size(bc%velocity_dirichlet%index)
      f = bc%velocity_dirichlet%index(n)
      if (f > mesh%nface_onP) cycle
      velocity_fn(f) = dot_product(mesh%unit_normal(:,f), bc%velocity_dirichlet%value(:,n))
      expected = expected - 0.25_r8*dot_product(density, material_fraction)*mesh%area(f)*velocity_fn(f) * &
          bc%velocity_dirichlet%value(:,n)
    end do
    call mesh%face_imap%gather_offp(velocity_fn)
    call compute_flux_volumes(mesh, 0.25_r8, velocity_fn, material_fraction, flux_volumes)
    rhs = 0.0_r8
    call momentum%add_advective_rhs(density, velocity_cc, flux_volumes, bc, rhs)
    result = [global_sum(sum(rhs(1,:))), global_sum(sum(rhs(2,:)))]
    expected = [global_sum(expected(1)), global_sum(expected(2))]
    call require(maxval(abs(result - expected)) < 1.0e-12_r8, &
        'velocity-boundary inflow does not supply donor momentum')
  end subroutine


  subroutine compute_flux_volumes(mesh, dt, velocity_fn, material_fraction, flux_volumes)
    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: dt, velocity_fn(:), material_fraction(:)
    real(r8), intent(out) :: flux_volumes(:,:)

    integer :: c, i, f

    flux_volumes = 0.0_r8
    do c = 1, mesh%ncell_onP
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        f = mesh%cface(i)
        flux_volumes(:,i) = material_fraction*dt*mesh%area(f)*velocity_fn(f)
        if (btest(mesh%cfpar(c), i-mesh%cstart(c)+1)) flux_volumes(:,i) = -flux_volumes(:,i)
      end do
    end do
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

end program test_flow_2d_momentum
