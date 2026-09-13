program test_t2d_flow_projection_update

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use simulation_environment_type
  use t2d_unstr_mesh_type
  use t2d_unstr_mesh_factory
  use t2d_flow_state_type
  use t2d_flow_operators_type
  use t2d_flow_bc_type
  use t2d_flow_projection_type
  use t2d_flow_projection_solver_type
  use t2d_flow_projection_update_type
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
  call env%simlog%init(env%comm, 'test_t2d_flow_projection_update.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_update
  call test_compliance

  call halt_parallel_communication
  stop status

contains

  subroutine test_update
    type(t2d_unstr_mesh), pointer :: mesh
    type(t2d_flow_operators), target :: operators
    type(t2d_flow_projection), target :: projection
    type(t2d_flow_projection_solver), target :: solver
    type(t2d_flow_projection_update) :: update
    type(t2d_flow_state) :: state
    type(t2d_flow_bc) :: bc
    type(parameter_list), target :: bc_params, solver_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: inv_density_c(:), inv_density_f(:), density_delta_c(:), flux(:)
    integer, allocatable :: cell_t(:), face_t(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call operators%init(mesh)
    call projection%init(mesh, operators)
    call state%init(mesh)
    allocate(inv_density_c(mesh%ncell), inv_density_f(mesh%nface), density_delta_c(mesh%ncell), &
        flux(mesh%ncell_onP), cell_t(mesh%ncell), face_t(mesh%nface))
    inv_density_c = 1.0_r8
    inv_density_f = 1.0_r8
    density_delta_c = 0.0_r8
    cell_t = regular_t
    face_t = regular_t

    plist => bc_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1,2,3,4])
    call bc%init(env, mesh, bc_params, stat, errmsg)
    call require(stat == 0, 'no-slip boundary condition initialization failed')
    call bc%compute(0.0_r8, 1.0_r8)

    call solver_params%set('rel-tol', 1.0e-10_r8)
    call solver_params%set('max-ds-iter', 100)
    call solver_params%set('max-amg-iter', 100)
    call solver%init(projection, solver_params)
    call update%init(mesh, operators, projection, solver)

    state%vel_cc = spread([1.0_r8, -0.5_r8], dim=2, ncopies=mesh%ncell)
    call update%correct(1.0_r8, inv_density_c, inv_density_f, density_delta_c, cell_t, face_t, bc, state, stat)
    call require(stat == 0, 'projection update did not converge')
    call operators%divergence(state%vel_fn, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, 'projection update did not make face velocity solenoidal')
  end subroutine


  subroutine test_compliance
    type(t2d_unstr_mesh), pointer :: mesh
    type(t2d_flow_operators), target :: operators
    type(t2d_flow_projection), target :: projection
    type(t2d_flow_projection_solver), target :: solver
    type(t2d_flow_projection_update) :: update
    type(t2d_flow_state) :: state
    type(t2d_flow_bc) :: bc
    type(parameter_list), target :: bc_params, solver_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: inv_c(:), inv_f(:), delta(:), flux(:), compliance(:)
    integer, allocatable :: cell_t(:), face_t(:)
    integer :: sign_p, stat
    character(:), allocatable :: errmsg
    real(r8), parameter :: dt = 0.2_r8

    mesh => new_unstr_2d_mesh(env, [0.0_r8,0.0_r8], [1.0_r8,1.0_r8], [8,8], 0.0_r8, 0.0_r8)
    call operators%init(mesh)
    call projection%init(mesh, operators)
    call state%init(mesh)
    allocate(inv_c(mesh%ncell), inv_f(mesh%nface), delta(mesh%ncell), flux(mesh%ncell_onP), &
        compliance(mesh%ncell_onP), cell_t(mesh%ncell), face_t(mesh%nface))
    inv_c = 1.0_r8
    inv_f = 1.0_r8
    delta = 0.0_r8
    compliance = 0.5_r8
    cell_t = regular_void_t
    face_t = regular_t
    call solver_params%set('rel-tol', 1.0e-12_r8)
    call solver_params%set('max-ds-iter', 100)
    call solver_params%set('max-amg-iter', 100)
    call solver%init(projection, solver_params)
    call update%init(mesh, operators, projection, solver)
    plist => bc_params%sublist('walls')
    call plist%set('type', 'free-slip')
    call plist%set('face-set-ids', [2,3,4])
    plist => bc_params%sublist('loading')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    do sign_p = -1, 1, 2
      call plist%set('pressure', real(sign_p,r8))
      call bc%init(env, mesh, bc_params, stat, errmsg)
      call require(stat == 0, 'compliance boundary initialization failed')
      call bc%compute(0.0_r8, dt)
      state%vel_cc = 0.0_r8
      state%p_cc = 3.0_r8 ! Detect a missing old-pressure offset.
      call update%correct(dt, inv_c, inv_f, delta, cell_t, face_t, bc, state, stat, &
          env=env, compliance=compliance, initial=.true.)
      call require(stat == 0, 'compliance projection did not converge')
      call operators%divergence(state%vel_fn, flux)
      call require(maxval(abs(flux/mesh%volume(:mesh%ncell_onP) + &
          compliance*state%p_cc(:mesh%ncell_onP))) < 1.0e-8_r8, 'compliance divergence identity failed')
      call require(global_maxval(sign_p*flux) < 0.0_r8, 'two-sided compliance has incorrect divergence sign')
    end do
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_t2d_flow_projection_update
