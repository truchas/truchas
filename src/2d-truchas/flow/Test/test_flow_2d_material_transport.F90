program test_flow_2d_material_transport

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_material_transport_type
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
  allocate(env%timer)
  call env%simlog%init(env%comm, 'test_flow_2d_material_transport.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_single_liquid_flux_volumes('simple')
  call test_single_liquid_flux_volumes('geometric')
  if (env%nproc == 1) call test_fluid_solid_transport

  deallocate(env%timer)
  call halt_parallel_communication
  stop status

contains

  subroutine test_single_liquid_flux_volumes(algorithm)
    character(*), intent(in) :: algorithm

    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_material_transport) :: transport
    real(r8), allocatable :: velocity_fn(:), expected_flux(:), vfrac_n(:,:)
    real(r8), pointer :: vfrac_trial(:,:)
    real(r8) :: dt, expected
    integer :: c, f, i, ncf_onp
    logical :: passed

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call mesh%init_face_centroid
    call transport%init(env, mesh, 1, 1, 1, algorithm)
    dt = 0.125_r8
    if (algorithm == 'geometric') dt = 0.01_r8
    allocate(velocity_fn(mesh%nface), expected_flux(size(mesh%cface)), vfrac_n(1,mesh%ncell), source=0.0_r8)
    vfrac_n = 1.0_r8
    expected_flux = 0.0_r8
    do f = 1, mesh%nface_onP
      velocity_fn(f) = dot_product([0.2_r8, -0.1_r8], mesh%unit_normal(:,f))
    end do
    call mesh%face_imap%gather_offp(velocity_fn)
    call transport%advance(1.0_r8, 1.0_r8 + dt, velocity_fn, vfrac_n)
    call transport%get_trial_volume_fractions(vfrac_trial)
    ncf_onp = mesh%cstart(mesh%ncell_onP+1)-1
    do c = 1, mesh%ncell_onP
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        f = mesh%cface(i)
        expected = dt*mesh%area(f)*velocity_fn(f)
        if (btest(mesh%cfpar(c), i-mesh%cstart(c)+1)) expected = -expected
        expected_flux(i) = expected
      end do
    end do
    passed = maxval(abs(transport%flux_volumes(1,:ncf_onp) - expected_flux(:ncf_onp))) < 1.0e-14_r8
    call require(passed, trim(algorithm) // ': single-liquid flux volume has incorrect magnitude or orientation')
    call require(maxval(abs(vfrac_trial-1.0_r8)) < 1.0e-14_r8, &
        trim(algorithm) // ': single-liquid volume fraction changed')
  end subroutine


  subroutine test_fluid_solid_transport

    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_material_transport) :: transport
    real(r8), allocatable :: velocity_fn(:), vfrac_n(:,:)
    real(r8), pointer :: vfrac_trial(:,:)
    real(r8) :: q
    integer :: f, j1, j2

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1], 0.0_r8, 0.0_r8)
    call transport%init(env, mesh, 1, 1, 2, 'simple')
    allocate(velocity_fn(mesh%nface), vfrac_n(2,mesh%ncell), source=0.0_r8)
    vfrac_n = 0.5_r8

    f = 0
    do j1 = mesh%cstart(1), mesh%cstart(2)-1
      if (mesh%cnhbr(j1) == 2) then
        f = mesh%cface(j1)
        exit
      end if
    end do
    call require(f > 0, 'could not find interior face')
    j2 = mesh%cstart(2) - 1 + &
        findloc(mesh%cface(mesh%cstart(2):mesh%cstart(3)-1), f, dim=1)
    velocity_fn(f) = 0.1_r8
    if (btest(mesh%cfpar(1), j1-mesh%cstart(1)+1)) velocity_fn(f) = -velocity_fn(f)

    call transport%advance(0.0_r8, 0.1_r8, velocity_fn, vfrac_n)
    call transport%get_trial_volume_fractions(vfrac_trial)
    q = 0.5_r8 * 0.1_r8 * 0.1_r8 * mesh%area(f) / mesh%volume(1)
    call require(maxval(abs(vfrac_trial(1,:2) - [0.5_r8-q, 0.5_r8+q])) < 1.0e-14_r8, &
        'fluid slot was not transported')
    call require(maxval(abs(vfrac_trial(2,:2) - 0.5_r8)) < 1.0e-14_r8, &
        'solid slot was transported')
    call require(size(transport%flux_volumes,1) == 1, 'solid slot is present in flux volumes')
    call require(abs(transport%flux_volumes(1,j1)-q*mesh%volume(1)) < 1.0e-14_r8, &
        'fluid flux volume is incorrect')
    call require(abs(transport%flux_volumes(1,j2)+q*mesh%volume(1)) < 1.0e-14_r8, &
        'fluid flux volume orientation is incorrect')

  end subroutine test_fluid_solid_transport


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_material_transport
