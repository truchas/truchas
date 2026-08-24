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
  call env%simlog%init(env%comm, 'test_flow_2d_material_transport.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_single_liquid_flux_volumes

  call halt_parallel_communication
  stop status

contains

  subroutine test_single_liquid_flux_volumes
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_material_transport) :: transport
    real(r8), allocatable :: velocity_fn(:)
    real(r8) :: expected
    integer :: c, f, i
    logical :: passed

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call mesh%init_face_centroid
    call transport%init(mesh, 1)
    allocate(velocity_fn(mesh%nface), source=0.0_r8)
    do f = 1, mesh%nface_onP
      velocity_fn(f) = 2.0_r8 - 0.25_r8*mesh%face_centroid(1,f) + 0.5_r8*mesh%face_centroid(2,f)
    end do
    call mesh%face_imap%gather_offp(velocity_fn)

    call transport%advance(1.0_r8, 1.125_r8, velocity_fn)
    passed = .true.
    do c = 1, mesh%ncell_onP
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        f = mesh%cface(i)
        expected = 0.125_r8*mesh%area(f)*velocity_fn(f)
        if (btest(mesh%cfpar(c), i-mesh%cstart(c)+1)) expected = -expected
        passed = passed .and. abs(transport%flux_volumes(1,i) - expected) < 1.0e-14_r8
      end do
    end do
    call require(passed, 'single-liquid flux volume has incorrect magnitude or orientation')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_material_transport
