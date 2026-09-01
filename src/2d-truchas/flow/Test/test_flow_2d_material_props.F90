program test_flow_2d_material_props

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_domain_types
  use flow_2d_material_props_type
  implicit none

  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env
  type(unstr_2d_mesh), pointer :: mesh
  type(flow_2d_material_props) :: props
  real(r8) :: vfrac(1,2), vfrac_void(2,2)
  integer :: interior_face

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_flow_2d_material_props.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1], 0.0_r8, 0.0_r8)
  interior_face = 0
  do stat = mesh%cstart(1), mesh%cstart(2)-1
    if (mesh%cnhbr(stat) == 2) then
      interior_face = mesh%cface(stat)
      exit
    end if
  end do
  call require(interior_face > 0, 'could not find interior face')

  status = 0
  call props%init(mesh, [2.0_r8], .true., stat, errmsg)
  call require(stat == 0, 'initializing fluid/solid properties: ' // errmsg)
  vfrac = reshape([1.0_r8, 0.0_r8], shape(vfrac))
  call props%set_volume_fractions(vfrac)
  call require(all(props%cell_t == [regular_t, solid_t]), 'wrong fluid/solid cell classification')
  call require(props%face_t(interior_face) == solid_t, 'wrong fluid/solid face classification')
  call require(all(abs(props%vof - [1.0_r8, 0.0_r8]) < 1.0e-14_r8), 'wrong fluid/solid fluid fractions')
  call require(maxval(abs(props%inv_density_c - [0.5_r8, 0.0_r8])) < 1.0e-14_r8, &
      'wrong fluid/solid inverse cell density')
  call require(abs(props%inv_density_f(interior_face) - 0.5_r8) < 1.0e-14_r8, &
      'wrong fluid/solid inverse face density')

  vfrac = reshape([1.0_r8, 0.5_r8], shape(vfrac))
  call props%set_volume_fractions(vfrac)
  call require(all(props%cell_t == [regular_t, regular_t]), 'wrong partially solid cell classification')
  call require(props%face_t(interior_face) == regular_t, 'wrong partially solid face classification')
  call require(all(abs(props%vof - [1.0_r8, 0.5_r8]) < 1.0e-14_r8), 'wrong partially solid fluid fractions')
  call require(maxval(abs(props%density_c - [2.0_r8, 2.0_r8])) < 1.0e-14_r8, &
      'wrong partially solid cell density')

  call props%init(mesh, [2.0_r8], .true., stat, errmsg, nfluid=2)
  call require(stat == 0, 'initializing fluid/VOID properties: ' // errmsg)
  vfrac_void = reshape([1.0_r8, 0.0_r8, 0.0_r8, 1.0_r8], shape(vfrac_void))
  call props%set_volume_fractions(vfrac_void)
  call require(all(props%cell_t == [regular_t, void_t]), 'wrong fluid/VOID cell classification')
  call require(props%face_t(interior_face) == regular_void_t, 'wrong fluid/VOID face classification')
  call require(props%any_void, 'fluid/VOID case did not report VOID')

  call env%simlog%close()
  call halt_parallel_communication
  stop status

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not.condition) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine require

end program test_flow_2d_material_props
