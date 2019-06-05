program test_pressure_poisson
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use pgslib_module
  use truchas_env, only: prefix
  use truchas_logging_services
  use parameter_list_type
  use flow_operators
  use unstr_mesh_type
  use scalar_func_class
  use scalar_func_tools
  use bndry_func1_class
  implicit none

  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()
  type(parameter_list), pointer :: p, pp
  !type(flow_bc), pointer :: bc
  type(unstr_mesh), pointer :: mesh
  integer :: i, in

  call parallel_init(argv)
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  call TLS_initialize()
  call TLS_set_verbosity(TLS_VERB_NOISY)
  call init_mesh(mesh)
  call flow_operators_init(mesh)

  call interpolate_fc_test(mesh)
  call halt_parallel_communication()

contains

  subroutine init_mesh(mesh)

    use mesh_manager
    use unstr_mesh_type

    type(unstr_mesh), pointer, intent(out) :: mesh
    character(:), allocatable :: indir
    type(parameter_list) :: pmesh
    type(parameter_list), pointer :: p

    call process_command_line(indir)

    p => pmesh%sublist("mesh")
    call p%set('mesh-file', indir//"/square-32x32.g")
    call init_mesh_manager(pmesh)

    mesh => unstr_mesh_ptr("MESH")
    call mesh%init_cell_centroid
    call mesh%init_face_centroid
    call mesh%init_face_normal_dist
    print *, "INTIALIZED MESH"
  end subroutine init_mesh


  subroutine interpolate_fc_test(mesh)
    type(unstr_mesh), intent(in) :: mesh
    !-
    integer :: j
    real(r8), allocatable :: xf(:,:), xc(:,:)
    real(r8) :: err

    associate(fc => mesh%face_centroid, cc => mesh%cell_centroid)
      allocate(xc(3,mesh%ncell))
      allocate(xf(3,mesh%nface))
      do j = 1, mesh%nface
        xf(:,j) = fc(:,j)
      end do

      call interpolate_fc(xc, xf)

      err = maxval(abs(cc(:,1:mesh%ncell_onP)-xc(:,1:mesh%ncell_onP)))
      print *, "max error for interpolate_fc: ", err
    end associate

  end subroutine interpolate_fc_test


  subroutine gradient_cf_test(id, mesh, init_field, flux_bc, dirichlet_bc, &
      solution_x, solution_y, solution_z)
    integer, intent(in) :: id
    type(unstr_mesh), intent(in) :: mesh
    class(scalar_func), intent(in) :: init_field, solution_x, solution_y, solution_z
    class(bndry_func1), intent(inout) :: dirichlet_bc, flux_bc
    !-
    real(r8), allocatable :: x(:), g(:,:)
    real(r8) :: args(0:3), linf_x, linf_y, linf_z
    integer :: j, i

    args(0) = 0.0_r8

    associate (cc => mesh%cell_centroid)
      allocate(x(mesh%ncell))
      allocate(g(3,mesh%nface))
      do i = 1, mesh%ncell
        args(1:3) = cc(:,i)
        x(i) = init_field%eval(args)
      end do
    end associate

    call gradient_cf(g, x, flux_bc, dirichlet_bc)

    ! check error
    linf_x = -huge(1.0_r8)
    linf_y = -huge(1.0_r8)
    linf_z = -huge(1.0_r8)

    associate (fc => mesh%face_centroid)
      do j = 1, mesh%nface_onP
        args(1:3) = fc(:,j)
        linf_x = max(linf_x, abs(g(1,j)-solution_x%eval(args)))
        linf_y = max(linf_y, abs(g(2,j)-solution_y%eval(args)))
        linf_z = max(linf_z, abs(g(3,j)-solution_z%eval(args)))
      end do
    end associate

    print *, "id: ", id
    print *, "linf_x: ", linf_x
    print *, "linf_y: ", linf_y
    print *, "linf_z: ", linf_z

!!$    call MPI_Reduce(linf, linf_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
!!$
!!$    if (this_rank==0) then
!!$      print *, 'linf: ', linf_global
!!$      if (linf_global > 1e-11_r8) status = 1
!!$    end if
!!$   call MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  end subroutine gradient_cf_test


    !! Handle a single optional argument which is the directory where the mesh
  !! files are located.  This is needed for cmake/ctest which puts and runs
  !! the executable from a different directory that this source code and mesh
  !! files.  Also handle '-h' or '--help' as an option to write out the usage.

  subroutine process_command_line (indir)

    use,intrinsic :: iso_fortran_env, only: output_unit

    character(:), allocatable, intent(out) :: indir

    character(:), allocatable :: prog
    character(256):: arg
    integer :: n, stat

    call get_command_argument (0, arg)
    n = scan(arg, '/', back=.true.) ! remove the leading path component, if any
    prog = trim(arg(n+1:))

    stat = 0
    select case (command_argument_count())
    case (0)
      indir = '.'
    case (1)
      call get_command_argument (1, arg)
      select case (arg)
      case ('-h','--help')
        stat = 1
      case default
        indir = trim(arg)
      end select
    case default
      stat = 1
    end select

    if (stat /= 0) then
      if (is_IOP) then
        write(output_unit,'(a)') 'Usage: ' // prog // ' [INDIR]'
        write(output_unit,'(a)') 'INDIR is the directory containing the input files (default ".")'
      end if
      call exit (stat) ! do not want this to count as a successful test
    end if

  end subroutine process_command_line

end program test_pressure_poisson
