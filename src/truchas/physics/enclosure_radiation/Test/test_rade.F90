!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

program test_rade

  use kinds, only: r8
  use pgslib_module
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use truchas_env, only: prefix, input_dir
  use truchas_logging_services
  use unstr_mesh_type
  use rad_problem_type
  implicit none

  type(unstr_mesh), pointer :: mesh
  character(len=PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:)
  real(r8), allocatable :: hc_temp(:)

  !! Enclosures to be tested
  character(len=31), parameter :: ENCL_EXT_FACE  = "exterior_face"
  character(len=31), parameter :: ENCL_INT_FACE  = "interior_face"
  character(len=31), parameter :: ENCL_EXT_PATCH = "exterior_patch"
  character(len=31), parameter :: ENCL_INT_PATCH = "interior_patch"

  !! Essential Truchas initialization.
  call parallel_init (argv)
  call init_parallel_communication ()

  !! Initialize Logging
  prefix='run'  ! TLS will write to 'run.log'
  call TLS_initialize ()
  call TLS_set_verbosity (TLS_VERB_NOISY)

  !! Initialize mesh and other Truchas state
  call init_state (mesh)

  !! Define an interesting temperature field over the mesh.
  call set_temp (mesh, hc_temp)

  !!!!!!!!!!!!!!!
  !! Run tests !!
  !!!!!!!!!!!!!!!

  !! Exterior enclosure
  call test_converge (mesh, ENCL_EXT_FACE, hc_temp)
  call test_converge (mesh, ENCL_EXT_PATCH, hc_temp)

  !! Interior enclosure
  call test_converge (mesh, ENCL_INT_FACE, hc_temp)
  call test_converge (mesh, ENCL_INT_PATCH, hc_temp)

  call halt_parallel_communication ()

contains


  !! Read Truchas input, set Truchas state, and initialize mesh.
  subroutine init_state (mesh)

    use enclosure_radiation_namelist, only: read_enclosure_radiation_namelists
    use mesh_manager
    use physical_constants, only: read_physical_constants
    use function_namelist, only: read_function_namelists

    type(unstr_mesh), pointer, intent(out) :: mesh

    character(:), allocatable :: indir
    integer, parameter :: lun = 10
    logical :: exists

    call process_command_line(indir)
    input_dir=indir // '/'
    if (is_IOP) open(lun,file=indir//'/test_rade.inp',position='rewind',action='read')

    !! Read essential Truchas namelists.
    call read_truchas_mesh_namelists (lun)
    call read_physical_constants (lun)
    call read_function_namelists (lun)

    !! Enclosure radiation namelists.
    call read_enclosure_radiation_namelists(lun)

    !! The mesh must be named 'main'
    call enable_mesh ('main', exists)
    INSIST(exists)

    !! Read and initialize the mesh.
    call init_mesh_manager ()
    mesh => unstr_mesh_ptr('main')
    INSIST(associated(mesh))

  end subroutine init_state


  !! Define an interesting temperature field over the mesh.
  subroutine set_temp(mesh, hc_temp)

    type(unstr_mesh), pointer, intent(in) :: mesh
    real(r8), allocatable :: hc_temp(:)

    real(r8) :: centroid(3)
    integer :: j

    allocate(hc_temp(mesh%nface))
    do j = 1, mesh%nface
      associate (face_nodes => mesh%fnode(mesh%xfnode(j):mesh%xfnode(j+1)-1))
        centroid = sum(mesh%x(:,face_nodes),dim=2) / size(face_nodes,dim=1)
        hc_temp(j) = 1.0_r8 + 2.0_r8 * sum(centroid**2)
      end associate
    end do

  end subroutine set_temp


  !! Test that radiosity solver converges within given error tolerance
  subroutine test_converge (mesh, encl_name, hc_temp)

    use enclosure_radiation_namelist, only: params
    use parameter_list_type

    type(unstr_mesh), pointer, intent(in) :: mesh
    character(len=31), intent(in) :: encl_name
    real(r8), intent(in) :: hc_temp(:)

    type(rad_problem) :: prob
    real(r8), allocatable :: temp(:), qrad(:), res(:)
    real(r8) :: res_norm, error, tol
    integer :: stat, numitr
    type(parameter_list), pointer :: plist

    if (is_IOP) write(*,'(/,a)') 'Solving radiosity system "' // trim(encl_name) // '" ...'

    plist => params%sublist(encl_name)
    call plist%get('error-tol', tol)
    call prob%init (mesh, encl_name, plist, tinit=0.0d0)  ! dummy time

    allocate(temp(prob%nface_hc), qrad(prob%nface_hc), res(prob%nface_hc))

    !! The given enclosure surface temperature.
    temp = hc_temp(prob%faces)

    !! Compute the surface radiosity.
    qrad = 0.0_r8
    call prob%solve_radiosity (0.0_r8, temp, qrad, stat, numitr, error)
    ASSERT(stat == 0)

    !! Compute the residual of the radiosity system
    call prob%residual  (0.0_r8, qrad, temp, res)
    res_norm = sqrt(global_sum(res**2))
    if (is_IOP) then
      print '(/,2x,a)',    'Radiosity solver results:'
      print '(4x,a,i13)',    'iterations    =', numitr
      print '(4x,a,es13.5)', 'error         =', error
      print '(4x,a,es13.5)', 'tolerance     =', tol
      print '(4x,a,es13.5)', 'residual norm =', res_norm
    end if

    !! Check error
    if (error > tol) then
      print '(/,a)', 'FAILURE: error for system "' // trim(encl_name) // '" exceeds tolerance.'
      stop 2
    end if

    deallocate(temp, res)

  end subroutine test_converge


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
      stop 1  ! do not want this to count as a successful test
    end if

  end subroutine process_command_line

end program test_rade
