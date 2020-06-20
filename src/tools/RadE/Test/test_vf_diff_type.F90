!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_vf_diff_type

  use kinds, only: r8
  use test_rade_tools_common
  use re_utilities
  use re_dist_vf_type
  use re_patch_type
  use re_vf_diff_type
  use string_utilities
  use scl
  implicit none

  integer :: ios
  logical :: is_IOP
  type(vf_diff) :: vfd
  character(:), allocatable :: infile
  real(r8) :: tolerance
  logical :: vf_file

  call scl_init()
  is_IOP = (scl_rank()==1)
  tolerance = 1E-14

  !! Read command line
  if (is_IOP) then
    call process_command_line(infile, vf_file)
    if (vf_file) open(unit=10,file=infile,status='old',action='read',iostat=ios)
  else
    allocate(infile, source='')
  endif
  call scl_bcast(vf_file)
  call scl_bcast(ios)
  if (vf_file .and. ios /= 0) call re_halt('Unable to open input file: ' // infile)

  !! Get test data
  if (vf_file) then
    call re_info('Running tests with input file: ' // infile)
    call init_data_from_file(vfd, infile)
  else
    call re_info('Running tests with synthetic data')
    if (scl_size() > 1) call re_halt('Only serial runs are supported.')
    call init_data(vfd)
  end if

  !! Run tests
  call test_zero_norms(vfd, tolerance)

  call re_info('All tests completed.')
  call scl_finalize()

contains


  !! Initialize a VF_DIFF type for testing. The internal matrices A and B are
  !! patch and face based, respectively. The VF_DIFF data is chosen so that
  !! all computations are exact.
  subroutine init_data (vfd)

    type(vf_diff), intent(out) :: vfd

    integer, parameter :: NFACE = 128
    integer, parameter :: NPATCH = NFACE/4
    real, allocatable :: vf_face(:,:), vf_patch(:,:)  ! Dense VF matrices
    real, allocatable :: parea(:)  ! patch areas
    real, allocatable :: farea(:)  ! face areas
    integer, allocatable :: f2p_map(:)
    integer :: i, j, n
    real :: r

    allocate(vf_face(NFACE,NFACE), vf_patch(NPATCH,NPATCH))
    allocate(farea(NFACE), parea(NPATCH))
    allocate(f2p_map(NFACE))

    !! F2P map with four faces per patch
    do i = 1, NFACE
      f2p_map(i) = (i+3)/4
    end do
    call shuffle(f2p_map)

    !! Compute patch and face areas
    call random_exact(parea)
    do i = 1, NFACE
      farea(i) = parea(f2p_map(i)) / 4
    end do

    !! Generate tridiagonal patch-based matrix.
    !! For efficient packing into DVF we set vf_patch(:,i) as the i-th row
    vf_patch = 0.0
    do i = 1, NPATCH-1
      call random_exact(r)
      vf_patch(i+1,i) = r
      vf_patch(i,i+1) = r * (parea(i) / parea(i+1))
    end do

    !! Construct face-based matrix.
    do i = 1, NFACE
      do j = 1, NFACE
        vf_face(j,i) = vf_patch(f2p_map(j),f2p_map(i)) * (farea(j)/parea(f2p_map(j)))
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!
    !! Initialize VF_DIFF !!
    !!!!!!!!!!!!!!!!!!!!!!!!
    vfd%nface_tot = NFACE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Initialize patch matrix A !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Enclosure patches component
    vfd%epA%npatch = NPATCH
    vfd%epA%nface  = NFACE
    vfd%epA%has_patches = .true.
    allocate(vfd%epA%f2p_map, source=f2p_map)

    !! Packed VF matrix component
    vfd%A%offset = 0
    vfd%A%npatch = NPATCH
    vfd%A%npatch_tot = NPATCH
    vfd%A%has_ambient = .false.
    vfd%A%has_area = .true.
    allocate(vfd%A%area(NPATCH))
    vfd%A%area = parea

    !! Face weights
    allocate(vfd%A%w(NFACE))
    do i = 1, NFACE
      vfd%A%w(i) = farea(i) / parea(f2p_map(i))
    end do

    !! Generate IA indexing array
    allocate(vfd%A%ia(NPATCH+1))
    vfd%A%ia(1) = 1
    do i = 1, NPATCH
       vfd%A%ia(i+1) = vfd%A%ia(i) + count(vf_patch(:,i) /= 0.0)
    end do

    n = vfd%A%ia(NPATCH+1)-1
    allocate(vfd%A%ja(n), vfd%A%val(n), vfd%A%ambient(0))

    !! Generate sparse matrix rows
    n = 1
    do i = 1, NPATCH
      do j = 1, NPATCH
        if (vf_patch(j,i) /= 0.0) then
          vfd%A%val(n) = vf_patch(j,i)
          vfd%A%ja(n) = j
          n = n + 1
        end if
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Initialize face matrix B !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Enclosure patches component
    vfd%epB%npatch = NFACE
    vfd%epB%nface  = NFACE
    vfd%epB%has_patches = .false.
    vfd%epB%f2p_map = [(i,i=1,NFACE)]

    !! Packed VF matrix component
    vfd%B%offset = 0
    vfd%B%npatch = NFACE
    vfd%B%npatch_tot = NFACE
    vfd%B%has_ambient = .false.
    vfd%B%has_area = .true.
    allocate(vfd%B%area(NFACE))
    vfd%B%area = farea

    !! Face weights
    allocate(vfd%B%w(NFACE), source=1.0_r8)

    !! Generate IB indexing array
    allocate(vfd%B%ia(NFACE+1))
    vfd%B%ia(1) = 1
    do i = 1, NFACE
       vfd%B%ia(i+1) = vfd%B%ia(i) + count(vf_face(:,i) /= 0.0)
    end do

    n = vfd%B%ia(NFACE+1)-1
    allocate(vfd%B%ja(n), vfd%B%val(n), vfd%B%ambient(0))

    !! Generate sparse matrix rows
    n = 1
    do i = 1, NFACE
      do j = 1, NFACE
        if (vf_face(j,i) /= 0.0) then
          vfd%B%val(n) = vf_face(j,i)
          vfd%B%ja(n) = j
          n = n + 1
        end if
      end do
    end do

  end subroutine init_data


  !! Initialize a VF_DIFF type for testing. The internal matrices A and B are
  !! patch two copies of the same input file.
  subroutine init_data_from_file (vfd, infile)

    type(vf_diff), intent(out) :: vfd
    character(:), allocatable, intent(in) :: infile

    call vfd%init(infile, infile)

  end subroutine init_data_from_file


  !! Test that all norms of the difference matrix D = B - A are zero.
  subroutine test_zero_norms (vfd, tol)

    type(vf_diff), intent(in) :: vfd
    real(r8), intent(in) :: tol

    real :: max_norm, one_norm, two_norm, max_diff
    integer :: rloc, cloc, cnt
    character(len=64) :: string

    call vfd%compute_norms(max_norm, one_norm, two_norm, max_diff, rloc, cloc)

    if (is_IOP) print '(3(a,es10.3))', 'max norm =', max_norm, '; one norm =', one_norm, '; two norm =', two_norm
    if (is_IOP) print '(a,es10.3,a)', 'max difference =', max_diff, &
          ' at row ' // i_to_c(rloc) // ', col ' // i_to_c(cloc)

    cnt = 0
    if (max_norm > tol) cnt = cnt + 1
    if (one_norm > tol) cnt = cnt + 1
    if (two_norm > tol) cnt = cnt + 1
    if (max_diff > tol) cnt = cnt + 1

    if (cnt /= 0) then
      write (string,'(i7, " norm(s) are above the tolerance of", es10.3)') cnt, tol
      call re_halt(trim(string))
    else
      write (string,'("SUCCESS: All norms are below tolerance of", es10.3)') tol
      call re_info(trim(string))
    end if

  end subroutine test_zero_norms

  !! Handle a optional single argument which is the radiation enclosure file to use
  !! for testing. Also handle '-h' or '--help' as an option to write out the usage.
  subroutine process_command_line (infile, read_file)

    character(:), allocatable, intent(out) :: infile
    logical, intent(out) :: read_file

    character(:), allocatable :: prog
    character(len=256):: arg
    integer :: n, stat

    call get_command_argument (0, arg)
    n = scan(arg, '/', back=.true.) ! remove the leading path component, if any
    prog = trim(arg(n+1:))

    stat = 0
    select case (command_argument_count())
    case (0)
      read_file = .false.
    case (1)
      call get_command_argument (1, arg)
      select case (arg)
      case ('-h','--help')
        stat = 1
      case default
        read_file = .true.
        infile = trim(arg)
      end select
    case default
      stat = 1
    end select

    if (stat /= 0) then
      if (is_IOP) then
        print *, 'Usage: ' // prog // ' [INFILE]'
        print *, 'INFILE is the enclosure radiation file to use for testing.'
        print *, 'If INFILE is not specified, the test is run on synthetic data.'
      end if
      stop 1  ! do not want this to count as a successful test
    end if

  end subroutine process_command_line

end program test_vf_diff_type
