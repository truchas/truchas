!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_rade_tools

  use kinds, only: i4, r8
  use re_utilities
  use re_dist_vf_type
  use re_patch_type
  use re_encl_type
  use string_utilities
  use scl
  implicit none

  integer :: ios
  logical :: is_IOP
  type(encl) :: e
  type(dist_vf) :: dvf
  type(re_patch) :: ep
  character(len=512) :: infile
  real(r8) :: tolerance
  logical :: vf_file

  call scl_init ()
  is_IOP = (scl_rank()==1)

  !! Read command line
  if (is_IOP) then
    call process_command_line(infile, vf_file)
    if (vf_file) open(unit=10,file=trim(infile),status='old',action='read',iostat=ios)
  endif
  call scl_bcast(vf_file)
  call scl_bcast(ios)
  if (vf_file .and. ios /= 0) call re_halt('Unable to open input file: ' // trim(infile))

  !! Get test data
  if (vf_file) then
    call re_info('Running tests with input file: ' // trim(infile))
    call init_data_from_file(dvf, ep, infile)
    tolerance = 1E-6
  else
    call re_info('Running tests with synthetic data')
    if (scl_size() > 1) call re_halt('Only serial runs are supported.')
    call init_data(dvf, ep)
    tolerance = 1E-7
  end if

  !! Run tests
  call test_get_column(dvf, tolerance)

  call re_info('All tests completed.')
  call scl_finalize()

contains


  !! Initialize data for TEST_GET_COLUMN.  The data is chosen
  !! so that that all TEST_GET_COLUMN computations are exact.
  subroutine init_data (dvf, ep)

    type(dist_vf), intent(out) :: dvf
    type(re_patch), intent(out) :: ep

    integer, parameter :: NFACE = 128
    integer, parameter :: NPATCH = NFACE/4
    real, allocatable :: vf_face(:,:), vf_patch(:,:)  ! Dense VF matrices
    real, allocatable :: parea(:)  ! patch areas
    integer :: i, j, n
    real :: r

    allocate(vf_face(NFACE,NFACE), vf_patch(NPATCH,NPATCH))
    allocate(dvf%area(NFACE), parea(NPATCH))

    !! F2P map with four faces per patch
    allocate(ep%f2p_map(NFACE))
    ep%npatch = NPATCH
    ep%has_patches = .true.
    do i = 1, NFACE
      ep%f2p_map(i) = (i+3)/4
    end do
    call shuffle(ep%f2p_map)

    !! Compute patch and face areas
    call random_exact(parea)
    do i = 1, NFACE
      dvf%area(i) = parea(ep%f2p_map(i)) / 4
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
        vf_face(j,i) = vf_patch(ep%f2p_map(j),ep%f2p_map(i)) * (dvf%area(j)/parea(ep%f2p_map(j)))
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Pack face-based data !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    dvf%offset = 0
    dvf%npatch = NFACE
    dvf%npatch_tot = NFACE

    !! Generate IA indexing array
    allocate(dvf%ia(NFACE+1))
    dvf%ia(1) = 1
    do i = 1, NFACE
       dvf%ia(i+1) = dvf%ia(i) + count(vf_face(:,i) /= 0.0)
    end do

    n = dvf%ia(NFACE+1)-1
    allocate(dvf%ja(n), dvf%val(n), dvf%ambient(0))

    !! Generate sparse matrix rows
    n = 1
    do i = 1, NFACE
      do j = 1, NFACE
        if (vf_face(j,i) /= 0.0) then
          dvf%val(n) = vf_face(j,i)
          dvf%ja(n) = j
          n = n + 1
        end if
      end do
    end do

  end subroutine init_data


  !! Initialize data for TEST_GET_COLUMN.  The data is chosen
  !! so that that all TEST_GET_COLUMN computations are exact.
  subroutine init_data_from_file (dvf, ep, infile)

    type(dist_vf), intent(out) :: dvf
    type(re_patch), intent(out) :: ep
    character(len=512), intent(in) :: infile

    logical :: has_vf

    call e%read(trim(infile), has_vf)
    if (.not. has_vf) call re_halt ('No view factor data in file: ' // trim(infile))

    if (is_IOP) call ep%read_patch_data(trim(infile))
    call read_dist_vf(dvf, trim(infile))

  end subroutine init_data_from_file


  !! Test that dist_vf gets the correct column values
  subroutine test_get_column (dvf, tol)

    type(dist_vf), intent(in) :: dvf
    real(r8), intent(in) :: tol

    real, allocatable :: col1(:), col2(:)
    real :: abs_err, max_abs_err, min_abs_err, rel_err, max_rel_err, min_rel_err
    integer :: i, j, cnt
    character(len=64) :: string

    cnt = 0
    max_abs_err = 0.0
    min_abs_err = huge(0.0)
    max_rel_err = 0.0
    min_rel_err = huge(0.0)

    do i = 1, dvf%npatch_tot
      col1 = unpack_dvf_col(dvf, i)
      col2 = unpack_dvf_col_explicit(dvf, i)

      if (is_IOP) then
        do j = 1, size(col1)
          abs_err = abs(col1(j) - col2(j))
          if (abs_err == 0.0) cycle
          rel_err = abs_err / col2(j)
          if (abs_err > max_abs_err) max_abs_err = abs_err
          if (abs_err < min_abs_err) min_abs_err = abs_err
          if (rel_err > max_rel_err) max_rel_err = rel_err
          if (rel_err < min_rel_err) min_rel_err = rel_err
          if (rel_err > tol) cnt = cnt + 1
        end do
      end if
    end do

    if (is_IOP .and. cnt /= 0) then
      print '("Mismatched entries: ", i7)', cnt
      print '("  Percent of nonzeros: ", f7.3,"% (", i8," total)")', cnt/REAL(size(dvf%val))*100.0, size(dvf%val)
      print '("  Percent of matrix:   ", f7.3,"% (", i8," total)")', cnt/REAL(dvf%npatch_tot*dvf%npatch_tot)*100.0, dvf%npatch_tot*dvf%npatch_tot
      print '("Max absolute error: ", es11.4)', max_abs_err
      print '("Min absolute error: ", es11.4)', min_abs_err
      print '("Max relative error: ", es11.4)', max_rel_err
      print '("Min relative error: ", es11.4)', min_rel_err
    end if

    call scl_bcast(cnt)
    if (cnt /= 0) then
      write (string,'(i7, " mismatched entries with tolerance", es10.3)') cnt, tol
      call re_halt(trim(string))
    else
      write (string,'("SUCCESS: No mismatched entries with tolerance", es10.3)') tol
      call re_info(trim(string))
    end if

  end subroutine test_get_column


  !! Randomly permutes an array using th Fisher-Yates algorithm
  !! Source: https://www.rosettacode.org/wiki/Knuth_shuffle#Fortran
  subroutine shuffle(a)

    integer, intent(inout) :: a(:)

    integer :: i, randpos, temp
    real :: r

    do i = size(a), 2, -1
      call random_number(r)
      randpos = int(r * i) + 1
      temp = a(randpos)
      a(randpos) = a(i)
      a(i) = temp
    end do

  end subroutine shuffle


  !! Returns a random exactly-representable real number.
  !! Rationals of the form a/2^k, where k is less than the
  !! number of the exponent bits, are exactly representatble.
  elemental impure subroutine random_exact(rand)

    real, intent(out) :: rand

    !! We only want fractions of the form 1/2^k
    integer(i4), parameter :: mask = Z'3F800000'
    integer(i4) :: i
    real :: r

    call random_number(r)
    i = transfer(r, i)
    i = iand(i, mask)
    rand = transfer(i, r)

  end subroutine random_exact


  !! Handle a optional single argument which is the radiation enclosure file to use
  !! for testing. Also handle '-h' or '--help' as an option to write out the usage.
  subroutine process_command_line (infile, read_file)

    character(len=512), intent(out) :: infile
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

end program test_rade_tools
