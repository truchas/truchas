subroutine TEST_DIST(error)
  ! Test distribute routines in PGSLib

  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  ! local
  logical(pgslib_log_type) :: local_error

  local_error = .false.
  call test_dist_integer(local_error)
  error = error .or. local_error
  call test_dist_real(local_error)
  error = error .or. local_error
  call test_dist_double(local_error)
  error = error .or. local_error
  call test_dist_log(local_error)
  error = error .or. local_error

  return
end subroutine TEST_DIST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test distribute of integers
subroutine test_dist_integer(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error


  integer :: int_given, int_expected
  integer, pointer, dimension(:):: int_src, src_len, int_vec_given, int_vec_expected
  integer :: npe, n, thispe, pe, offset, local_size, src_size
  logical :: local_error,abort
  CHARACTER (LEN=1024) :: pgslib_out_string

  nPE = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()

  if (PGSLib_Inquire_thisPE() .eq. PGSLib_Inquire_IO_ROOT_PE()) then
     allocate(int_src(npe))
     int_src = (/ (n, n=1, npe) /)
  else
     allocate(int_src(1))
  endif

  int_expected = thispe
  int_given = -1

  call pgslib_dist(int_given, int_src)

  local_error = .false.
  if (int_given .ne. int_expected) then
     local_error = .true.
     print *, thispe, int_given, int_expected
     call pgslib_error("Failed integer scalar distribute test")
     call pgslib_output("Failed integer scalar distribute test")
  endif


  deallocate(int_src)
  ! Test distribute integer vector

  local_size = nPE - thisPE + 1
  allocate(int_vec_given(local_size))
  allocate(int_vec_expected(local_size))

  ! Set up expected result
  int_vec_given = -1
  int_vec_expected = (/ (n, n=local_size,1,-1) /)

  ! Set up source
  if (PGSLib_Inquire_thisPE() .eq. PGSLib_Inquire_IO_ROOT_PE()) then
     allocate(src_len(npe))
     src_len = (/ (npe - pe + 1, pe=1,npe) /)

     src_size = SUM( src_len ) ! = (npe*(npe+3))/2
     allocate(int_src(src_size))

     offset = 1
     do pe = 1, npe
        do n = npe-pe+1,1,-1
           int_src(offset) = n
           offset = offset + 1
        enddo
     enddo
  else
     allocate(src_len(1))
     allocate(int_src(1))
  endif

  call pgslib_dist(int_vec_given, int_src, src_len)

  if (ANY(int_vec_given .ne. int_vec_expected)) then
     local_error = .true.
     call pgslib_error("Failed integer vector distribute test")
     call pgslib_output("Failed integer vector distribute test")
  endif

  deallocate(int_src)
  deallocate(src_len)
  deallocate(int_vec_expected)
  deallocate(int_vec_given)

  if (.not. local_error) then
     call pgslib_output("Passed distribute integer test")
  endif
  error = error .or. local_error

  return
end subroutine test_dist_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test distribute of real
subroutine test_dist_real(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  real (pgslib_real_type) :: real_given, real_expected
  real (pgslib_real_type), pointer, dimension(:):: real_src, real_vec_given, real_vec_expected
  integer, pointer, dimension(:):: src_len
  integer :: npe, n, thispe, pe, offset, local_size, src_size
  logical :: local_error,abort

  nPE = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()

  allocate(real_src(npe))
  real_src = (/ (n, n=1, npe) /)
  real_src = sin(real_src)

  real_expected = thispe
  real_expected = sin(real_expected)
  real_given = -2

  call pgslib_dist(real_given, real_src)
  local_error = .false.
  if (real_given .ne. real_expected) then
     local_error = .true.
     call pgslib_error("Failed real scalar distribute test")
     call pgslib_output("Failed real scalar distribute test")
  endif

  deallocate(real_src)

  ! Test distribute real vector

  local_size = nPE - thisPE + 1
  allocate(real_vec_given(local_size))
  allocate(real_vec_expected(local_size))

  ! Set up expected result
  real_vec_given = -2
  real_vec_expected = (/ (n, n=local_size,1,-1) /)
  real_vec_expected = sin(real_vec_expected)

  ! Set up source
  allocate(src_len(npe))
  src_len = (/ (npe - pe + 1, pe=1,npe) /)
  
  src_size = SUM( src_len )      ! = (npe*(npe+3))/2
  allocate(real_src(src_size))

  offset = 1
  do pe = 1, npe
     do n = npe-pe+1,1,-1
        real_src(offset) = n
        offset = offset + 1
     enddo
  enddo
  real_src = sin(real_src)

  call pgslib_dist(real_vec_given, real_src, src_len)

  if (ANY(real_vec_given .ne. real_vec_expected)) then
     local_error = .true.
     call pgslib_error("Failed real vector distribute test")
     call pgslib_output("Failed real vector distribute test")
  endif

  if (ANY(ABS(real_vec_given - real_vec_expected) .gt. 1e-3)) then
     local_error = .true.
     call pgslib_error("Failed real vector distribute weak test")
     call pgslib_output("Failed real vector distribute weak test")
  endif

  deallocate(real_src)
  deallocate(src_len)
  deallocate(real_vec_expected)
  deallocate(real_vec_given)

  if (.not. local_error) then
     call pgslib_output("Passed distribute real test")
  endif

  error = error .or. local_error
  return
end subroutine test_dist_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test distribute of double
subroutine test_dist_double(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type):: error

  real (pgslib_double_type) :: double_given, double_expected
  real (pgslib_double_type), pointer, dimension(:):: double_src, double_vec_given, double_vec_expected
  integer, pointer, dimension(:):: src_len
  integer :: npe, n, thispe, pe, offset, local_size, src_size
  logical :: local_error,abort

  nPE = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()

  allocate(double_src(npe))
  double_src = (/ (n, n=1, npe) /)
  double_src = sin(double_src)

  double_expected = thispe
  double_expected = sin(double_expected)
  double_given = -2

  call pgslib_dist(double_given, double_src)

  local_error = .false.
  if (double_given .ne. double_expected) then
     local_error = .true.
     call pgslib_error("Failed double scalar distribute test")
     call pgslib_output("Failed double scalar distribute test")
  endif

  deallocate(double_src)
  ! Test distribute double vector

  local_size = nPE - thisPE + 1
  allocate(double_vec_given(local_size))
  allocate(double_vec_expected(local_size))

  ! Set up expected result
  double_vec_given = -2
  double_vec_expected = (/ (n, n=local_size,1,-1) /)
  double_vec_expected = sin(double_vec_expected)

  ! Set up source
  allocate(src_len(npe))
  src_len = (/ (npe - pe + 1, pe=1,npe) /)
  
  src_size = SUM( src_len )    ! = (npe*(npe+3))/2
  allocate(double_src(src_size))
  
  offset = 1
  do pe = 1, npe
     do n = npe-pe+1,1,-1
        double_src(offset) = n
        offset = offset + 1
     enddo
  enddo
  double_src = sin(double_src)

  call pgslib_dist(double_vec_given, double_src, src_len)

  if (ANY(double_vec_given .ne. double_vec_expected)) then
     local_error = .true.
     call pgslib_error("Failed double vector distribute test")
     call pgslib_output("Failed double vector distribute test")
  endif

  deallocate(double_src)
  deallocate(src_len)

  deallocate(double_vec_expected)
  deallocate(double_vec_given)

  if (.not. local_error) then
     call pgslib_output("Passed distribute double test")
  endif

  error = error .or. local_error
  return
end subroutine test_dist_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test distribute of log
subroutine test_dist_log(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type):: error

  logical (pgslib_log_type) :: log_given, log_expected
  logical (pgslib_log_type), pointer, dimension(:):: log_src, log_vec_given, log_vec_expected
  integer, pointer, dimension(:):: src_len
  integer :: npe, n, thispe, pe, offset, local_size, src_size
  logical :: local_error,abort

  nPE = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()

  allocate(log_src(npe))
  log_src = MOD((/ (n, n=1, npe) /),3) == 0

  log_expected = MOD(thispe,3) == 0
  log_given = .FALSE.

  call pgslib_dist(log_given, log_src)
  local_error = .false.
  if ((log_given .or. log_expected) .and. .not.(log_given .and. log_expected)) then ! == (log_g .xor. log_e)
     local_error = .true.
     call pgslib_error("Failed log scalar distribute test")
     call pgslib_output("Failed log scalar distribute test")
  endif

  deallocate(log_src)

  ! Test distribute log vector

  local_size = nPE - thisPE + 1
  allocate(log_vec_given(local_size))
  allocate(log_vec_expected(local_size))

  ! Set up expected result
  log_vec_given = .FALSE.
  log_vec_expected = MOD((/ (n, n=local_size,1,-1) /),3) == 0

  ! Set up source
  allocate(src_len(npe))
  src_len = (/ (npe - pe + 1, pe=1,npe) /)
  
  src_size = SUM( src_len )               ! = (npe*(npe+3))/2
  allocate(log_src(src_size))
  
  offset = 1
  do pe = 1, npe
     do n = npe-pe+1,1,-1
        log_src(offset) = MOD(n,3) == 0
        offset = offset + 1
     enddo
  enddo
  

  call pgslib_dist(log_vec_given, log_src, src_len)

  if ( ANY((log_vec_given .or. log_vec_expected) .and. &
       &    .not.(log_vec_given .and. log_vec_expected))  )then ! == (log_g .xor. log_e)     local_error = .true.
     local_error = .true.
     call pgslib_error("Failed log vector distribute test")
     call pgslib_output("Failed log vector distribute test")
  endif

  deallocate(log_src)
  deallocate(src_len)
  deallocate(log_vec_expected)
  deallocate(log_vec_given)

  if (.not. local_error) then
     call pgslib_output("Passed distribute log test")
  endif

  error = error .or. local_error
  return
end subroutine test_dist_log
