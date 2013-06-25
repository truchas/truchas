subroutine TEST_BCAST(error)
  ! Test broadcast routines in PGSLib

  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error
  
  ! local variables
  logical (pgslib_log_type) :: local_error

  local_error = .false.
  call test_bcast_integer(local_error)
  error = error .or. local_error
  call test_bcast_real(local_error)
  error = error .or. local_error
  call test_bcast_double(local_error)
  error = error .or. local_error
  call test_bcast_log(local_error)
  error = error .or. local_error
  call test_bcast_char(local_error)
  error = error .or. local_error

  return
end subroutine TEST_BCAST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test broadcast of integers
subroutine test_bcast_integer(routine_error)
  USE PGSLib_MODULE
  implicit none
  logical (pgslib_log_type) :: routine_error

  integer :: int_given, int_expected
  integer, pointer, dimension(:) :: int_vec_given, int_vec_expected
  integer, allocatable, dimension(:,:,:) :: int_3d_given, int_3d_expected
  integer :: nPE, n, i, j, k
  logical :: local_error, abort
  

  ! Test broadcast integer scalar
  int_expected = 100
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     int_given = int_expected
  else
     int_given = - int_expected
  endif
  
  call pgslib_bcast(int_given)
  local_error = .false.
  if (int_given .ne. int_expected) then
     local_error = .true.
     call pgslib_output("Failed integer scalar broadcast test")
  endif

  
  abort = pgslib_global_any(local_error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif
  routine_error = abort

  ! Test broadcast integer vector
  nPE = PGSLib_Inquire_nPE()
  allocate(int_vec_given(nPE))
  allocate(int_vec_expected(nPE))

  int_vec_expected = (/ (n, n=0,npe-1,1) /)
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     int_vec_given = int_vec_expected
  else
     int_vec_given = - int_vec_expected
  endif

  call pgslib_bcast(int_vec_given)
  local_error = .false.
  if ( ANY(int_vec_given .ne. int_vec_expected) ) then
     local_error = .true.
     call pgslib_output("Failed integer vector broadcast test")
  endif
  routine_error = routine_error .or. local_error

  deallocate(int_vec_given)
  deallocate(int_vec_expected)

  ! Test 3d vector broadcast
  nPE = PGSLib_Inquire_nPE()
  allocate(int_3d_given(npe, 11, npe*3))
  allocate(int_3d_expected(npe, 11, npe*3))

  do k=1,size(int_3d_expected,3)
     do j=1,size(int_3d_expected,2)
        do i = 1,size(int_3d_expected,1)
           int_3d_expected(i,j,k) = i+2*j*5 +k
        end do
     end do
  end do
  
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     int_3d_given = int_3d_expected
  else
     int_3d_given = - int_3d_expected
  end if
  
  call pgslib_bcast(int_3d_given)
  local_error = .false.
  if (ANY(int_3d_given .ne. int_3d_expected) ) then
     local_error = .true.
     call pgslib_output("Faile integer 3d-vector broadcast test")
  endif
  routine_error = routine_error .or. local_error
  
  deallocate(int_3d_expected)
  deallocate(int_3d_given)
  
  abort = pgslib_global_any(local_error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif

  if (.not. abort) then
     call pgslib_output("Passed broadcast integer test")
  endif

  return
end subroutine test_bcast_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_bcast_real(routine_error)
  USE PGSLib_MODULE
  implicit none
  logical (pgslib_log_type) :: routine_error


  real (pgslib_real_type) :: real_given, real_expected
  real (pgslib_real_type), pointer, dimension(:) :: real_vec_given, real_vec_expected
  real (pgslib_real_type), pointer, dimension(:,:,:) :: real_3d_given, real_3d_expected
  integer :: nPE, n, i, j, k
  logical :: error, abort
  

  ! Test broadcast real scalar
  real_expected = 1.2389
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     real_given = real_expected
  else
     real_given = - real_expected
  endif
  
  call pgslib_bcast(real_given)
  error = .false.
  if (real_given .ne. real_expected) then
     error = .true.
     call pgslib_output("Failed real scalar broadcast test")
  endif

  
  abort = pgslib_global_any(error)
  routine_error = abort
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif

  ! Test broadcast real vector
  nPE = PGSLib_Inquire_nPE()
  allocate(real_vec_given(nPE))
  allocate(real_vec_expected(nPE))

  real_vec_expected = (/ (n, n=0,npe-1,1) /)
  real_vec_expected = sin(real_vec_expected)
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     real_vec_given = real_vec_expected
  else
     real_vec_given = - real_vec_expected
  endif

  call pgslib_bcast(real_vec_given)
  error = .false.
  if ( ANY(real_vec_given .ne. real_vec_expected) ) then
     error = .true.
     call pgslib_output("Failed real vector broadcast test")
  endif
  routine_error = routine_error .or. error

  deallocate(real_vec_given)
  deallocate(real_vec_expected)

  ! Test 3d vector broadcast
  nPE = PGSLib_Inquire_nPE()
  allocate(real_3d_given(npe, 11, npe*3))
  allocate(real_3d_expected(npe, 11, npe*3))

  do k=1,size(real_3d_expected,3)
     do j=1,size(real_3d_expected,2)
        do i = 1,size(real_3d_expected,1)
           real_3d_expected(i,j,k) = i+2*j*5 +k
        end do
     end do
  end do
  
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     real_3d_given = real_3d_expected
  else
     real_3d_given = - real_3d_expected
  end if
  
  call pgslib_bcast(real_3d_given)
  error = .false.
  if (ANY(real_3d_given .ne. real_3d_expected) ) then
     error = .true.
     call pgslib_output("Faile real 3d-vector broadcast test")
  endif
  routine_error = routine_error .or. error

  deallocate(real_3d_expected)
  deallocate(real_3d_given)
  
  abort = pgslib_global_any(error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif


  if (.not. abort) then
     call pgslib_output("Passed broadcast real test")
  endif

  return
end subroutine test_bcast_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_bcast_double(routine_error)
  USE PGSLib_MODULE
  implicit none

  logical (pgslib_log_type) :: routine_error

  real (pgslib_double_type) :: double_given, double_expected
  real (pgslib_double_type), pointer, dimension(:) :: double_vec_given, double_vec_expected
  real (pgslib_double_type), allocatable, dimension(:,:,:) :: double_3d_given, double_3d_expected
  integer :: nPE, n, i, j,k
  logical :: error, abort
  

  ! Test broadcast double scalar
  double_expected = 1.2389
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     double_given = double_expected
  else
     double_given = - double_expected
  endif
  
  call pgslib_bcast(double_given)
  error = .false.
  if (double_given .ne. double_expected) then
     error = .true.
     call pgslib_output("Failed double scalar broadcast test")
  endif

  routine_error = error

  abort = pgslib_global_any(error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif

  ! Test broadcast double vector
  nPE = PGSLib_Inquire_nPE()
  allocate(double_vec_given(nPE))
  allocate(double_vec_expected(nPE))

  double_vec_expected = (/ (n, n=0,npe-1,1) /)
  double_vec_expected = sin(double_vec_expected)
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     double_vec_given = double_vec_expected
  else
     double_vec_given = - double_vec_expected
  endif

  call pgslib_bcast(double_vec_given)
  error = .false.
  if ( ANY(double_vec_given .ne. double_vec_expected) ) then
     error = .true.
     call pgslib_output("Failed double vector broadcast test")
  endif

  routine_error = routine_error .or. error

  deallocate(double_vec_given)
  deallocate(double_vec_expected)

  ! Test 3d vector broadcast
  nPE = PGSLib_Inquire_nPE()
  allocate(double_3d_given(npe, 11, npe*3))
  allocate(double_3d_expected(npe, 11, npe*3))

  do k=1,size(double_3d_expected,3)
     do j=1,size(double_3d_expected,2)
        do i = 1,size(double_3d_expected,1)
           double_3d_expected(i,j,k) = i+2*j*5 +k
        end do
     end do
  end do
  
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     double_3d_given = double_3d_expected
  else
     double_3d_given = - double_3d_expected
  end if
  
  call pgslib_bcast(double_3d_given)
  error = .false.
  if (ANY(double_3d_given .ne. double_3d_expected) ) then
     error = .true.
     call pgslib_output("Faile double 3d-vector broadcast test")
  endif
  routine_error = routine_error .or. error
  
  deallocate(double_3d_expected)
  deallocate(double_3d_given)
  
  abort = pgslib_global_any(error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif


  if (.not. abort) then
     call pgslib_output("Passed broadcast double test")
  endif

  return
end subroutine test_bcast_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_bcast_log(routine_error)
  USE PGSLib_MODULE
  implicit none
  logical (pgslib_log_type) :: routine_error


  logical (pgslib_log_type) :: log_given, log_expected
  logical (pgslib_log_type), pointer, dimension(:) :: log_vec_given, log_vec_expected
  integer :: nPE, n
  logical :: error, abort
  

  ! Test broadcast log scalar
  log_expected = .false.
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     log_given = log_expected
  else
     log_given = .not. log_expected
  endif
  
  call pgslib_bcast(log_given)
  error = .false.
  if ((log_given .or. log_expected) .and. .not.(log_given .and. log_expected)) then ! == (log_g .xor. log_e)
     error = .true.
     call pgslib_output("Failed log scalar broadcast test")
  endif

  routine_error = error
  
  abort = pgslib_global_any(error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif

  ! Test broadcast log vector
  nPE = PGSLib_Inquire_nPE()
  allocate(log_vec_given(nPE))
  allocate(log_vec_expected(nPE))

  log_vec_expected = MOD((/ (n, n=0,npe-1,1) /),2) == 0
  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     log_vec_given = log_vec_expected
  else
     log_vec_given = .not.log_vec_expected
  endif

  call pgslib_bcast(log_vec_given)
  error = .false.
  if ( ANY((log_vec_given .or. log_vec_expected) .and. &
       &    .not.(log_vec_given .and. log_vec_expected))  )then ! == (log_g .xor. log_e)     error = .true.
     error = .true.
     call pgslib_output("Failed log vector broadcast test")
  endif
  routine_error = routine_error .or. error


  deallocate(log_vec_given)
  deallocate(log_vec_expected)

  abort = pgslib_global_any(error)
  if (abort) then
     call pgslib_fatal_error("Failed Broadcast test")
  endif


  if (.not. abort) then
     call pgslib_output("Passed broadcast log test")
  endif

  return
end subroutine test_bcast_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_bcast_char(routine_error)
  USE PGSLib_MODULE
  implicit none
  logical (pgslib_log_type) :: routine_error

  character (LEN=100) :: char_vec_given, char_vec_expected
  character (LEN=256), allocatable, dimension(:) :: char_2d_given, char_2d_expected
  character (LEN=1024) :: root_string
  integer :: nPE, n, nstrings, i
  logical :: error, abort
  

  ! Test broadcast char string
  char_vec_given = ''
  char_vec_expected = ''
  char_vec_expected = 'This is the expected string.'

  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     char_vec_given = char_vec_expected
  else
     char_vec_given = ''
  endif

  call pgslib_bcast(char_vec_given)
  error = .false.
  if (char_vec_given .ne. char_vec_expected) then
     call pgslib_output("Failed char vector broadcast test")
     error = .true.
  endif
  routine_error = error

  abort = pgslib_global_any(error)

  ! Test bcast of vector of strings
  nstrings = 3
  allocate(char_2d_given(nstrings))
  allocate(char_2d_expected(nstrings))

  root_string = 'The root of all evil'
  do i=1,nstrings
     char_2d_expected(i) = REPEAT(TRIM(root_string),NCOPIES=i)
  end do

  if (PGSLib_Inquire_thisPE() == PGSLib_Inquire_IO_ROOT_PE() ) then
     char_2d_given = char_2d_expected
  else
     char_2d_given = char_2d_expected(nstrings)
  endif

  call pgslib_bcast(char_2d_given)

  error = ANY(char_2d_given .ne. char_2d_expected) 
  if (error) then
     call pgslib_output("Failed char vector of strings broadcast test.")
  end if
  routine_error = routine_error .or. error
  ! Completed testing character broadcasts

  abort = pgslib_global_any(error)
  if (abort) then
!!$     call pgslib_fatal_error("Failed Broadcast test")
     call pgslib_output("Failed Broadcast test")
  endif

  if (.not. abort) then
     call pgslib_output("Passed broadcast char test")
  endif

  return
end subroutine test_bcast_char

