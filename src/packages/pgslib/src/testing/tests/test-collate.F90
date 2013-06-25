subroutine TEST_COLLATE(error)
  ! Test collate routines in PGSLib

  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  !local variables
  logical(pgslib_log_type) :: local_error

  local_error = .false.
  call test_collate_int(local_error)
  error = error .or. local_error

  local_error = .false.
  call test_collate_real(local_error)
  error = error .or. local_error

  local_error = .false.
  call test_collate_double(local_error)
  error = error .or. local_error

  local_error = .false.
  call test_collate_char(local_error)
  error = error .or. local_error

  local_error = .false.
  call test_collate_log(local_error)
  error = error .or. local_error

  return
end subroutine TEST_COLLATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test collate of integers
subroutine test_collate_int(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  integer(pgslib_int_type) :: int_given, int_src
  integer(pgslib_int_type), pointer, dimension(:):: int_vec_src, int_vec_given, int_vec_expected
  integer(pgslib_int_type):: npe, n, thispe, pe, offset, local_size, dest_size
  logical(pgslib_log_type) :: local_error, abort

  npe = PGSLib_Inquire_nPE()
  thisPE = PGSLib_Inquire_thisPE_Actual()

  ! Set up for scalar test
  allocate(int_vec_given(nPE))
  int_vec_given = -1
  allocate(int_vec_expected(npe))
  int_vec_expected = (/ (n, n=npe, 1,-1) /)

  int_src = npe - (thispe -1)

  call pgslib_collate(int_vec_given, int_src)

  local_error = .false.
  if (PGSLib_Inquire_IO_P() ) then
     if (ANY(int_vec_given .ne. int_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed integer scalar collate test")
        call pgslib_output("Failed integer scalar collate test")
     endif
  endif
      
  local_error = PGSLib_Global_ANY(local_error)

  deallocate(int_vec_expected)
  deallocate(int_vec_given)
  ! Set up for vector test
  local_size = thispe*2 + 1
  allocate(int_vec_src(local_size))
  int_vec_src = (/ (n, n=thispe*2+1,1,-1) /)

  dest_size = SUM( (/ (pe*2+1, pe=1,npe) /) )
  allocate(int_vec_given(dest_size))
  allocate(int_vec_expected(dest_size))
  
  int_vec_given = -1
  
  offset = 1
  do pe=1,npe
     do n=pe*2+1,1,-1
        int_vec_expected(offset) = n
        offset = offset + 1
     enddo
  enddo

  call pgslib_collate(int_vec_given, int_vec_src)

  if (PGSLib_Inquire_IO_P() ) then
     if (ANY(int_vec_given .ne. int_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed integer vector collate test")
        call pgslib_output("Failed integer vector collate test")
     endif
  endif

  local_error = PGSLib_Global_ANY(local_error)

  deallocate(int_vec_expected)
  deallocate(int_vec_given)
  deallocate(int_vec_src)

  if (.not. local_error) then
     call pgslib_output("Passed collate integer test")
  endif
  error = error .or. local_error

  return
end subroutine test_collate_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test collate of reals
subroutine test_collate_real(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  real (pgslib_real_type) :: real_given, real_src
  real (pgslib_real_type), pointer, dimension(:):: real_vec_src, real_vec_given, real_vec_expected
  integer(pgslib_int_type):: npe, n, thispe, pe, offset, local_size, dest_size
  logical(pgslib_log_type) :: local_error, abort

  npe = PGSLib_Inquire_nPE()
  thisPE = PGSLib_Inquire_thisPE_Actual()

  ! Set up for scalar test
  allocate(real_vec_given(nPE))
  real_vec_given = -1
  allocate(real_vec_expected(npe))
  real_vec_expected = (/ (n, n=npe, 1,-1) /)
  real_vec_expected = sin(real_vec_expected)

  real_src = npe - (thispe - 1)
  real_src = sin(real_src)

  call pgslib_collate(real_vec_given, real_src)

  local_error = .false.
  if (PGSLib_Inquire_thisPE() .eq. PGSLib_Inquire_IO_ROOT_PE()) then
     if (ANY(real_vec_given .ne. real_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed real scalar collate test")
        call pgslib_output("Failed real scalar collate test")
     endif
  endif

  deallocate(real_vec_expected)
  deallocate(real_vec_given)

  ! Set up for vector test
  local_size = thispe*2 + 1
  allocate(real_vec_src(local_size))
  real_vec_src = (/ (n, n=thispe*2+1,1,-1) /)
  real_vec_src = sin(real_vec_src)

  dest_size = SUM( (/ (pe*2+1, pe=1,npe) /) )
  allocate(real_vec_given(dest_size))
  allocate(real_vec_expected(dest_size))
  
  real_vec_given = -1
  
  offset = 1
  do pe=1,npe
     do n=pe*2+1,1,-1
        real_vec_expected(offset) = n
        offset = offset + 1
     enddo
  enddo
  real_vec_expected = sin(real_vec_expected)

  call pgslib_collate(real_vec_given, real_vec_src)

  if (PGSLib_Inquire_thisPE() .eq. PGSLib_Inquire_IO_ROOT_PE()) then
     if (ANY(real_vec_given .ne. real_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed real vector collate test")
        call pgslib_output("Failed real vector collate test")
     endif
  endif

  deallocate(real_vec_expected)
  deallocate(real_vec_given)
  deallocate(real_vec_src)

  if (.not. local_error) then
     call pgslib_output("Passed collate real test")
  endif
  error = error .or. local_error

  return
end subroutine test_collate_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test collate of doubles
subroutine test_collate_double(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  real (pgslib_double_type) :: double_given, double_src
  real (pgslib_double_type), pointer, dimension(:):: double_vec_src, double_vec_given, double_vec_expected
  integer(pgslib_int_type):: npe, n, thispe, pe, offset, local_size, dest_size
  logical(pgslib_log_type) :: local_error, abort

  npe = PGSLib_Inquire_nPE()
  thisPE = PGSLib_Inquire_thisPE_Actual()

  ! Set up for scalar test
  allocate(double_vec_given(nPE))
  double_vec_given = -1
  allocate(double_vec_expected(npe))
  double_vec_expected = (/ (n, n=npe, 1,-1) /)
  double_vec_expected = sin(double_vec_expected)

  double_src = npe - (thispe - 1)
  double_src = sin(double_src)

  call pgslib_collate(double_vec_given, double_src)

  local_error = .false.
  if (PGSLib_Inquire_thisPE() .eq. PGSLib_Inquire_IO_ROOT_PE()) then
     if (ANY(double_vec_given .ne. double_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed double scalar collate test")
        call pgslib_output("Failed double scalar collate test")
     endif
  endif

  deallocate(double_vec_expected)
  deallocate(double_vec_given)

  ! Set up for vector test
  local_size = thispe*2 + 1
  allocate(double_vec_src(local_size))
  double_vec_src = (/ (n, n=thispe*2+1,1,-1) /)
  double_vec_src = sin(double_vec_src)

  dest_size = SUM( (/ (pe*2+1, pe=1,npe) /) )
  allocate(double_vec_given(dest_size))
  allocate(double_vec_expected(dest_size))
  
  double_vec_given = -1
  
  offset = 1
  do pe=1,npe
     do n=pe*2+1,1,-1
        double_vec_expected(offset) = n
        offset = offset + 1
     enddo
  enddo
  double_vec_expected = sin(double_vec_expected)

  call pgslib_collate(double_vec_given, double_vec_src)

  if (PGSLib_Inquire_thisPE() .eq. PGSLib_Inquire_IO_ROOT_PE()) then
     if (ANY(double_vec_given .ne. double_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed double vector collate test")
        call pgslib_output("Failed double vector collate test")
     endif
  endif

  deallocate(double_vec_expected)
  deallocate(double_vec_given)
  deallocate(double_vec_src)  

  if (.not. local_error) then
     call pgslib_output("Passed collate double test")
  endif
  error = error .or. local_error

  return
end subroutine test_collate_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test collate of chars
subroutine test_collate_char(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  integer, parameter :: string_len = 611
  character (LEN=string_len) :: char_give, char_src
  character (LEN=string_len), pointer, dimension(:) :: char_vec_src, char_vec_given, char_vec_expected
  integer(pgslib_int_type):: npe, n, thispe, pe, offset, local_size, dest_size
  logical(pgslib_log_type) :: local_error, abort

  npe = PGSLib_Inquire_nPE()
  thisPE = PGSLib_Inquire_thisPE_Actual()

  ! Set up for scalar test
  allocate(char_vec_given(nPE))
  char_vec_given = 'X'
  allocate(char_vec_expected(npe))
  char_vec_expected = ACHAR( MOD( (/ (n, n=npe, 1,-1) /),128) )

  char_src = ACHAR( MOD( (npe - (thispe - 1)) ,128) )

  call pgslib_collate(char_vec_given, char_src)

  local_error = .false.
  if (PGSLib_Inquire_IO_P()) then
     if (ANY(char_vec_given .ne. char_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed char scalar collate test")
        call pgslib_output("Failed char scalar collate test")
     endif
  endif

  deallocate(char_vec_expected)
  deallocate(char_vec_given)

  ! Set up for vector test
  local_size = thispe*2 + 1
  allocate(char_vec_src(local_size))
  char_vec_src = 'X'
  char_vec_src = ACHAR( MOD( (/ (n, n=thispe*2+1,1,-1) /),128) )

  dest_size = SUM( (/ (pe*2+1, pe=1,npe) /) )
  allocate(char_vec_given(dest_size))
  allocate(char_vec_expected(dest_size))
  
  char_vec_given = 'Y'
  char_vec_expected = 'Z'
  
  offset = 1
  do pe=1,npe
     do n=pe*2+1,1,-1
        char_vec_expected(offset) = ACHAR( MOD(n,128) )
        offset = offset + 1
     enddo
  enddo

  call pgslib_collate(char_vec_given, char_vec_src)

  if (PGSLib_Inquire_IO_P()) then
     if (ANY(char_vec_given .ne. char_vec_expected) ) then
        local_error = .true.
        call pgslib_error("Failed char vector collate test")
        call pgslib_output("Failed char vector collate test")
     endif
  endif

  deallocate(char_vec_expected)
  deallocate(char_vec_given)
  deallocate(char_vec_src)  

  if (.not. local_error) then
     call pgslib_output("Passed collate char test")
  endif
  error = error .or. local_error

  return
end subroutine test_collate_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test collate of logs
subroutine test_collate_log(error)
  USE PGSLib_MODULE
  implicit none
  logical(pgslib_log_type) :: error

  logical (pgslib_log_type) :: log_given, log_src
  logical (pgslib_log_type), pointer, dimension(:):: log_vec_src, log_vec_given, log_vec_expected
  integer(pgslib_int_type):: npe, n, thispe, pe, offset, local_size, dest_size
  logical(pgslib_log_type) :: local_error, abort

  npe = PGSLib_Inquire_nPE()
  thisPE = PGSLib_Inquire_thisPE_Actual()

  ! Set up for scalar test
  allocate(log_vec_given(nPE))
  log_vec_given = .FALSE.
  allocate(log_vec_expected(npe))
  log_vec_expected = MOD((/ (n, n=npe, 1,-1) /),3) == 0

  log_src = MOD(npe - (thispe - 1), 3) == 0

  call pgslib_collate(log_vec_given, log_src)

  local_error = .false.
  if (PGSLib_Inquire_IO_P() ) then
     if (ANY((log_vec_given .or. log_vec_expected) .and. &
          &   .not.(log_vec_given .and. log_vec_expected) ) ) then ! == (log_g .xor. log_e)
        local_error = .true.
        call pgslib_error("Failed log scalar collate test")
        call pgslib_output("Failed log scalar collate test")
     endif
  endif

  local_error = PGSLib_Global_ANY(local_error)

  deallocate(log_vec_expected)
  deallocate(log_vec_given)

  ! Set up for vector test
  local_size = thispe*2 + 1
  allocate(log_vec_src(local_size))
  log_vec_src = MOD((/ (n, n=thispe*2+1,1,-1) /), 5) == 0

  dest_size = SUM( (/ (pe*2+1, pe=1,npe) /) )
  allocate(log_vec_given(dest_size))
  allocate(log_vec_expected(dest_size))
  
  log_vec_given = .FALSE.
  
  offset = 1
  do pe=1,npe
     do n=pe*2+1,1,-1
        log_vec_expected(offset) = MOD(n,5) == 0
        offset = offset + 1
     enddo
  enddo

  call pgslib_collate(log_vec_given, log_vec_src)

  if (PGSLib_Inquire_IO_P() ) then
  if ( ANY((log_vec_given .or. log_vec_expected) .and. &
       &    .not.(log_vec_given .and. log_vec_expected))  )then ! == (log_g .xor. log_e)     
        local_error = .true.
        call pgslib_error("Failed log vector collate test")
        call pgslib_output("Failed log vector collate test")
     endif
  endif

  local_error = PGSLib_Global_ANY(local_error)

  deallocate(log_vec_expected)
  deallocate(log_vec_given)
  deallocate(log_vec_src)  

  if (.not. local_error) then
     call pgslib_output("Passed collate log test")
  endif
  error = error .or. local_error

  return
end subroutine test_collate_log

