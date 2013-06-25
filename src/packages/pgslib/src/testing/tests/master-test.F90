! This is the master test program for PGSLib.  It is used to test both the serial
! and the parallel versions of the library.

PROGRAM Master_Test
  
  USE PGSLib_MODULE
  ! The value input for the IO_ROOT_PE indicates whether the library is to support 
  ! parallel I/O or not.

  implicit none

  integer, parameter :: IO_ROOT_PE_SUGGESTED = 1
  integer :: IO_ROOT_PE
  
  logical(pgslib_log_type) :: error, fatal_error, error_detected, error_this_test

  character (LEN=256) :: OutString

  IO_ROOT_PE = IO_ROOT_PE_SUGGESTED
  call pgslib_initialize(IO_ROOT_PE, FILE_PER_PE=.FALSE.)

  OutString = ''
  call pgslib_output(OutString)
  WRITE(OutString,*) 'Starting Master Tests'
  call pgslib_output(OutString)

  WRITE(OutString,10) PGSLib_Inquire_nPE(), PGSLib_Inquire_thisPE(), PGSLib_Inquire_IO_ROOT_PE()
  call pgslib_output(OutString)
10 FORMAT(1x,'nPE = ', I6, ' thisPE = ', I6, ' IO_ROOT_PE = ', I6)


  IF (PGSLib_Inquire_IO_ROOT_PE() .NE. IO_ROOT_PE) THEN
     OutString = ''
     WRITE(OutString,*) 'Initialize changed IO_ROOT_PE to ', PGSLib_Inquire_IO_ROOT_PE()
     call pgslib_output(OutString)
  ENDIF

  IF (PGSLib_Inquire_thisPE() .NE. PGSLib_Inquire_thisPE_Actual() ) THEN
     OutString = ''
     WRITE(OutString,*) 'User thinks this PE = ', PGSLib_Inquire_thisPE()
     call pgslib_output(OutString)
     OutString = ''
     WRITE(OutString,*) '  Actually, this PE = ', PGSLib_Inquire_thisPE_Actual()
     call pgslib_output(OutString)
  ENDIF


  error_detected = .false.
  call test_redux_prelim()
  
  WRITE(Outstring, *) ' At start of main tests, VMSize = ', PGSLib_Memory_Size()
  call pgslib_output(Outstring)

  error_this_test = .false.
  call test_bcast(error_this_test)
  error_detected = error_detected .OR. error_this_test

  error_this_test = .false.
  call test_dist(error_this_test)
  error_detected = error_detected .OR. error_this_test
  
  error_this_test = .false.
  call test_collate(error_this_test)
  error_detected = error_detected .OR. error_this_test
  

  error_this_test = .false.
  call test_redux(error_this_test)
  error_detected = error_detected .OR. error_this_test
  
  WRITE(Outstring, *) ' At completion of main tests, VMSize = ', PGSLib_Memory_Size()
  call pgslib_output(Outstring)


  WRITE(Outstring,*) 'Time used for global dot_product = ', &
       &                 PGSLib_Read_Maximum_Time(GLOBAL_DOT_PRODUCT_STATISTICS())
  call pgslib_output(OutString)

  WRITE(Outstring,*) 'Time used for global sum = ', &
       &                 PGSLib_Read_Maximum_Time(GLOBAL_SUM_STATISTICS())
  call pgslib_output(OutString)

  if (error_detected) then
     Outstring = ''
     WRITE(Outstring,*) 'FAILED: Error detected in some tests'
     call pgslib_output(Outstring)
     call PGSLib_Error(Outstring)
  else
     Outstring = ''
     WRITE(Outstring, *) 'SUCCESS: All tests passed'
     call pgslib_output(Outstring)
  endif

  call pgslib_finalize()

  STOP
  END
