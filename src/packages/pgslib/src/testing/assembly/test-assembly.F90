PROGRAM Test_Assembly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE
  !   Test routines used for matrix assembly.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: test-assembly.F,v 1.1.1.1 2000/10/11 22:44:26 ferrell Exp $

  USE PGSLib_MODULE
  USE TEST_SCAN
  USE Test_Sort
  implicit none

  logical :: error_this_test, error_detected
  character (LEN=256):: out_string


  call pgslib_initialize(IO_PE=1, FILE_PREFIX='ASSEMBLY_Test')

!!$  if (PGSLib_Inquire_IO_P()) then
!!$     print *, '<Enter> to continue'
!!$     read (*,*)
!!$  end if
  call pgslib_barrier()

  error_detected = .false.

  error_this_test = .false.
  call test_scans(error_this_test)
  error_detected = error_detected .or. error_this_test

  error_this_test = .false.
  call pgslib_output("Starting sort tests")
  call test_sorting(error_this_test)
  error_detected = error_detected .or. error_this_test

  if (error_detected) then
     OUT_STRING = ''
     WRITE(OUT_STRING, *) 'FAILED: Error detected in some assembly tests'
     call pgslib_output(OUT_STRING)
     call pgslib_error(out_string)
  else
     out_string = ''
     WRITE(OUT_STRING, *) 'SUCCESS: All assembly tests passed'
     call pgslib_output(OUT_STRING)
  end if

  call pgslib_finalize()

  STOP
  END
