MODULE TEST_SCAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE -
  !   Test the scan routines
  !
  ! INPUT
  !   None
  ! OUTPUT
  !   Results of tests, warning_error and fatal_error
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  ! $Id: test_scan.F,v 1.1.1.1 2000/10/11 22:44:26 ferrell Exp $

  USE PGSLib_MODULE
  implicit none
  
  PRIVATE
  PUBLIC  :: Test_Scans

CONTAINS
  subroutine test_scans(local_error)
    implicit none
    logical :: local_error

    logical :: fatal_error, warning_error

    local_error = .false.

    fatal_error = .false.
    warning_error = .false.
    call test_scan_prefix_int_1D(fatal_error, warning_error)
    local_error = local_error .or. (fatal_error .or. warning_error)

    fatal_error = .false.
    warning_error = .false.
    call test_scan_prefix_single_1D(fatal_error, warning_error)
    local_error = local_error .or. (fatal_error .or. warning_error)

    fatal_error = .false.
    warning_error = .false.
    call test_scan_prefix_double_1D(fatal_error, warning_error)
    local_error = local_error .or. (fatal_error .or. warning_error)

    fatal_error = .false.
    warning_error = .false.
    call test_scan_suffix_int_1D(fatal_error, warning_error)
    local_error = local_error .or. (fatal_error .or. warning_error)

    fatal_error = .false.
    warning_error = .false.
    call test_scan_suffix_single_1D(fatal_error, warning_error)
    local_error = local_error .or. (fatal_error .or. warning_error)

    fatal_error = .false.
    warning_error = .false.
    call test_scan_suffix_double_1D(fatal_error, warning_error)
    local_error = local_error .or. (fatal_error .or. warning_error)

    return
  end subroutine test_scans

  subroutine test_scan_prefix_int_1D(fatal_error, warning_error)
#define _DATA_TYPE_ integer
#define _OP_ID_ 0
#define _DATA_TYPE_STRING_ " integer "
#define _RAND_OP_   
#define _START_ 1
#define _STOP_ SIZE(Dest_Expected_Tot, 1)
#define _STEP_  1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_PREFIX
#define _SCAN_TEST_STRING_    " PGSLib_SUM_PREFIX"

#include "test_scan.fpp"

  end subroutine test_scan_prefix_int_1D

  subroutine test_scan_prefix_single_1D(fatal_error, warning_error)
#define _DATA_TYPE_ REAL(PGSLib_REAL_TYPE)
#define _OP_ID_ 0.0_PGSLib_Real_Type
#define _DATA_TYPE_STRING_ " single "
#define _RAND_OP_   COS
#define _START_ 1
#define _STOP_  SIZE(Dest_Expected_Tot, 1)
#define _STEP_  1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_PREFIX
#define _SCAN_TEST_STRING_    " PGSLib_SUM_PREFIX"

#include "test_scan.fpp"

  end subroutine test_scan_prefix_single_1D

  subroutine test_scan_prefix_double_1D(fatal_error, warning_error)
#define _DATA_TYPE_ REAL(KIND(1.0D0))
#define _OP_ID_ 0.0D0
#define _DATA_TYPE_STRING_ ' double '
#define _RAND_OP_   COS
#define _START_ 1
#define _STOP_  SIZE(Dest_Expected_Tot, 1)
#define _STEP_  1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_PREFIX
#define _SCAN_TEST_STRING_    'PGSLib_SUM_PREFIX'

#include "test_scan.fpp"

  end subroutine test_scan_prefix_double_1D

  subroutine test_scan_suffix_int_1D(fatal_error, warning_error)
#define _DATA_TYPE_ integer
#define _OP_ID_ 0
#define _DATA_TYPE_STRING_ ' integer '
#define _RAND_OP_   
#define _START_ SIZE(Dest_Expected_Tot, 1)
#define _STOP_  1
#define _STEP_  -1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_SUFFIX
#define _SCAN_TEST_STRING_    'PGSLib_SUM_SUFFIX'

#include "test_scan.fpp"

  end subroutine test_scan_suffix_int_1D

  subroutine test_scan_suffix_single_1D(fatal_error, warning_error)
#define _DATA_TYPE_ REAL
#define _OP_ID_ 0.0
#define _DATA_TYPE_STRING_ ' single '
#define _RAND_OP_   COS
#define _START_ SIZE(Dest_Expected_Tot, 1) 
#define _STOP_  1
#define _STEP_  -1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_SUFFIX
#define _SCAN_TEST_STRING_    'PGSLib_SUM_SUFFIX'

#include "test_scan.fpp"

  end subroutine test_scan_suffix_single_1D

  subroutine test_scan_suffix_double_1D(fatal_error, warning_error)
#define _DATA_TYPE_ REAL(KIND(1.0D0))
#define _OP_ID_ 0.0D0
#define _DATA_TYPE_STRING_ ' double '
#define _RAND_OP_   COS
#define _START_ SIZE(Dest_Expected_Tot, 1)
#define _STOP_  1
#define _STEP_  -1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_SUFFIX
#define _SCAN_TEST_STRING_    'PGSLib_SUM_SUFFIX'

#include "test_scan.fpp"

  end subroutine test_scan_suffix_double_1D

END MODULE TEST_SCAN
