!!CPP!! This file is included in the reduction routines.
!!CPP!! _DATA_TYPE_, _ARG_, _ROUTINE_ must be defined before including this file

! $Id: red_sum.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _ARRAY_
#error "_ARRAY_ must be defined before including this file"
#endif

#ifndef _ARG_
#error "_ARG_ must be defined before including this file"
#endif

#ifndef _ROUTINE_
#error "_ROUTINE_ must be defined before including this file"
#endif

!!CPP!! Need this strange set of macros to get proper expansion.
#define _STRINGIZE_(S)           #S
#define _STRING_(S)   _STRINGIZE_(S)

function _ROUTINE_(A, MASK, SCOPE)
  USE PGSLib_Type_MODULE
  use pgslib_globals_module
  use pgslib_stats,         only: GLOBAL_SUM_STATISTICS, Enter_Routine, Exit_Routine
  IMPLICIT NONE
  _DATA_TYPE_  :: _ROUTINE_
  _DATA_TYPE_ ,              intent(IN   ) _ARRAY_ :: A
  logical (PGSLib_Log_Type), intent(IN   ) _ARRAY_ , OPTIONAL :: MASK
  type (PGSLib_SCOPE),       intent(IN   ),          OPTIONAL :: SCOPE

  ! Local temporaries
  _DATA_TYPE_ :: SumA
  logical     :: Global

#ifdef USE_TIMERS_2
  call Enter_Routine(GLOBAL_SUM_STATISTICS())
#endif

  ! Local or Global operation
#include "red_global_test.fpp"

  ! On-PE reduction
  IF (PRESENT(MASK)) THEN
     SumA = SUM ( _ARG_ , MASK=MASK)
  ELSE
     SumA = SUM ( _ARG_ )
  END IF

  ! Across-PE reduction if global
  IF (Global) Call PGSLib_global_sum_c(SumA)
  _ROUTINE_ = SumA

#ifdef USE_TIMERS_2
  call Exit_Routine(GLOBAL_SUM_STATISTICS())
#endif

  RETURN
END Function _ROUTINE_

#undef _DATA_TYPE_
#undef _ARRAY_
#undef _ARG_
#undef _ROUTINE_
#undef _STRING_
#undef _STRINGIZE_
