!!CPP!! This file is included in the reduction routines.
!!CPP!! _ROUTINE_, _DATA_TYPE_ must be defined before including this file

! $Id: red_dtp.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _ROUTINE_
#error "_ROUTINE_ must be defined before including this file"
#endif

!!CPP!! Need this strange set of macros to get proper expansion.
#define _STRINGIZE_(S)           #S
#define _STRING_(S)   _STRINGIZE_(S)

function _ROUTINE_(Vector_A, Vector_B, SCOPE)
  USE PGSLib_Type_MODULE
  use pgslib_globals_module
  use pgslib_stats,         only: GLOBAL_DOT_PRODUCT_STATISTICS, Enter_Routine, Exit_Routine
  IMPLICIT NONE
  _DATA_TYPE_ _ROUTINE_
  _DATA_TYPE_, INTENT(in   ),                &
               DIMENSION(:)   :: Vector_A
  _DATA_TYPE_, INTENT(in   ),                &
               DIMENSION(:)   :: Vector_B
  type (PGSLib_SCOPE),                       &
               intent(in   ),                &
               optional       :: SCOPE

  ! Local temporaries
  _DATA_TYPE_ :: AdotB
  logical     :: Global

#ifdef USE_TIMERS_2
  call Enter_Routine(GLOBAL_DOT_PRODUCT_STATISTICS())
#endif

  ! Local or Global operation
#include "red_global_test.fpp"
     
  ! On-PE reduction
  ADotB = DOT_PRODUCT(Vector_A, Vector_B)

  ! Off-PE reduction, if requested
  IF (Global) then
     _ROUTINE_ = PGSLib_Global_SUM(AdotB)
  ELSE
     _ROUTINE_ = ADotB
  END IF
  
#ifdef USE_TIMERS_2
  call Exit_Routine(GLOBAL_DOT_PRODUCT_STATISTICS())
#endif

  RETURN

END FUNCTION _ROUTINE_

#undef _ROUTINE_
#undef _DATA_TYPE_
#undef _STRINGIZE_
#undef _STRING_
