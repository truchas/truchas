!!CPP!! This file is included in the reduction routines.
!!CPP!! _MASK_, _ROUTINE_  must be defined before including this file

! $Id: red_any.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

#ifndef _MASK_
#error "_MASK_ must be defined before including this file"
#endif

#ifndef _MASK_SHAPE_
#error "_MASK_SHAPE_ must be defined before including this file"
#endif

#ifndef _ROUTINE_
#error "_ROUTINE_ must be defined before including this file"
#endif

function _ROUTINE_ (MASK, SCOPE)
  USE pgslib_globals_module
  use pgslib_stats,         only: GLOBAL_ANY_STATISTICS, Enter_Routine, Exit_Routine
  USE PGSLib_Type_MODULE
  IMPLICIT NONE
  logical (PGSLib_LOG_TYPE) _ROUTINE_
  logical (PGSLib_Log_Type),                 &
       intent(IN   )                         &
       _MASK_SHAPE_                          :: MASK
  type (PGSLib_SCOPE),                       &
       intent(in   ),                        &
       optional                              :: SCOPE

  ! Local temporaries
  logical     :: Global

  ! Notice that this must be integer type to pass through to C properly
  integer (PGSLib_Int_TYPE) AllMask

#ifdef USE_TIMERS_2
  call Enter_Routine(GLOBAL_ANY_STATISTICS())
#endif

  ! Local or Global operation
#include "red_global_test.fpp"
  
  ! On-PE reduction
  AllMask = MERGE(PGSLib_TRUE, PGSLib_FALSE, ANY( _MASK_ ) )

  ! Global reduction, with broadcast of result
  If (Global) Call PGSLib_global_any_C(AllMask)

  ! C routine returns integer, need to convert to logical
  ! Anything .NOT. PGSLib_FALSE is .true.
  _ROUTINE_ = .NOT.(ALLMask == PGSLib_FALSE) 

#ifdef USE_TIMERS_2
  call Exit_Routine(GLOBAL_ANY_STATISTICS())
#endif

  RETURN
END FUNCTION _ROUTINE_

#undef _MASK_
#undef _MASK_SHAPE_
#undef _ROUTINE_
