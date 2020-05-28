!!CPP!! This file is included in the reduction routines.
!!CPP!! _ARG_, _ROUTINE_  must be defined before including this file

! $Id: red_max.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

function _ROUTINE_(A, MASK, SCOPE)
  IMPLICIT NONE
  _DATA_TYPE_ :: _ROUTINE_
  _DATA_TYPE_,               intent(IN   ) _ARRAY_           :: A
  logical (PGSLib_Log_Type), intent(IN   ) _ARRAY_, OPTIONAL :: MASK
  type (PGSLib_SCOPE),       intent(IN   ),          OPTIONAL :: SCOPE

  ! Local temporaries
  _DATA_TYPE_ MaxA
  logical     :: Global

#ifdef USE_TIMERS_2
  call Enter_Routine(GLOBAL_MAX_STATISTICS())
#endif

  ! Local or Global operation
#include "red_global_test.fpp"


  ! On-PE reduction
  IF (PRESENT(MASK)) THEN
     MAXA = MAXVAL ( _ARG_ , MASK=MASK)
  ELSE
     MAXA = MAXVAL ( _ARG_ )
  END IF

  ! Across-PE reduction
  If (Global) Call PGSLib_global_MAX_c(MAXA)
  _ROUTINE_ = MAXA

#ifdef USE_TIMERS_2
  call Exit_Routine(GLOBAL_MAX_STATISTICS())
#endif

  RETURN
END FUNCTION _ROUTINE_

#undef _DATA_TYPE_
#undef _ARRAY_
#undef _ARG_
#undef _ROUTINE_
