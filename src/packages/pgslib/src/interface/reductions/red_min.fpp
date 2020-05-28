!!CPP!! This file is included in the reduction routines.
!!CPP!! _ARG_, _ROUTINE_  must be defined before including this file

! $Id: red_min.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

function _ROUTINE_(A, MASK, SCOPE)
  IMPLICIT NONE
  _DATA_TYPE_ :: _ROUTINE_
  _DATA_TYPE_,               intent(IN   ) _ARRAY_           :: A
  logical (PGSLib_Log_Type), intent(IN   ) _ARRAY_, OPTIONAL :: MASK
  type (PGSLib_SCOPE),       intent(IN   ),          OPTIONAL :: SCOPE

  ! Local temporaries
  _DATA_TYPE_ MinA
  logical     :: Global

#ifdef USE_TIMERS_2
  call Enter_Routine(GLOBAL_MIN_STATISTICS())
#endif

  ! Local or Global operation
#include "red_global_test.fpp"

  ! On-PE reduction
  IF (PRESENT(MASK)) THEN
     MINA = MINVAL ( _ARG_ , MASK=MASK)
  ELSE
     MINA = MINVAL ( _ARG_ )
  END IF

  ! Across-PE reduction
  If (Global) Call PGSLib_global_MIN_c(MINA)
  _ROUTINE_ = MINA

#ifdef USE_TIMERS_2
  call Exit_Routine(GLOBAL_MIN_STATISTICS())
#endif

  RETURN
END FUNCTION _ROUTINE_

#undef _DATA_TYPE_
#undef _ARRAY_
#undef _ARG_
#undef _ROUTINE_
