!!CPP!! This file is included in the reduction routines.
!!CPP!! This file contains the code fragment for determining the 
!!CPP!! SCOPE of the operation.
!!CPP!! _GENERIC_ROUTINE_NAME_ must be defined before including this file

! $Id: red_global_test.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

#ifndef _GENERIC_ROUTINE_NAME_
#error "_GENERIC_ROUTINE_NAME_ must be defined before including this file"
#endif

!!CPP!! Global is a logical defined in the host.  It is set to 
!!CPP!! .TRUE. or .FALSE. by this code fragment

  ! Local or Global operation
  if (PRESENT(SCOPE)) then
     OP_SCOPE: SELECT CASE(SCOPE%SCOPE)
     CASE(PGSLib_Local%SCOPE)
        Global = .FALSE.
     CASE(PGSLib_Global%SCOPE)
        Global = .TRUE.
     CASE DEFAULT
        call pgslib_error("SCOPE argument must be one of PGSLib_Local or PGSLib_Global", &
             &            _STRING_(_GENERIC_ROUTINE_NAME_))
     END SELECT OP_SCOPE
  ELSE
     Global = .TRUE.
  END if
