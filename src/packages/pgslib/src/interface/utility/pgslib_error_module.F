MODULE PGSLib_ERROR_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE
  !    Error routines to support PGSLib.  These are in a distinct module
  !    so they can get compiled before most other modules.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_error_module.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $
  USE PGSLib_Type_MODULE

  IMPLICIT NONE
  SAVE
  PUBLIC:: PGSLib_Error, PGSLib_Abort, PGSLib_Fatal_Error

  INTERFACE PGSLib_Error
     MODULE PROCEDURE PGSLib_Error_F
  END INTERFACE!PGSLib_Error_F

  INTERFACE PGSLib_Abort
     MODULE PROCEDURE PGSLib_Abort_F
  END INTERFACE!PGSLib_Abort

  INTERFACE PGSLib_Fatal_Error
     MODULE PROCEDURE PGSLib_Fatal_Error_F
  END INTERFACE!PGSLib_Fatal_Error

  INTERFACE PGSLib_Error_C
     SUBROUTINE PGSLib_Error_C(EString)
       IMPLICIT NONE
       CHARACTER (*):: Estring
     END SUBROUTINE PGSLib_Error_C

  END INTERFACE! PGSLib_Error_C

  INTERFACE PGSLib_Abort_C
     SUBROUTINE PGSLib_Abort_C()
     END SUBROUTINE PGSLib_Abort_C

  END INTERFACE! PGSLib_Abort_C

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Error_F(IString, ROUTINE)
    IMPLICIT NONE
    CHARACTER (LEN=*):: IString
    character (LEN=*), intent(in), optional :: ROUTINE

    ! Local variables
    INTEGER :: lf1, lf
    CHARACTER (LEN=1024):: Estring

    EString = ''
    ! IF ROUTINE is present, prepend it to output
    if (PRESENT(ROUTINE)) then
       lf1 = MIN(len_trim(ROUTINE), 1023)
       EString(1:lf1) = ROUTINE(1:lf1)
    ELSE
       lf1 = 0
    ENDIF
    ! Append a NULL to the end of the string, then call a C routine
    !      to do the output.
    lf = MIN(len_trim(IString),1023-lf1)
    EString(lf1+1:lf1+lf) = IString(1:lf)
    EString(lf1+lf+1:) = ACHAR(0)
    Call PGSLib_Error_C(Estring)
    RETURN
  END SUBROUTINE PGSLib_Error_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Abort_F()
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    Call PGSLib_Abort_C()
    RETURN
  END SUBROUTINE PGSLib_Abort_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Fatal_Error_F(IString)
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    CHARACTER (LEN=*):: IString

    Call PGSLib_Error(IString)
    Call PGSLib_Abort()

    RETURN
  END SUBROUTINE PGSLib_Fatal_Error_F

END MODULE PGSLib_ERROR_MODULE
