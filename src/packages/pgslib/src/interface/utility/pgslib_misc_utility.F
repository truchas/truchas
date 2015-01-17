MODULE PGSLib_Misc_Utility
  USE PGSLib_Type_MODULE
  use pgslib_globals_module
  use pgslib_error_module
  use pgslib_c_binding

  IMPLICIT NONE
  SAVE
  PUBLIC:: PGSLib_Inquire_nPE, PGSLib_Inquire_thisPE, PGSLib_Inquire_IO_ROOT_PE
  PUBLIC:: PGSLib_Inquire_thisPE_Actual, PGSLib_Inquire_IO_P
  PUBLIC:: PGSLib_UseGlobalServices
  PUBLIC:: PGSLib_check_error
  PUBLIC:: PGSLib_Output, PGSLib_Flush_Output
  PUBLIC:: PGSLib_Barrier
  PUBLIC:: PGSLib_Scope_Check


  ! $Id: pgslib_misc_utility.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $

  ! Variables for use by source which USEs this module
  CHARACTER (LEN=1024), PUBLIC :: pgslib_out_string
  integer (PGSLib_Int_TYPE), PRIVATE, PARAMETER :: FILE_NAME_LENGTH = 2048

  INTERFACE PGSLib_Inquire_nPE
     MODULE PROCEDURE PGSLib_Inquire_nPE_F
  END INTERFACE!PGSLib_Inquire_nPE

  INTERFACE PGSLib_Inquire_thisPE
     MODULE PROCEDURE PGSLib_Inquire_thisPE_F
  END INTERFACE!PGSLib_Inquire_thisPE

  INTERFACE PGSLib_Inquire_thisPE_Actual
     MODULE PROCEDURE PGSLib_Inquire_thisPE_Actual_F
  END INTERFACE!PGSLib_Inquire_thisPE

  INTERFACE PGSLib_Inquire_IO_ROOT_PE
     MODULE PROCEDURE PGSLib_Inquire_IO_ROOT_PE_F
  END INTERFACE!PGSLib_Inquire_IO_ROOT_PE

  INTERFACE PGSLib_Inquire_IO_P
     MODULE PROCEDURE PGSLib_Inquire_IO_P_F
  END INTERFACE!PGSLib_Inquire_IO_P

  INTERFACE PGSLib_UseGlobalServices
     MODULE PROCEDURE PGSLib_UseGlobalServices_F
  END INTERFACE

  INTERFACE PGSLib_Output
     MODULE PROCEDURE PGSLib_Output_F
     MODULE PROCEDURE PGSLib_Output_Array_F
  END INTERFACE!PGSLib_Output_F

  INTERFACE PGSLib_Flush_Output
     MODULE PROCEDURE PGSLib_Flush_Output_F
  END INTERFACE!PGSLib_Flush_Output_F

  INTERFACE PGSLib_Check_Error
     MODULE PROCEDURE PGSLib_Check_Error_F
  END INTERFACE!PGSLib_Check_Error

  INTERFACE PGSLib_Barrier
     MODULE PROCEDURE PGSLib_Barrier_F
  END INTERFACE!PGSLib_Barrier

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function PGSLib_Inquire_nPE_F( )
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE, ONLY : PGSLib_PEInfo
    IMPLICIT NONE
    INTEGER (PGSLib_Int_Type) PGSLib_Inquire_nPE_F

    PGSLib_Inquire_nPE_F = PGSLib_PEInfo%nPE

    RETURN
  END FUNCTION PGSLib_Inquire_nPE_F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Inquire_thisPE_F()
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE, ONLY : PGSLib_PEInfo
    IMPLICIT NONE
    INTEGER (PGSLib_Int_Type) PGSLib_Inquire_thisPE_F

    IF (PGSLib_PEInfo%IO_ROOT_PE < 0) THEN
       PGSLib_Inquire_thisPE_F = PGSLib_Inquire_IO_ROOT_PE()
    ELSE
       PGSLib_Inquire_thisPE_F = PGSLib_PEInfo%thisPE
    ENDIF

    RETURN
  END FUNCTION PGSLib_Inquire_thisPE_F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Inquire_thisPE_Actual_F()
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE, ONLY : PGSLib_PEInfo
    IMPLICIT NONE
    INTEGER (PGSLib_Int_Type) PGSLib_Inquire_thisPE_Actual_F

    PGSLib_Inquire_thisPE_Actual_F = PGSLib_PEInfo%thisPE

    RETURN
  END FUNCTION PGSLib_Inquire_thisPE_Actual_F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Inquire_IO_ROOT_PE_F()
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE, ONLY : PGSLib_PEInfo
    IMPLICIT NONE
    INTEGER (PGSLib_Int_Type) PGSLib_Inquire_IO_ROOT_PE_F

    PGSLib_Inquire_IO_ROOT_PE_F = PGSLib_PEInfo%IO_ROOT_PE

    RETURN
  END FUNCTION PGSLib_Inquire_IO_ROOT_PE_F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Inquire_IO_P_F()
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE, ONLY : PGSLib_PEInfo
    IMPLICIT NONE
    logical (PGSLib_Log_Type) :: PGSLib_Inquire_IO_P_F

    PGSLib_Inquire_IO_P_F = (PGSLib_PEInfo%IO_ROOT_PE < 0) .OR. &
         &                  (PGSLib_PEInfo%IO_ROOT_PE == PGSLib_PEInfo%thisPE)

    RETURN
    END FUNCTION PGSLib_Inquire_IO_P_F
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function PGSLib_UseGlobalServices_F()
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE, ONLY : PGSLib_UseGlobalServicesFlag
    IMPLICIT NONE 
    INTEGER (PGSLib_Int_Type) PGSLib_UseGlobalServices_F

    PGSLib_UseGlobalServices_F = PGSLib_UseGlobalServicesFlag

    RETURN
  END FUNCTION PGSLib_UseGlobalServices_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Output_F(IString)
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    CHARACTER (LEN=*):: IString

    ! Local variables
    INTEGER lf
    CHARACTER (LEN=1024):: Estring

    ! Append a NULL to the end of the string, then call a C routine
    !      to do the output.
    lf = MIN(len_trim(IString),1023)
    EString(1:lf) = IString(1:lf)
    EString(lf+1:) = ACHAR(0)
    Call PGSLib_Output_C(Estring)
    RETURN
  END SUBROUTINE PGSLib_Output_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Output_Array_F(IString_Array)
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    CHARACTER (LEN=*), DIMENSION(:):: IString_Array

    ! Local variables
    INTEGER n
    do n = 1, SIZE(IString_Array)
       ! If string is empty, do not print
       if (LEN_TRIM(ISTring_Array(n)) <= 0) cycle
       call pgslib_output(Istring_Array(n))
    end do

    RETURN
  END SUBROUTINE PGSLib_Output_Array_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Flush_Output_F()
    IMPLICIT NONE

    Call PGSLib_Flush_Output_C()
    RETURN
  END SUBROUTINE PGSLib_Flush_Output_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Check_Error_F(ErrorFlag, IString)
    USE PGSLib_Type_MODULE
    USE PGSLib_Red_Numeric_MODULE,  ONLY : PGSLib_Global_ANY
    IMPLICIT NONE
    logical (pgslib_log_type):: ErrorFlag
    CHARACTER (LEN=*):: IString

    IF (pgslib_global_any(ErrorFlag)) Then
       call PGSLib_Fatal_Error(IString)
    endif

    RETURN
  END SUBROUTINE PGSLib_Check_Error_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Barrier_F()
    implicit none

    call PGSLib_Barrier_C()

    RETURN
  END SUBROUTINE PGSLib_Barrier_F
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Scope_Check(Scope)
    implicit none
    type (PGSLib_SCOPE), OPTIONAL :: Scope
    type (PGSLib_Scope)           :: PGSLib_Scope_Check
    
    ! Local or Global operation
    if (PRESENT(SCOPE)) then
       OP_SCOPE: SELECT CASE(SCOPE%SCOPE)
       CASE(PGSLib_Local%SCOPE)
          PGSLib_Scope_Check = PGSLib_Local
       CASE(PGSLib_Global%SCOPE)
          PGSLib_Scope_Check = PGSLib_Global
       CASE DEFAULT
          call pgslib_error("SCOPE argument must be one of PGSLib_Local or PGSLib_Global")
       END SELECT OP_SCOPE
    ELSE
       PGSLib_Scope_Check = PGSLib_Global
    END if
  
  RETURN
  END function PGSLib_Scope_Check
  


END MODULE PGSLib_Misc_Utility

