MODULE PGSLib_GS_Util_MODULE
  use PGSlib_Stats
  use pgslib_timing_module

  implicit none 

  save
  private
  public :: PGSLib_GS_Init_Trace
  public :: PGSLib_GS_Release_Trace
  public :: PGSLib_Deallocate_Trace
  public :: PGSLib_GS_Trace_Setup_P

  INTERFACE PGSLib_GS_Init_Trace
     MODULE PROCEDURE PGSLib_GS_Init_Trace_F
  END INTERFACE

  INTERFACE PGSLib_GS_Release_Trace
     MODULE PROCEDURE PGSLib_GS_Release_Trace_F
  END INTERFACE

  ! PGSLib_Deallocate_Trace is just another name for PGSLib_GS_Release_Trace
  INTERFACE PGSLib_Deallocate_Trace
     MODULE PROCEDURE PGSLib_GS_Release_Trace_F
  END INTERFACE
  
     

! Interfaces to the C routines
  INTERFACE 
     subroutine PGSLib_GS_Init_Trace_C(GS_Trace)
       implicit none
       integer, DIMENSION(*) :: GS_Trace
     END subroutine PGSLib_GS_Init_Trace_C
  END INTERFACE
  INTERFACE 
     subroutine PGSLib_GS_Release_Trace_C(GS_Trace)
       implicit none
       integer, DIMENSION(*) :: GS_Trace
     END subroutine PGSLib_GS_Release_Trace_C
  END INTERFACE
       

CONTAINS
  
  !=====================================================================
  !          PGSLib_GS_Init_Trace_F
  ! PURPOSE: 
  !          Initialize a trace.  This sets the components of the F90 structure
  !          and it also calls a C routine to set the C stuff.
  !=====================================================================

  function PGSLib_GS_Init_Trace_F() RESULT(Trace)
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE
    USE PGSLib_Index_Module
    implicit none
    type (PGSLib_GS_Trace), POINTER :: Trace


#ifdef USE_TIMERS_1
      call Enter_Routine(GS_INIT_TRACE_STATISTICS())
#endif

    ! ALLOCATE the trace
    ALLOCATE(Trace)
      
    ! Set various F90 visible things to default values.
    Trace%SetupP = .FALSE.
    Trace%N_Supplement = -1
    Trace%N_Duplicate  = -1

    NULLIFY(Trace%Supplement_Global_Indices)
    NULLIFY(Trace%Supplement_Global_PEs)
    NULLIFY(Trace%Supplement_I_Data)
    NULLIFY(Trace%Supplement_R_Data)
    NULLIFY(Trace%Supplement_D_Data)
    NULLIFY(Trace%Supplement_L_Data)

    NULLIFY(Trace%Supplement_2_I_Data)
    NULLIFY(Trace%Supplement_2_R_Data)
    NULLIFY(Trace%Supplement_2_D_Data)
    NULLIFY(Trace%Supplement_2_L_Data)

    NULLIFY(Trace%Duplicate_Indices)
    NULLIFY(Trace%Duplicate_I_Data)
    NULLIFY(Trace%Duplicate_R_Data)
    NULLIFY(Trace%Duplicate_D_Data)
    NULLIFY(Trace%Duplicate_L_Data)

    NULLIFY(Trace%Duplicate_2_I_Data)
    NULLIFY(Trace%Duplicate_2_R_Data)
    NULLIFY(Trace%Duplicate_2_D_Data)
    NULLIFY(Trace%Duplicate_2_L_Data)

    NULLIFY(Trace%Global_Index_1)

    call pgslib_reset_timer(Trace%GatherTotal_Timer)
    call pgslib_reset_timer(Trace%ScatterTotal_Timer)
    call pgslib_reset_timer(Trace%GatherBuffer_Timer)
    call pgslib_reset_timer(Trace%ScatterBuffer_Timer)
    call pgslib_reset_timer(Trace%Setup_Timer)

#ifdef USE_TIMERS_1
      call PGSLib_Start_Timer(Trace%Setup_Timer)
#endif

    ! Allocate enough space to hold a pointer
    ALLOCATE(Trace%GS_Trace(4))
    Call PGSLib_GS_Init_Trace_C(Trace%GS_Trace)

!    Trace%Access_Table => PGSLib_Init_offPE_Access_Table()

#ifdef USE_TIMERS_1
      call Exit_Routine(GS_INIT_TRACE_STATISTICS())
      call PGSLib_Stop_Timer(Trace%Setup_Timer)
#endif


    RETURN
  end function PGSLib_GS_Init_Trace_F

  !======================================================================
  !          PGSLib_GS_Release_Trace_F
  ! PURPOSE:
  !          Release the storage used by trace and reset the trace.
  !
  !======================================================================

  subroutine PGSLib_GS_Release_Trace_F(Trace)

    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE
    Use PGSLib_Index_MODULE
    implicit none
    type (PGSLib_GS_Trace), POINTER :: Trace


#ifdef USE_TIMERS_1
    call Enter_Routine(GS_RELEASE_TRACE_STATISTICS())
    call PGSLib_Start_Timer(Trace%Setup_Timer)
#endif


    ! Release stuff internal to C structures.
    call pgslib_gs_release_trace_c(Trace%GS_Trace)

    ! Release the access table stuff
    ! call pgslib_free_offPE_access_table(Trace%Access_Table)

    ! Release each component
    
    if (ASSOCIATED(Trace%Supplement_Global_Indices)) DEALLOCATE(Trace%Supplement_Global_Indices)
    if (ASSOCIATED(Trace%Supplement_Global_PEs))     DEALLOCATE(Trace%Supplement_Global_PEs)
    if (ASSOCIATED(Trace%Supplement_I_Data))         DEALLOCATE(Trace%Supplement_I_Data)
    if (ASSOCIATED(Trace%Supplement_R_Data))         DEALLOCATE(Trace%Supplement_R_Data)
    if (ASSOCIATED(Trace%Supplement_D_Data))         DEALLOCATE(Trace%Supplement_D_Data)
    if (ASSOCIATED(Trace%Supplement_L_Data))         DEALLOCATE(Trace%Supplement_L_Data)

    if (ASSOCIATED(Trace%Supplement_2_I_Data))       DEALLOCATE(Trace%Supplement_2_I_Data)
    if (ASSOCIATED(Trace%Supplement_2_R_Data))       DEALLOCATE(Trace%Supplement_2_R_Data)
    if (ASSOCIATED(Trace%Supplement_2_D_Data))       DEALLOCATE(Trace%Supplement_2_D_Data)
    if (ASSOCIATED(Trace%Supplement_2_L_Data))       DEALLOCATE(Trace%Supplement_2_L_Data)

    if (ASSOCIATED(Trace%Duplicate_Indices))         DEALLOCATE(Trace%Duplicate_Indices)
    if (ASSOCIATED(Trace%Duplicate_I_Data))          DEALLOCATE(Trace%Duplicate_I_Data)
    if (ASSOCIATED(Trace%Duplicate_R_Data))          DEALLOCATE(Trace%Duplicate_R_Data)
    if (ASSOCIATED(Trace%Duplicate_D_Data))          DEALLOCATE(Trace%Duplicate_D_Data)
    if (ASSOCIATED(Trace%Duplicate_L_Data))          DEALLOCATE(Trace%Duplicate_L_Data)

    if (ASSOCIATED(Trace%Duplicate_2_I_Data))        DEALLOCATE(Trace%Duplicate_2_I_Data)
    if (ASSOCIATED(Trace%Duplicate_2_R_Data))        DEALLOCATE(Trace%Duplicate_2_R_Data)
    if (ASSOCIATED(Trace%Duplicate_2_D_Data))        DEALLOCATE(Trace%Duplicate_2_D_Data)
    if (ASSOCIATED(Trace%Duplicate_2_L_Data))        DEALLOCATE(Trace%Duplicate_2_L_Data)

    if (ASSOCIATED(Trace%Global_Index_1))            DEALLOCATE(Trace%Global_Index_1)

    if (ASSOCIATED(Trace%GS_Trace))                  DEALLOCATE(Trace%GS_Trace)

#ifdef USE_TIMERS_1
    call PGSLib_Stop_Timer(Trace%Setup_Timer)
#endif
    DEALLOCATE(Trace)
#ifdef USE_TIMERS_1
    call Exit_Routine(GS_RELEASE_TRACE_STATISTICS())
#endif

    RETURN
  END subroutine PGSLib_GS_Release_Trace_F

  function PGSLib_GS_Trace_Setup_P(Trace)
    ! Returns .TRUE. if Trace has been setup
    ! Returns .FALSE. otherwise.
    use PGSLib_Type_MODULE
    implicit none
    
    type (PGSLib_GS_Trace)   :: Trace
    logical (PGSLib_Log_Type) :: PGSLib_GS_Trace_Setup_P

    PGSLib_GS_Trace_Setup_P = Trace%SetupP
  end function PGSLib_GS_Trace_Setup_P
  
END MODULE PGSLib_GS_Util_MODULE
