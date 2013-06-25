!!CPP!! This file is included by pgslib_permute_module.F
!!CPP!! This file provides the redistribute functionality
!!CPP!! which is a special case of permutation

! $Id: redist.fpp,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file."
#endif

#ifndef _SCATTER_OP_
#error "_SCATTER_OP_ must be defined before including this file"
#endif

#ifndef _OP_ID_
#error "_OP_ID_ must be defined before including this file"
#endif

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif


  Subroutine _ROUTINE_NAME_(Dest, Source, TRACE) 
    USE PGSLib_Type_Module,      ONLY : PGSLib_Int_Type,   &
                                        PGSLib_Real_Type,  &
                                        PGSLib_Double_Type,  &
                                        PGSLib_Log_Type,   &
                                        PGSLib_GS_Trace
    
    USE PGSLib_GS_MODULE,        ONLY : PGSLib_Setup_Trace,      &
         &                              PGSLib_Deallocate_Trace
    USE PGSLIB_User_GS_MODULE,   ONLY : PGSLib_Scatter_SUM, PGSLib_Scatter_OR
    USE PGSLib_Reductions_MODULE,ONLY : PGSLib_Global_SUM
    USE PGSLib_Utility_MODULE,   ONLY : pgslib_error
    USE PGSLib_Scan_No_Seg_MODULE,ONLY : PGSLib_SUM_PREFIX


    implicit none
    _DATA_TYPE_,     &
         &   intent(  OUT)   ,     &
         &   dimension(           :)  :: Dest
    _DATA_TYPE_,     &
         &   intent(IN   )   ,     &
         &   dimension(           :)  :: Source

    type (PGSLib_GS_Trace),        &
             OPTIONAL,             &
             POINTER               :: TRACE

    ! Local variables
    type (PGSLib_GS_Trace), pointer :: Local_Trace
    logical (PGSLib_Log_Type)       :: New_Trace

    ! Local_Index is used, but gets saved if we save the trace.
    integer (PGSLib_Int_Type),     &
         &   pointer,              &
         &   dimension(        :) :: Local_Index

    ! If a TRACE was passed in, use it, otherwise use a local trace
    if (PRESENT(Trace)) then
       New_Trace = .NOT. ASSOCIATED(Trace)
    else
       New_Trace = .TRUE.
    end if

    IF (New_Trace) then
       ! The total sizes of Source and Dest must be the same, even
       ! though the local sizes may not be the same.  (That is, they
       ! will have different distributions.)
       
       if (PGSLib_Global_SUM(SIZE(Dest,1)) < PGSLib_Global_SUM(SIZE(Source,1))) then
          call pgslib_error("Total size Dest not large enough in PGSLib_Redistribute")
       end if

       ! We need an index.  That is just the identity, since we are not reordering
       ALLOCATE(Local_Index(SIZE(SOURCE,1)))
       Local_Index = 1
       Local_Index = PGSLib_SUM_PREFIX(Local_Index)

       Local_Trace => PGSLib_Setup_Trace(Local_Index, SIZE(Dest,1))
       if (PRESENT(Trace)) Trace => Local_Trace

       Local_Trace%Global_Index_1 => Local_Index
    ELSE
       Local_Trace => Trace
       Local_Index => Local_Trace%Global_Index_1
    END IF


    Dest = _OP_ID_
    call _SCATTER_OP_      (DEST   = Dest,        &
                            SOURCE = Source,      &
                            INDEX  = Local_Index, &
                            TRACE  = Local_Trace)


    IF (.NOT. PRESENT(TRACE)) THEN
       call PGSLib_Deallocate_Trace(Local_Trace)
       NULLIFY(Local_Index)
    END IF
    return

  END SUBROUTINE _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _SCATTER_OP_
#undef _OP_ID_
