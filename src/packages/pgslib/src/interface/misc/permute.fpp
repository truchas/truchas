!!CPP!! This file is included by pgslib_permute_module.F
!!CPP!! This file provides the permute functionality

! $Id: permute.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file."
#endif

#ifndef _OP_ID_
#error "_OP_ID_ must be defined before including this file"
#endif

#ifndef _SCATTER_OP_
#error "_SCATTER_OP_ must be defined before including this file"
#endif

#ifndef _COMP_OP_
#error "_COMP_OP_ must be defined before including this file"
#endif

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif


  Subroutine _ROUTINE_NAME_(Dest, Source, Index, Mask, TRACE) 
    USE PGSLib_Type_Module,      ONLY : PGSLib_Int_Type,   &
                                        PGSLib_Real_Type,  &
                                        PGSLib_Double_Type,  &
                                        PGSLib_Log_Type,   &
                                        PGSLib_GS_Trace
    
    USE PGSLib_GS_MODULE,        ONLY : PGSLib_Setup_Trace,      &
         &                              PGSLib_Deallocate_Trace
    USE PGSLIB_User_GS_MODULE,   ONLY : PGSLib_Scatter_SUM, PGSLib_Scatter_OR
    USE PGSLib_Reductions_MODULE,ONLY : PGSLib_Global_SUM, PGSLib_Global_COUNT
    USE PGSLib_Utility_MODULE,   ONLY : pgslib_error


    implicit none
    _DATA_TYPE_,     &
         &   intent(  OUT)   ,     &
         &   dimension(           :)  :: Dest
    _DATA_TYPE_,     &
         &   intent(IN   )   ,     &
         &   dimension(           :)  :: Source
    integer (PGSLib_Int_Type),     &
         &   intent(IN   )   ,     &
         &   dimension(SIZE(Source,1)):: Index

    logical (PGSLib_Log_Type),     &
         &   intent(IN   )   ,     &
         &   OPTIONAL        ,     &
         &   dimension(SIZE(Source,1)):: Mask

    type (PGSLib_GS_Trace),        &
             OPTIONAL,             &
             POINTER               :: TRACE
    ! Local variables
    type (PGSLib_GS_Trace), pointer :: Local_Trace
    logical (PGSLib_Log_Type)       :: New_Trace

    ! Local_Index is used so that Index doesn''t change
    integer (PGSLib_Int_Type),     &
         &   pointer,              &
         &   dimension(        :) :: Local_Index

    ! Dest_Temp is used so that Source and Dest may overlap.
    _DATA_TYPE_,   &
         &   dimension(SIZE(Dest,1))  :: Dest_Temp
    
    ! Source_Surity is used to set Dest_Mask which 
    ! is used test whether a dest got modified or not.
    integer (PGSLib_Int_Type),     &
         &   dimension(SIZE(Index,1)) :: Source_Surity

    integer (PGSLib_Int_Type),     &
         &   dimension(SIZE(Dest,1))  :: Dest_Mask
    
    ! If a TRACE was passed in, use it, otherwise use a local trace
    if (PRESENT(Trace)) then
       New_Trace = .NOT. ASSOCIATED(Trace)
    else
       New_Trace = .TRUE.
    end if

    IF (New_Trace) then
       ! The total sizes of Source and Dest must be the same, even
       ! though the local sizes may not be the same.  (That is, they
       ! may have different distributions.)
       
       if (PRESENT(Mask)) then
          if (PGSLib_Global_SUM(SIZE(Dest,1)) < PGSLib_Global_COUNT(Mask)) then
             call pgslib_error("Total size Dest not large enough in PGSLib_Permute")
          endif
       else
          if (PGSLib_Global_SUM(SIZE(Dest,1)) < PGSLib_Global_SUM(SIZE(Source,1))) then
             call pgslib_error("Total size Dest not large enough in PGSLib_Permute")
          end if
       endif

       ALLOCATE(Local_Index(SIZE(Index,1)))
       Local_Index = Index

       Local_Trace => PGSLib_Setup_Trace(Local_Index, SIZE(Dest,1), MASK=Mask)
       if (PRESENT(Trace)) Trace => Local_Trace

       Local_Trace%Global_Index_1 => Local_Index
    ELSE
       Local_Trace => Trace
       Local_Index => Local_Trace%Global_Index_1
    END IF


    Dest_Temp = _OP_ID_
    Source_Surity = -1
    Dest_Mask = 0
    call _SCATTER_OP_      (DEST = Dest_Temp,     &
                            SOURCE = Source,      &
                            INDEX  = Local_Index, &
                            TRACE=Local_Trace,    &
                            MASK=Mask)
    call PGSLib_scatter_sum(DEST   = Dest_Mask,     &
                            SOURCE = Source_Surity, &
                            INDEX  = Local_Index,   &
                            TRACE  = Local_Trace,   &
                            MASK   = Mask)
    ! Pick up original value for dest, unless was actually changed.
    Dest = MERGE(Dest, Dest_Temp, Dest_Mask == 0)

    IF (.NOT. PRESENT(TRACE)) THEN
       call PGSLib_Deallocate_Trace(Local_Trace)
    END IF
    return

  END SUBROUTINE _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _OP_ID_
#undef _SCATTER_OP_
#undef _COMP_OP_
