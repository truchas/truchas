!!CPP!! This file is included for the gather buffer routines in pgslib_gs_comm_module.F

!!CPP!! $Id: buff_gather.fpp,v 1.2 2001/03/22 00:26:13 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _GATHER_BUF_C_
#error "_GATHER_BUF_C_ must be defined before including this file"
#endif

#ifndef _DUP_DIMENSION_
#error "_DUP_DIMENSION_ must be defined before including this file"
#endif

#ifndef _SUP_DIMENSION_
#error "_SUP_DIMENSION_ must be defined before including this file"
#endif

#ifndef _BLOCKSIZE_
#error "_BLOCKSIZE_ must be defined before including this file"
#endif

    subroutine _ROUTINE_NAME_ (Supplement, Duplicate, Trace) 
    Use pgslib_timing_module
    USE PGSLib_Type_MODULE
    USE PGSLib_Utility_MODULE, ONLY : PGSLib_Error
    IMPLICIT NONE

    _DATA_TYPE_,         INTENT(IN   ), &
                        _DUP_DIMENSION_ :: Duplicate
    type (PGSLib_GS_Trace), INTENT(IN)  :: Trace
    _DATA_TYPE_,         INTENT(  OUT), &
                        _SUP_DIMENSION_ :: Supplement


    ! Local variables
    integer (PGSLib_INT_TYPE) :: NSupplement, NDuplicate, BlockSize

#ifdef USE_TIMERS_1
    call Enter_Routine(GATHER_BUFFER_STATISTICS())
    call pgslib_start_timer(Trace%GatherBuffer_Timer)
#endif

#ifdef SAFETY_1

    ! First check that the arrays are the right size
    ! This funky piece of code gets the length of the right most dimension
    NSupplement = SIZE(Supplement, SIZE(SHAPE(Supplement)))
    IF (NSupplement .NE. PGSLIB_Size_Of_Sup(Trace)) THEN
       Call PGSLib_Error("Wrong size for Supplement in PGSLib_Gather_Buffer")
    ENDIF
    
    ! This funky piece of code gets the length of the right most dimension
    NDuplicate = SIZE(Duplicate, SIZE(SHAPE(Duplicate)))
    IF (NDuplicate .NE. PGSLib_Size_Of_Dup(Trace)) THEN
       Call PGSLib_Error("Wrong size for Duplicate in PGSLib_Gather_Buffer")
    ENDIF

#endif

    ! Blocks are always on the left most axes.
    BlockSize = _BLOCKSIZE_

    ! Now do the communication
    Call _GATHER_BUF_C_ (Supplement, Duplicate, BlockSize, Trace%GS_Trace)

#ifdef USE_TIMERS_1
    call Exit_Routine(GATHER_BUFFER_STATISTICS())
    call pgslib_stop_timer(Trace%GatherBuffer_Timer)
#endif
    RETURN

    end subroutine _ROUTINE_NAME_

#undef _BLOCKSIZE_
#undef _DATA_TYPE_
#undef _GATHER_BUF_C_
#undef _DUP_DIMENSION_
#undef _SUP_DIMENSION_
#undef _ROUTINE_NAME_
