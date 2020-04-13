!!CPP!! This file is included for the scatter buffer routines in pgslib_gs_comm_module.F

!!CPP!! $Id: scatter_buff.fpp,v 1.2 2001/03/22 00:26:13 ferrell Exp $

  function _ROUTINE_NAME_ (Supplement, Trace) RESULT(Duplicate)
    IMPLICIT NONE
    _DATA_TYPE_,            INTENT(IN), &
                       _SUP_DIMENSION_  :: Supplement
    type (PGSLib_GS_Trace), INTENT(IN)  :: Trace

    _DATA_TYPE_,                        &
                       _DUP_DIMENSION_ :: Duplicate

    ! Local variables
    integer (PGSLib_Int_TYPE) :: NSupplement, NDuplicate, BlockSize

#ifdef USE_TIMERS_1
    call Enter_Routine(SCATTER_BUFFER_STATISTICS())
    call pgslib_start_timer(Trace%ScatterBuffer_Timer)
#endif
    ! First check that the arrays are the right size
    ! This funky piece of code gets the length of the right most dimension
    NDuplicate = SIZE(Duplicate, SIZE(SHAPE(Duplicate)))
    IF (NDuplicate .NE. PGSLib_Size_Of_Dup(Trace)) THEN
       Call PGSLib_Error("Wrong size for Duplicate in PGSLib_Scatter_Buffer")
    ENDIF

    ! This funky piece of code gets the length of the right most dimension
    NSupplement = SIZE(Supplement, SIZE(SHAPE(Supplement)))
    IF (NSupplement .NE. PGSLib_Size_Of_Sup(Trace)) THEN
       Call PGSLib_Error("Wrong size for Supplement in PGSLib_Scatter_Buffer")
    ENDIF

    ! Now do the communication
    BlockSize = _BLOCKSIZE_
    Call _SCATTER_BUF_C_ (Duplicate, Supplement, BlockSize, Trace%GS_Trace)

#ifdef USE_TIMERS_1
    call Exit_Routine(SCATTER_BUFFER_STATISTICS())
    call pgslib_stop_timer(Trace%ScatterBuffer_Timer)
#endif
    RETURN

  end function _ROUTINE_NAME_
#undef _BLOCKSIZE_
#undef _DATA_TYPE_
#undef _SCATTER_BUF_C_
#undef _DUP_DIMENSION_
#undef _SUP_DIMENSION_
#undef _ROUTINE_NAME_
