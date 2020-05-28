!!CPP!! This file is included by pgslib_permute_module.F
!!CPP!! This file provides the PACK functionality

! $Id: pack.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

    implicit none

    _DATA_TYPE_,  &
         &   intent(INOUT),     &
         &   dimension(:)       :: _RESULT_
    _DATA_TYPE_,  &
         &   intent(IN   ),     &
         &   dimension(:)       :: Source
    logical (PGSLib_Log_Type),  &
         &   intent(IN   ),     &
         &   dimension(SIZE(Source,1)) :: MASK

    ! local variables
    integer (PGSLib_Int_Type) :: N_PAcked, N_Packed_Tot
    integer (PGSLib_Int_Type),  &
         &   dimension(SIZE(Source,1)) :: Index, Enum
    type (PGSLib_GS_Trace), POINTER :: Trace


    ! Check that the size is large enough
    N_Packed = SIZE(_RESULT_,1)
    if (PGSLib_Global_SUM(N_Packed) < PGSLib_Global_COUNT(MASK)) then
       call pgslib_error("Size of Result is not large enough for PGSLib_PACK")
    end if

    ! Enumerate the values to packe
    Enum = MERGE(1,0,MASK)
    Index = PGSLib_SUM_PREFIX(Enum)

    ! Setup the trace
    Trace => PGSLib_Setup_Trace(Index, N_Packed, MASK=MASK)

    ! Do the packing
    _RESULT_ = _OP_ID_
    call _SCATTER_OP_(DEST   = _RESULT_, &
                      SOURCE = Source, &
                      INDEX  = Index,   &
                      TRACE  = Trace,   &
                      MASK   = MASK)

    ! Release the trace and we are done
    call PGSLib_Deallocate_Trace(Trace)

    return

#undef _DATA_TYPE_
#undef _SCATTER_OP_
#undef _OP_ID_
#undef _RESULT_
