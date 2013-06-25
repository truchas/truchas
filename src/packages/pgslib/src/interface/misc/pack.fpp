!!CPP!! This file is included by pgslib_permute_module.F
!!CPP!! This file provides the PACK functionality

! $Id: pack.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file."
#endif

#ifndef _OP_ID_
#error "_OP_ID_ must be defined before including this file"
#endif

#ifndef _SCATTER_OP_
#error "_SCATTER_OP_ must be defined before including this file"
#endif

#ifndef _RESULT_
#error "_RESULT_ must be defined before including this file"
#endif

    USE PGSLib_Type_Module,      ONLY : PGSLib_Int_Type,   &
                                        PGSLib_Real_Type,  &
                                        PGSLib_Double_Type,&
                                        PGSLib_Log_Type,   &
                                        PGSLib_GS_Trace
    USE PGSLib_GS_MODULE,        ONLY : PGSLib_Setup_Trace, PGSLib_Deallocate_Trace
    USE PGSLIB_User_GS_MODULE,   ONLY : _SCATTER_OP_
    USE PGSLib_Reductions_MODULE,ONLY : PGSLib_Global_SUM, PGSLib_Global_COUNT
    USE PGSLib_Scan_No_Seg_MODULE,ONLY : PGSLib_Sum_Prefix
    USE PGSLib_Utility_MODULE,   ONLY : pgslib_error

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
