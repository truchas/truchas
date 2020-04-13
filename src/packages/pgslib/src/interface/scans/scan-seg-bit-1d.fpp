!!CPP!! This file provides the core routines for the segmented
!!CPP!! scans.  The versions here use a segment-start-bit.  The version
!!CPP!! which uses parity segments calls this after converting
!!CPP!! the parity-segment to a segment-start-bit.

!!CPP!! $Id: scan-seg-bit-1d.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $
  function _ROUTINE_NAME_ (INPUT_ARRAY, SEGMENT_BIT, MASK, SCOPE)
    implicit none

    ! Subroutine arguments
    _SCAN_DATA_TYPE_ , intent(IN   ),          &
             TARGET,                           &
             dimension(            :)       :: INPUT_ARRAY
    logical (PGSLib_Log_Type), intent(IN   ),  &
             TARGET,                           &
             dimension(SIZE(INPUT_ARRAY,1)) :: SEGMENT_BIT
    logical (PGSLib_Log_Type), intent(IN   ),  &
             optional,                         &
             dimension(SIZE(INPUT_ARRAY,1)) :: MASK
    type (PGSLib_SCOPE),       intent(IN   ),  &
             optional                       :: SCOPE

    ! Function type
    _SCAN_DATA_TYPE_ ,dimension(SIZE(INPUT_ARRAY,1)) :: _ROUTINE_NAME_

    ! Local variables

    integer                                 :: Local_N
    integer                                 :: Src_Seg, Dest_Seg
    _SCAN_DATA_TYPE_                        :: Src_Data, Dest_Data
    _SCAN_DATA_TYPE_               ,  &
         &   dimension(SIZE(INPUT_ARRAY,1)) :: Dest_ARRAY
    logical (PGSLib_Log_Type),                 &
             dimension(SIZE(INPUT_ARRAY,1)) :: Dest_BIT

    _SCAN_DATA_TYPE_               ,           &
             POINTER,                          &
             dimension(            :)       :: Src_Array




    Local_N      = SIZE(INPUT_ARRAY, 1)


    if (PRESENT(MASK)) then
       ALLOCATE(Src_Array(Local_N))
       Src_Array = MERGE(INPUT_Array, _OP_ID_, MASK=MASK)
    else
       Src_Array => INPUT_ARRAY
    end if

    ! First we need the on-node scan (either PREFIX or SUFFIX, depending on specific routine)
    ! This returns Dest_Array = SCAN(Dest_ARRAY, Src_ARRAY)
    call _ON_NODE_SCAN_(Dest_ARRAY, Dest_Bit, Src_ARRAY, SEGMENT_Bit)

    if (PRESENT(MASK)) DEALLOCATE(Src_Array)

    Global: If (PGSLib_Scope_Check(SCOPE) == PGSLib_Global) then
       ! Next we do the accross-node sum_prefix
       ! The global scan takes the high or low item, depending on
       ! whether this is up (PREFIX) or down (SUFFIX).
       ! The scan is done on that item.  The return data
       ! is then scanned into the rest of the data.

       ! If this PE has 0 sized array, then we use 0 for the
       ! value to be scanned and FALSE for the segment bit.

       if (Local_N > 0) then
          Src_Data = Dest_Array(_LAST_)
          Src_Seg  = MERGE(PGSLib_TRUE, PGSLib_FALSE, Dest_Bit(_LAST_) )
       else
          ! For 0 sized arrays, Upper is identity
          Src_Data = _OP_ID_
          Src_Seg  = PGSLib_FALSE
       end if


       call _OFF_NODE_SCAN_(Dest_Data, Dest_Seg, Src_Data, Src_Seg)

       ! Finally we need to patch Dest_Array with the contribution "donated"
       ! by the lower PEs.
       ! The data to use is from Dest_Data and Dest_Seg.

       call _SCAN_FIXUP_ (_ROUTINE_NAME_, Dest_Array, Dest_Bit, Dest_Data)

    END If Global

    ! We are done

    RETURN
    END FUNCTION _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _SCAN_DATA_TYPE_
#undef _OP_ID_
#undef _ON_NODE_SCAN_
#undef _OFF_NODE_SCAN_
#undef _SCAN_FIXUP_
#undef _LAST_
