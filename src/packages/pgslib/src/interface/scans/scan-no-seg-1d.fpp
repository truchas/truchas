!!CPP!! $Id: scan-no-seg-1d.fpp,v 1.2 2001/03/22 00:26:13 ferrell Exp $

  function _ROUTINE_NAME_ (INPUT_ARRAY, DIM, MASK, SCOPE)
    implicit none

    ! Subroutine arguments
    _SCAN_DATA_TYPE_ , intent(IN   ),          &
             TARGET,                           &
             dimension(            :)       :: INPUT_ARRAY
    integer (PGSLib_Int_Type), intent(IN   ),  &
             optional                       :: DIM
    logical (PGSLib_Log_Type), intent(IN   ),  &
             optional,                         &
             dimension(SIZE(INPUT_ARRAY,1)) :: MASK
    type (PGSLib_SCOPE),       intent(IN   ),  &
             optional                       :: SCOPE

    ! Function type
    _SCAN_DATA_TYPE_ ,dimension(SIZE(INPUT_ARRAY,1)) :: _ROUTINE_NAME_

    ! Local variables

    integer                                 :: Local_N, Src_Seg, Dest_Seg
    _SCAN_DATA_TYPE_                        :: Src_Data, Dest_Data
    _SCAN_DATA_TYPE_               ,  &
         &   dimension(SIZE(INPUT_ARRAY,1)) :: Dest_ARRAY

    _SCAN_DATA_TYPE_               ,           &
         &   POINTER,                          &
         &   dimension(            :)       :: Src_Array



    Local_N      = SIZE(INPUT_ARRAY, 1)

    ! Check that we got valid argument combinations
    if (PRESENT(DIM)) then
       if (DIM /= 1) then
          call PGSLib_Fatal_ERROR('In SUM_PREFIX if DIM is present it must == 1')
       end if
    end if

! In principle should be able to just point at INPUT_Array if we do not have a MASK.
! Unfortunately, Lahey does not like that when Input_Array is make up on the fly
! with something like (/ (i, i = 1, n) /).

    ALLOCATE(Src_Array(Local_N))
    if (PRESENT(MASK)) then
       Src_Array = MERGE(INPUT_Array, _OP_ID_, MASK=MASK)
    else
       Src_Array = INPUT_ARRAY
    end if

    ! First we need the on-node scan (either PREFIX or SUFFIX, depending on specific routine)
    ! This returns Dest_Array = SCAN(Dest_ARRAY, Src_ARRAY)
    call _ON_NODE_SCAN_(Dest_ARRAY, Src_ARRAY)

    Global: If (PGSLib_Scope_Check(SCOPE) == PGSLib_Global) then
       ! Next we do the accross-node sum_prefix
       ! We do not have to worry about segments, since this is the
       ! no segment case.  We can use FALSE for the segment part
       ! and still use the segmented off-pe sum-prefix routine.
       ! If this PE has 0 sized array, then we use 0 for the
       ! value to be scanned.

       if (Local_N > 0) then
          Src_Data = Dest_Array(_LAST_)
       else
          ! For 0 sized arrays, Upper is identity
          Src_Data = _OP_ID_
       end if
       Src_Seg  = PGSLib_FALSE

       call _OFF_NODE_SCAN_(Dest_Data, Dest_Seg, Src_Data, Src_Seg)

       ! Finally we need to patch Dest_Array with the contribution "donated"
       ! by the lower PEs.
       ! The data to use is from Dest_Data and Dest_Seg

       _ROUTINE_NAME_ = Dest_Array + Dest_Data

    END If Global

    ! We are done
    DEALLOCATE(Src_Array)

    RETURN
    END FUNCTION _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _SCAN_DATA_TYPE_
#undef _OP_ID_
#undef _ON_NODE_SCAN_
#undef _OFF_NODE_SCAN_
#undef _LAST_
