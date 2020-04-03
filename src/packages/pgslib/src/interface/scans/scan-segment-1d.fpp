!!CPP!! This file provides the routines for the segmented
!!CPP!! scans.  The versions here use a parity segment.
!!CPP!! These routines convert the parity segment into
!!CPP!! a start-bit segment, and then call the start-bit
!!CPP!! versions of the segmented scans.

!!CPP!! Because these versions use SHIFTS for the conversions,
!!CPP!! they depend on PGSLib_Shift_Module.

!!CPP!! $Id: scan-segment-1d.fpp,v 1.1.1.1 2000/10/11 22:44:30 ferrell Exp $

  function _ROUTINE_NAME_ (INPUT_ARRAY, SEGMENT, DIM, MASK, SCOPE)
    implicit none

    ! Subroutine arguments
    _SCAN_DATA_TYPE_ , intent(IN   ),          &
             TARGET,                           &
             dimension(            :)       :: INPUT_ARRAY
    logical (PGSLib_Log_Type), intent(IN   ),  &
             dimension(SIZE(INPUT_ARRAY,1)) :: SEGMENT
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

    integer                                 :: Local_N, n, Wrapped_Index
    logical (PGSLib_Log_Type),                 &
             dimension(SIZE(INPUT_ARRAY,1)) :: Segment_BIT
    integer (PGSLib_Int_Type),                 &
             dimension(SIZE(INPUT_ARRAY,1)) :: Array_Index
    logical (PGSLib_Log_Type),                 &
             dimension(SIZE(INPUT_ARRAY,1)) :: Shifted_Segment


    Local_N      = SIZE(INPUT_ARRAY, 1)

    ! Check that we got valid argument combinations
    if (PRESENT(DIM)) then
       if (DIM /= 1) then
          call PGSLib_Fatal_ERROR('In SUM_???FIX if DIM is present it must == 1')
       end if
    end if

    ! Construct a segment bit from a segment.
    ! This varies slightly depending on whether the scope is global or
    ! local.
    if (PGSLib_Scope_Check(SCOPE) == PGSLib_Global) then
       Shifted_Segment = PGSLib_Global_CSHIFT(Segment, SHIFT = _SHIFT_DIRECTION_)
    else
       Shifted_Segment = CSHIFT(Segment, SHIFT = _SHIFT_DIRECTION_)
    end if

    ! Need to make sure that first or last bit gets set correctly.
    ! This is a lot of work just to patch up the first or last element,
    ! but cannot think of a robust way of making this easier.
    Array_Index   = PGSLib_SUM_PREFIX( (/ (1, n=1, Local_N) /), SCOPE=SCOPE)
    Wrapped_Index = _GLOBAL_MAX_OR_MIN_VAL_ (Array_Index, SCOPE=SCOPE)

    WHERE(Array_Index == Wrapped_Index)
       Shifted_Segment = .NOT. Segment
    END WHERE

    SEGMENT_BIT = .NOT. (Shifted_Segment .EQV. Segment)

    ! Now that we have a segment bit, we can call the segment-bit scans
    _ROUTINE_NAME_ = _SEG_BIT_SCAN_ (INPUT_ARRAY, SEGMENT_BIT = Segment_Bit, MASK=MASK, SCOPE=SCOPE)


    ! We are done

    RETURN
    END FUNCTION _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _SCAN_DATA_TYPE_
#undef _SHIFT_DIRECTION_
#undef _GLOBAL_MAX_OR_MIN_VAL_
#undef _SEG_BIT_SCAN_
