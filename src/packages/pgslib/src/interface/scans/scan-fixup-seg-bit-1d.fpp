!!CPP!! Patch the data by combining the donated item, which came
!!CPP!! from a global scan, individually with each item in the array.
!!CPP!! Once a segment bit is found we can just copy the rest of the
!!CPP!! input array into the destination, since the donation will not
!!CPP!! cross a segment boundary.

!!CPP!!  $Id: scan-fixup-seg-bit-1d.fpp,v 1.1.1.1 2000/10/11 22:44:30 ferrell Exp $

  subroutine _ROUTINE_NAME_(Dest_Array, Src_Array, Src_Bit, Donor_Data)
    ! Combine donor from lower processors into dest array
    implicit none


    ! Subroutine arguments
    _SCAN_DATA_TYPE_, intent(IN   ),         &
         &   dimension(            :)     :: SRC_ARRAY
    logical (PGSLib_Log_Type), intent(IN   ),&
         &   dimension(SIZE(SRC_ARRAY,1)) :: Src_Bit
    _SCAN_DATA_TYPE_, intent(IN   )          &
         &                                :: Donor_Data
    _SCAN_DATA_TYPE_, intent(  OUT),         &
         &   dimension(SIZE(SRC_ARRAY,1)) :: Dest_ARRAY

    ! Local variables

    integer (PGSLib_Int_Type) :: Local_N, i

    Local_N      = SIZE(SRC_ARRAY, 1)

    ! Nothing to do for 0 sized arrays
    if (Local_N < 1) RETURN

    ! Loop trough elements of Src_Array, combine with donor and put into Dest_Array
    do i = _FIRST_ , _LAST_, _INDEX_INCREMENT_
       ! If we find a start bit which is set, then we have a new segment, so skip out of loop.
       if (Src_Bit(i)) exit
       ! If we got past that, this isn''t start of new segment, so add in donation
       ! We don''t need to keep track of the segment bit because we are dropping out
       ! as soon as we find a new segment
       Dest_ARRAY(i) = Src_ARRAY(i) + Donor_Data
    end do

    ! Now copy top part of src_array into the dest
    Dest_Array(i: _LAST_ : _INDEX_INCREMENT_ ) = Src_Array(i: _LAST_ : _INDEX_INCREMENT_ )

    RETURN
  END SUBROUTINE _ROUTINE_NAME_


#undef _ROUTINE_NAME_
#undef _SCAN_DATA_TYPE_
#undef _FIRST_
#undef _LAST_
#undef _INDEX_INCREMENT_
