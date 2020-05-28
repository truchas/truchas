!!CPP!! $Id: on-node-scan-no-seg-1d.fpp,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

  subroutine _ROUTINE_NAME_(Dest_ARRAY, Src_ARRAY)
    ! On Node scan operation
    implicit none

    ! Subroutine arguments
    _SCAN_DATA_TYPE_, intent(IN   ),         &
         &   dimension(            :)     :: Src_ARRAY
    _SCAN_DATA_TYPE_, intent(  OUT),         &
         &   dimension(SIZE(Src_ARRAY,1)) :: Dest_ARRAY
    ! Local variables

    integer (PGSLib_Int_TYPE) :: Local_N, i

    Local_N      = SIZE(Src_ARRAY, 1)

    ! Nothing to do for 0 sized arrays
    if (Local_N < 1) RETURN

    Dest_Array( _FIRST_ ) =     Src_Array( _FIRST_ )
    DO i = _FIRST_ + _INDEX_INCREMENT_, _LAST_, _INDEX_INCREMENT_
       Dest_ARRAY(i) = Src_ARRAY(i) + Dest_ARRAY(i - _INDEX_INCREMENT_ )
    END DO

    RETURN
  END subroutine _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _SCAN_DATA_TYPE_
#undef _OP_ID_
#undef _FIRST_
#undef _LAST_
#undef _INDEX_INCREMENT_
