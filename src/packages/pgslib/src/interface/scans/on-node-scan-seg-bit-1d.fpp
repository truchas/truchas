!!CPP!! Provides the on-node, segmented (w/ start-bit) scans

!!CPP!! $Id: on-node-scan-seg-bit-1d.fpp,v 1.1.1.1 2000/10/11 22:44:30 ferrell Exp $

#ifndef _SCAN_DATA_TYPE_
#error "_SCAN_DATA_TYPE_ must be defined before including this file."
#endif

#ifndef _OP_ID_
#error "_OP_ID_ must be defined before including this file"
#endif

#ifndef _FIRST_
#error "_FIRST_ must be defined before including this file"
#endif

#ifndef _LAST_
#error "_LAST_ must be defined before including this file"
#endif

#ifndef _INDEX_INCREMENT_
#error "_INDEX_INCREMENT_ must be defined before including this file"
#endif

  subroutine _ROUTINE_NAME_(Dest_ARRAY, Dest_Bit, Src_ARRAY, Src_Bit)
    use PGSLib_Type_Module
    ! On Node segmented scan operation
    implicit none

    ! Subroutine arguments
    _SCAN_DATA_TYPE_, intent(IN   ),         &
         &   dimension(            :)     :: Src_ARRAY
    LOGICAL,          intent(IN   ),         &
         &   dimension(SIZE(Src_ARRAY,1)) :: Src_Bit
    _SCAN_DATA_TYPE_, intent(  OUT),         &
         &   dimension(SIZE(Src_ARRAY,1)) :: Dest_ARRAY
    LOGICAL,          intent(  OUT),         &
         &   dimension(SIZE(Src_ARRAY,1)) :: Dest_Bit
    ! Local variables
    integer (PGSLib_Int_TYPE) :: Local_N, i
    _SCAN_DATA_TYPE_          :: Src_Data, Dest_Data
    logical                   :: Src_Seg, Dest_Seg

    Local_N      = SIZE(Src_ARRAY, 1)

    ! Nothing to do for 0 sized arrays
    if (Local_N < 1) RETURN

    Dest_Array( _FIRST_ ) = Src_Array( _FIRST_ )
    Dest_Bit  ( _FIRST_ ) = Src_Bit  ( _FIRST_ )
    DO i = _FIRST_ + _INDEX_INCREMENT_, _LAST_, _INDEX_INCREMENT_

       Src_Data  = Src_Array(i)
       Src_Seg   = Src_Bit  (i)
       Dest_Data = Dest_Array( i - _INDEX_INCREMENT_ )
       Dest_Seg  = Dest_Bit  ( i - _INDEX_INCREMENT_ )
       
       Dest_Array(i) = Src_Data + MERGE(Dest_Data, _OP_ID_ , .NOT. Src_Seg)
       Dest_Bit  (i) = Src_Seg .OR. Dest_Seg

    END DO

    RETURN
  END subroutine _ROUTINE_NAME_
    
#undef _ROUTINE_NAME_
#undef _SCAN_DATA_TYPE_
#undef _OP_ID_
#undef _FIRST_
#undef _LAST_
#undef _INDEX_INCREMENT_
