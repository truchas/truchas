! This is included for the var_vector_create routines

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _VAR_DATA_TYPE_
#error "_VAR_DATA_TYPE_ must be defined before including this file"
#endif

  Subroutine _ROUTINE_NAME_ (ARRAY, SIZES)
    !====================================================================
    ! Purpose(s):
    !   Allocate ARRAY of varying vectors according to SIZES
    !   NOTE: This allocates each of the vectors.  The array 
    !   is already allocated.
    !====================================================================
#ifdef DEBUG_V_V
     use truchas_logging_services
#endif

    ! Arguments
    type (_VAR_DATA_TYPE_), &
          DIMENSION(:)             :: ARRAY
    integer (int_kind),    &
          intent(IN),      &
          DIMENSION(SIZE(ARRAY,1)) :: SIZES

    ! Local variables
    integer (int_kind) :: i, offset, TotalSize, l_size, upper
    _DATA_TYPE_,           &
             pointer,      &
             dimension(:) :: BigArray

    ! If ARRAY has 0 size, then there is nothing to do:
    if (SIZE(ARRAY) >= 1) then

       ! A ragged array is an array of varying vectors.
       ! For efficiency, we allocate all the varying vectors at once, in
       ! a single large array.  Then we point at the section we want.
       TotalSize = SUM(SIZES)
       Call ArrayCreate(BigArray, 1, TotalSize, 'Create Var_Vector')

       ! We keep a pointer to this BigArray in the first Var_Vector
       ARRAY(1)%Container => BigArray
       ! Now we point the varying vectors at sections of BigArray
       OffSet = 1
       do i = 1, SIZE(ARRAY)
          l_size = SIZES(i)
          upper  = Offset + l_size - 1
          Array(i)%v => BigArray(Offset: upper)
          Offset = Offset + l_size
          Array(I)%L = l_size
#ifdef DEBUG_V_V
          if (SIZE(Array(i)%v) /= l_size) then
             call TLS_panic ('CREATE: sizes incorrect in CREATE var_vector')
          endif
#endif
       end do

    endif

    return
  end Subroutine _ROUTINE_NAME_


#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_
  
