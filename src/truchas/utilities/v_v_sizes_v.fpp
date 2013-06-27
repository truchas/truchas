! This is included for the var_vector SIZES routines

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _VAR_DATA_TYPE_
#error "_VAR_DATA_TYPE_ must be defined before including this file"
#endif

  function _ROUTINE_NAME_(ARRAY) RESULT(SIZES)
    !====================================================================
    ! Purpose(s):
    !   Return the any array of the sizes of var_vectors
    !====================================================================
    
    ! Arguments
    type (_VAR_DATA_TYPE_), &
         INTENT(in   ),    &
         DIMENSION(:)             :: ARRAY

    ! Result
    integer, DIMENSION(SIZE(ARRAY,1)) :: SIZES

    ! Local variables
    integer :: i

    do i = 1, SIZE(ARRAY)
       SIZES(i) = SIZE(ARRAY(i)%V)
    end do

  end function _ROUTINE_NAME_
    
#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_

  
