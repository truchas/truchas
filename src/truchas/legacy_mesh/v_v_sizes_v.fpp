!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is included for the var_vector SIZES routines

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
