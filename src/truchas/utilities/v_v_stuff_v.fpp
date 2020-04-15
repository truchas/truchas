!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is included for the var_vector STUFF VECTOR routines

  subroutine _ROUTINE_NAME_(RAGGED_ARRAY, SOURCE)
    !====================================================================
    ! Purpose(s):
    !   STUFF the values of SOURCE into RAGGED_ARRAY.  In this imlementation
    !   this is efficient because RAGGED_ARRAY has a large container.
    !   If RAGGED_ARRAY is not allocated, then it is an error to call this routine.
    !   It is required that SIZE(SOURCE) == SIZES(RAGGED_ARRAY)
    !====================================================================

    ! Arguments
    type (_VAR_DATA_TYPE_), &
         INTENT(inout),    &
         DIMENSION(:)             :: RAGGED_ARRAY

    _DATA_TYPE_,    &
         INTENT(IN   ),    &
         DIMENSION(:)             :: SOURCE

    ! If ARRAY has 0 size, then there is nothing to do:
    if (SIZE(RAGGED_ARRAY) >= 1) then
       ! Check that sizes are commensurate
       if (SUM(SIZES(RAGGED_ARRAY)) /= SIZE(SOURCE)) then
          call TLS_panic ('varying vector STUFF: sizes not commensurate')
       end if

       ! Since in this module we are controlling the implementation, the following
       ! is okay.  If the implementation changes, this will have to change.
       RAGGED_ARRAY(1)%CONTAINER = SOURCE
    endif
    RETURN
  END subroutine _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_
