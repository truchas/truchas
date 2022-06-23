!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is included for the var_vector SIZES routines

  function _ROUTINE_NAME_(Scalar) RESULT(SIZES)
    !====================================================================
    ! Purpose(s):
    !   Return the size of the var_vector
    !====================================================================

    ! Arguments
    type (_VAR_DATA_TYPE_),  &
          INTENT(in    )          :: Scalar

    ! Result
    integer :: SIZES

    SIZES = Scalar%L

  end function _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_
