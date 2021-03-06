!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is included for the var_vector FLATTEN SCALAR routines


  function _ROUTINE_NAME_(V_V_SCALAR) RESULT(Values)
    !====================================================================
    ! Purpose(s):
    !   Return a pointer to vector of values in the varying vector
    !====================================================================

    ! Arguments
    type (_VAR_DATA_TYPE_), &
!         TARGET,            &
         INTENT(in   )     :: V_V_SCALAR

    _DATA_TYPE_,           &
         POINTER,          &
         dimension(:)      :: Values

    ! This returs the vector of values in the Varying vector
    Values => V_V_Scalar%V
    RETURN
  END function _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_
