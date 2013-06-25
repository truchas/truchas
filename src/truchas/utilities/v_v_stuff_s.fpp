! This is included for the var_vector STUFF VECTOR routines

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _VAR_DATA_TYPE_
#error "_VAR_DATA_TYPE_ must be defined before including this file"
#endif

  subroutine _ROUTINE_NAME_(V_V_SCALAR, SOURCE)
    !====================================================================
    ! Purpose(s):
    !   STUFF the values of SOURCE into V_V_SCALAR.  This hides the
    !   implementatioh of varying vectors.
    !   If V_V_SCALAR is not allocated, then it is an error to call this routine.
    !   It is required that SIZE(SOURCE) == SIZES(V_V_SCALAR)
    !====================================================================
    
    ! Arguments
    type (_VAR_DATA_TYPE_), &
         INTENT(INOUT)     :: V_V_SCALAR

    _DATA_TYPE_,           &
         INTENT(IN   ),    &
         DIMENSION(:)      :: SOURCE

    ! Check that sizes are commensurate
    if (SIZES(V_V_SCALAR) /= SIZE(SOURCE))then
       call TLS_panic ('varying vector STUFF: sizes not commensurate')
    end if

    ! Since in this module we are controlling the implementation, the following
    ! is okay.  If the implementation changes, this will hve to change.
    V_V_SCALAR%V = SOURCE
    RETURN
  END subroutine _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_

  
