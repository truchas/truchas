!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is included for the var_vector FLATTEN VECTOR routines

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _VAR_DATA_TYPE_
#error "_VAR_DATA_TYPE_ must be defined before including this file"
#endif

  function _ROUTINE_NAME_(ARRAY) RESULT(BigArray)
    !====================================================================
    ! Purpose(s):
    !   Return a pointer to a single large array which contains
    !   all the varying vectors in the ragged array.
    !====================================================================

    ! Arguments
    type (_VAR_DATA_TYPE_), &
         INTENT(in   ),     &
!         TARGET,            &
         DIMENSION(:)             :: ARRAY

    _DATA_TYPE_, POINTER, dimension(:) :: BigArray

    ! Local variable, for returning empty result if necessary
    _DATA_TYPE_, TARGET, SAVE, dimension(0) :: ZeroArray

    ! If ARRAY has 0 size, then return 0 array
    ! You might think "do nothing" would be appropriate, but returning 0 sized
    ! array seems to be route of least surprise.
    if (SIZE(ARRAY) >= 1) then
       ! The first varying vector has a pointer to the head of the BigArray.
       BigArray => Array(1)%CONTAINER
    else
       BigArray => ZeroArray
    end if
    RETURN
  END function _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_

  
