!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is included for the var_vector_destory routines

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _VAR_DATA_TYPE_
#error "_VAR_DATA_TYPE_ must be defined before including this file"
#endif

  Subroutine _ROUTINE_NAME_(ARRAY)
    !====================================================================
    ! Purpose(s):
    !   Deallocate ARRAY of varying vectors.  
    !   NOTE: This deallocates each of the varying vectors.  It does not 
    !   deallocate the ARRAY.
    !====================================================================
    
    ! Arguments
    type (_VAR_DATA_TYPE_), &
         TARGET,            &
         DIMENSION(:)              :: ARRAY

    ! Local variables
    integer :: i

    ! If ARRAY has 0 size, then there is nothing to do:
    if (SIZE(ARRAY) >= 1) then
       ! Since the ragged array was produced by pointing to sections of
       ! a single big array, we need merely deallocate that big array.
       ! The first Var_Vector contains a pointer to the container
       deallocate(Array(1)%Container)
       
       ! Now we want to NULLIFY the rest of pointers
       NULLIFY(Array(1)%Container)
       do i = 1, SIZE(ARRAY)
          NULLIFY(Array(i)%v)
          Array(I)%L = -1
       end do

    endif

  end Subroutine _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _VAR_DATA_TYPE_

  
