!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Scan Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! This file contains the DUMMY, serial versions.  This file
!     is replaced by a .c file containing MPI calls in the
!     actual library.
!
!


! $Id: !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file."
#endif

#ifndef PGSLIB_DATA_TYPE
#error "_PGSLIB_DATA_TYPE_ must be defined before including this file."
#endif

     subroutine _ROUTINE_NAME_ (Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       ! Since this is serial emulation, there is no contribution from other PEs
       implicit none
       PGSLIB_DATA_TYPE, intent(INOUT) :: Dest_Data, Src_Data
       integer,          intent(INOUT) :: Dest_Seg, Src_Seg

       Dest_Data = 0
       Dest_Seg  = 0
       
       return
     end subroutine _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef PGSLIB_DATA_TYPE
