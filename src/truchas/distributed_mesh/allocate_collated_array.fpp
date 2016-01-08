!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# /* allocate_collated_array -- include file for parallel_communication.F90 */

#ifdef _TYPE_
#undef _TYPE_
#endif

#ifdef _INTEGER_DATA_
#undef _INTEGER_DATA_
#define _TYPE_ integer
#define _PROC1_ allocate_CA_I1
#define _PROC2_ allocate_CA_I2
#endif

#ifdef _SINGLE_DATA_
#undef _SINGLE_DATA_
#define _TYPE_ real
#define _PROC1_ allocate_CA_S1
#define _PROC2_ allocate_CA_S2
#endif

#ifdef _DOUBLE_DATA_
#undef _DOUBLE_DATA_
#define _TYPE_ real(kind(1.0d0))
#define _PROC1_ allocate_CA_D1
#define _PROC2_ allocate_CA_D2
#endif

#ifdef _LOGICAL_DATA_
#undef _LOGICAL_DATA_
#define _TYPE_ logical
#define _PROC1_ allocate_CA_L1
#define _PROC2_ allocate_CA_L2
#endif

#ifndef _TYPE_
#error "one of LOGICAL_DATA, INTEGER_DATA, SINGLE_DATA, DOUBLE_DATA must be defined"
#endif

  subroutine _PROC1_ (array, size1, stat)
  
    _TYPE_, pointer    :: array(:)
    integer, intent(in) :: size1
    integer, intent(out), optional :: stat
    
    if (is_IOP) then
      if (present(stat)) then
        allocate(array(size1), stat=stat)
      else
        allocate(array(size1))
      end if
    else
      allocate(array(0))
    end if
    
    if (present(stat)) call broadcast (stat)
    
  end subroutine _PROC1_
  
  subroutine _PROC2_ (array, size1, size2, stat)
  
    _TYPE_, pointer    :: array(:,:)
    integer, intent(in) :: size1, size2
    integer, intent(out), optional :: stat
    
    if (is_IOP) then
      if (present(stat)) then
        allocate(array(size1,size2), stat=stat)
      else
        allocate(array(size1,size2))
      end if
    else
      allocate(array(size1,0))
    end if
    
    if (present(stat)) call broadcast (stat)
    
  end subroutine _PROC2_
  
#undef _TYPE_
#undef _PROC1_
#undef _PROC2_
