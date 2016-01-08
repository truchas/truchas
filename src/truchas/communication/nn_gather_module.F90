!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NN_GATHER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Supply the node<->node (NN) gather routines
  !
  !=======================================================================
  use kinds, only: r8
  use gs_info_module, only: NN_All_Ngbr_Trace
  use mesh_module,    only: Vertex_Ngbr_All
  use parameter_module
  use var_vector_module
  use truchas_logging_services
  implicit none
  private

  PUBLIC :: NN_Gather, NN_Gather_BoundaryData

  ! Interface blocks
  INTERFACE NN_GATHER
     !=======================================================================
     ! Purpose - 
     !   Gather node data from all neighboring nodes
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - Node-based quantity
     !   Output: 
     !          Dest - Node based var-vector
     !          
     !=======================================================================
     MODULE PROCEDURE NN_GATHER_ALL_V_S_INT
     MODULE PROCEDURE NN_GATHER_ALL_V_S_REAL
     MODULE PROCEDURE NN_GATHER_ALL_V_S_LOG
  END INTERFACE

CONTAINS
  !========== NN_GATHER_ALL_V_S_INT========================================

#define _ROUTINE_NAME_   NN_Gather_ALL_V_S_INT
#define _DEST_DATA_TYPE_ type (INT_VAR_VECTOR)
#define _DATA_TYPE_      integer
#define _SRC_DIMENSION_ _DIMENSION_((:))
#define _DST_DIMENSION_ _DIMENSION_((:))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "nn_gather_all.fpp"    
  
  !========== NN_GATHER_ALL_V_S_REAL========================================

#define _ROUTINE_NAME_   NN_Gather_ALL_V_S_REAL
#define _DEST_DATA_TYPE_ type (REAL_VAR_VECTOR)
#define _DATA_TYPE_      real(r8)
#define _SRC_DIMENSION_ _DIMENSION_((:))
#define _DST_DIMENSION_ _DIMENSION_((:))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "nn_gather_all.fpp"    
  
  !========== NN_GATHER_ALL_V_S_LOG========================================

#define _ROUTINE_NAME_   NN_Gather_ALL_V_S_LOG
#define _DEST_DATA_TYPE_ type (LOG_VAR_VECTOR)
#define _DATA_TYPE_      logical
#define _SRC_DIMENSION_ _DIMENSION_((:))
#define _DST_DIMENSION_ _DIMENSION_((:))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "nn_gather_all.fpp"    
  

  subroutine NN_Gather_BoundaryData (BOUNDARY, SOURCE)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Gather only the off-processor data which is needed for 
    !    a node-node gather.
    !
    !   The functionality of this bit of code is included in
    !   NN_GATHER, but NN_GATHER does more.  This bit
    !   just gets the neighboring processor data, but returns them
    !   in the BOUNDARY buffer.  
    !   NN_GATHER continues by filling them into an ncells  var-vector
    !   temporary.  Use NN_GATHER if you want to do vector operations
    !   with the temporary.  Use this if you have to do indirect
    !   addressing yourself.
    !
    !---------------------------------------------------------------------------
    use pgslib_module,  Only: PGSLib_Size_Of_Dup, PGSLib_Dup_Index, &
                              PGSLib_Size_Of_Sup, PGSLib_gather_Buffer
    ! arguments
    real(r8), Dimension(:), INTENT(IN   ) :: SOURCE
    real(r8), Dimension(:), POINTER       :: BOUNDARY
    
    ! local variables
    Integer :: status, i
    real(r8), POINTER, Dimension(:) :: Duplicate_Data

    !---------------------------------------------------------------------------
    ! if BOUNDARY doesn't point at anything, get some space
    ! if BOUNDARY is already allocated, don't do anything.

    If (.Not. Associated(BOUNDARY)) Then
       ALLOCATE(Duplicate_Data(PGSLib_Size_Of_Dup(NN_All_Ngbr_Trace)))
       do i = 1, PGSLib_Size_Of_Dup(NN_All_Ngbr_Trace)
          Duplicate_Data(i) = SOURCE(PGSLib_Dup_Index(NN_All_Ngbr_Trace, i))
       end do

       Allocate(BOUNDARY(PGSLib_Size_OF_Sup(NN_All_Ngbr_Trace)), STAT=status)
       call TLS_fatal_if_any (status /= 0, 'NN_Gather_BoundaryData: error allocating BOUNDARY')

       ! the communication, takes data from Duplicate buffer, puts data into BOUNDARY
       BOUNDARY = PGSLib_gather_Buffer(Duplicate_Data, NN_All_Ngbr_Trace)
       ! Clean up
       DEALLOCATE(Duplicate_Data)
    End If

  End subroutine NN_Gather_BoundaryData

END MODULE NN_GATHER_MODULE
