!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE EE_GATHER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Supply the EE gather routines
  !
  !=======================================================================
  use kinds, only: r8
  use gather_module,  only: GATHER
  use gs_info_module, only: EE_All_Ngbr_Trace,      &
                            EE_Mask_Initialized,    &
                            EL_Nbr_Mask
                            
                            
  use gs_util,        only: gs_init_ee_mask,        &
                            ee_gs_init
  use mesh_module,    only: Mesh
  use mesh_parameter_module, only: ncells, nfc
  use var_vector_module
  use truchas_logging_services

  implicit none
  private

  PUBLIC :: EE_Gather

  ! Interface blocks
  INTERFACE EE_GATHER
     !=======================================================================
     ! Purpose - 
     !   Gather integer cell-centered data from cell face neighbors
     !   into a integer cell-centered vector of length nfc
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - Cell-centered quantity 
     !   Output: 
     !          Dest - Cell-centered vector of length nfc containing
     !                 the value of Src for all face neighbor cells
     !          
     !=======================================================================
     MODULE PROCEDURE EE_GATHER_INT
     MODULE PROCEDURE EE_GATHER_SINGLE
     MODULE PROCEDURE EE_GATHER_DOUBLE
     MODULE PROCEDURE EE_GATHER_LOG
     MODULE PROCEDURE EE_GATHER_V_V_INT
     MODULE PROCEDURE EE_GATHER_V_V_SINGLE
     MODULE PROCEDURE EE_GATHER_V_V_DOUBLE
     MODULE PROCEDURE EE_GATHER_V_V_LOG
     MODULE PROCEDURE EE_GATHER_ALL_V_S_INT
     MODULE PROCEDURE EE_GATHER_ALL_V_S_REAL
     MODULE PROCEDURE EE_GATHER_ALL_V_S_LOG
  END INTERFACE

CONTAINS
    !========== EE_GATHER_INT========================================

#define _ROUTINE_NAME_  EE_GATHER_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         0
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc,ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather.fpp"

  !========== EE_GATHER_SINGLE========================================

#define _ROUTINE_NAME_  EE_GATHER_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         0.0
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc,ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather.fpp"

  !========== EE_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_  EE_GATHER_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         0.0_r8
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc,ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather.fpp"

  !========== EE_GATHER_LOG========================================

#define _ROUTINE_NAME_  EE_GATHER_LOG
#define _DATA_TYPE_     logical
#define _OP_ID_         .false.
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc,ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather.fpp"

  !========== EE_GATHER_V_V_INT ========================================

#define _ROUTINE_NAME_  EE_GATHER_V_V_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         0
#define _SRC_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:,:))

#include "ee_gather.fpp"

  !========== EE_GATHER_V_V_SINGLE ========================================

#define _ROUTINE_NAME_  EE_GATHER_V_V_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         0.0
#define _SRC_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:,:))

#include "ee_gather.fpp"

  !========== EE_GATHER_V_V_DOUBLE ========================================

#define _ROUTINE_NAME_  EE_GATHER_V_V_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         0.0_r8
#define _SRC_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:,:))

#include "ee_gather.fpp"

  !========== EE_GATHER_V_V_LOG ========================================

#define _ROUTINE_NAME_  EE_GATHER_V_V_LOG
#define _DATA_TYPE_     logical
#define _OP_ID_         .false.
#define _SRC_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nfc, ncells))
#define _BDY_DIMENSION_ _DIMENSION_((:,:))

#include "ee_gather.fpp"

  !========== EE_GATHER_ALL_V_S_INT========================================

#define _ROUTINE_NAME_   EE_Gather_ALL_V_S_INT
#define _DEST_DATA_TYPE_ type (INT_VAR_VECTOR)
#define _DATA_TYPE_      integer
#define _SRC_DIMENSION_ _DIMENSION_((:))
#define _DST_DIMENSION_ _DIMENSION_((:))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather_all.fpp"    
  
  !========== EE_GATHER_ALL_V_S_REAL========================================

#define _ROUTINE_NAME_   EE_Gather_ALL_V_S_REAL
#define _DEST_DATA_TYPE_ type (REAL_VAR_VECTOR)
#define _DATA_TYPE_      real(r8)
#define _SRC_DIMENSION_ _DIMENSION_((:))
#define _DST_DIMENSION_ _DIMENSION_((:))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather_all.fpp"    
  
  !========== EE_GATHER_ALL_V_S_LOG========================================

#define _ROUTINE_NAME_   EE_Gather_ALL_V_S_LOG
#define _DEST_DATA_TYPE_ type (LOG_VAR_VECTOR)
#define _DATA_TYPE_      logical
#define _SRC_DIMENSION_ _DIMENSION_((:))
#define _DST_DIMENSION_ _DIMENSION_((:))
#define _BDY_DIMENSION_ _DIMENSION_((:))

#include "ee_gather_all.fpp"    
  

END MODULE EE_GATHER_MODULE


