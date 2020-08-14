MODULE PGSLib_GS_Comm_MODULE
  ! Provide the communication support for the gather and scatter operations
  ! $Id: pgslib_gs_comm_module.F,v 1.2 2001/03/22 00:26:13 ferrell Exp $
  use PGSLib_Stats
  Use pgslib_timing_module
  USE PGSLib_Type_MODULE
  USE PGSLib_Utility_MODULE, ONLY : PGSLib_Error
  use pgslib_c_binding
  IMPLICIT NONE 
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_Gather_Buffer, PGSLib_Buffer_Gather
  PUBLIC :: PGSLib_Scatter_Buffer

  INTERFACE PGSLib_Gather_Buffer
     MODULE PROCEDURE PGSLib_Gather_Buf_1_1_INT
     MODULE PROCEDURE PGSLib_Gather_Buf_1_1_SINGLE
     MODULE PROCEDURE PGSLib_Gather_Buf_1_1_DOUBLE
     MODULE PROCEDURE PGSLib_Gather_Buf_1_1_LOG
     MODULE PROCEDURE PGSLib_Gather_Buf_2_2_INT
     MODULE PROCEDURE PGSLib_Gather_Buf_2_2_SINGLE
     MODULE PROCEDURE PGSLib_Gather_Buf_2_2_DOUBLE
     MODULE PROCEDURE PGSLib_Gather_Buf_2_2_LOG
     MODULE PROCEDURE PGSLib_Gather_Buf_3_3_INT
     MODULE PROCEDURE PGSLib_Gather_Buf_3_3_SINGLE
     MODULE PROCEDURE PGSLib_Gather_Buf_3_3_DOUBLE
     MODULE PROCEDURE PGSLib_Gather_Buf_3_3_LOG
  END INTERFACE

  ! These are subroutines.  Provided because SGI
  ! has trouble with the functions which return arrays.
  INTERFACE PGSLib_Buffer_Gather
     MODULE PROCEDURE PGSLib_Buf_Gather_1_1_INT
     MODULE PROCEDURE PGSLib_Buf_Gather_1_1_SINGLE
     MODULE PROCEDURE PGSLib_Buf_Gather_1_1_DOUBLE
     MODULE PROCEDURE PGSLib_Buf_Gather_1_1_LOG
     MODULE PROCEDURE PGSLib_Buf_Gather_2_2_INT
     MODULE PROCEDURE PGSLib_Buf_Gather_2_2_SINGLE
     MODULE PROCEDURE PGSLib_Buf_Gather_2_2_DOUBLE
     MODULE PROCEDURE PGSLib_Buf_Gather_2_2_LOG
     MODULE PROCEDURE PGSLib_Buf_Gather_3_3_INT
     MODULE PROCEDURE PGSLib_Buf_Gather_3_3_SINGLE
     MODULE PROCEDURE PGSLib_Buf_Gather_3_3_DOUBLE
     MODULE PROCEDURE PGSLib_Buf_Gather_3_3_LOG
  END INTERFACE

  INTERFACE PGSLib_Scatter_Buffer
     MODULE PROCEDURE PGSLib_Scatter_Buf_1_1_INT
     MODULE PROCEDURE PGSLib_Scatter_Buf_1_1_SINGLE
     MODULE PROCEDURE PGSLib_Scatter_Buf_1_1_DOUBLE
     MODULE PROCEDURE PGSLib_Scatter_Buf_1_1_LOG
     MODULE PROCEDURE PGSLib_Scatter_Buf_2_2_INT
     MODULE PROCEDURE PGSLib_Scatter_Buf_2_2_SINGLE
     MODULE PROCEDURE PGSLib_Scatter_Buf_2_2_DOUBLE
     MODULE PROCEDURE PGSLib_Scatter_Buf_2_2_LOG
     MODULE PROCEDURE PGSLib_Scatter_Buf_3_3_INT
     MODULE PROCEDURE PGSLib_Scatter_Buf_3_3_SINGLE
     MODULE PROCEDURE PGSLib_Scatter_Buf_3_3_DOUBLE
     MODULE PROCEDURE PGSLib_Scatter_Buf_3_3_LOG
  END INTERFACE

CONTAINS

!======================================================================
! Some macros for cpp to use for defining the first array element
#define _DUP_1D_ARRAY_ELEMENT_ (LBOUND(Duplicate, 1))
#define _SUP_1D_ARRAY_ELEMENT_ (LBOUND(Supplement,1))

#define _DUP_2D_ARRAY_ELEMENT_ (1, LBOUND(Duplicate, 2))
#define _SUP_2D_ARRAY_ELEMENT_ (1, LBOUND(Supplement,2))

#define _DUP_3D_ARRAY_ELEMENT_ (1, 1, LBOUND(Duplicate, 3))
#define _SUP_3D_ARRAY_ELEMENT_ (1, 1, LBOUND(Supplement,3))
!======================================================================
!                    GATHER ROUTINES
!======================================================================

!======================================================================
!                    FUNCTIONS (PGSLib_Gather_Buffer)
!======================================================================

!========== Gather Vector from Scalar =================================

!======================================================================
!          PGSLib_Gather_Buf_1_1_INT(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_  PGSLib_Gather_Buf_1_1_INT
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_INT_C
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_     integer (PGSLib_Int_TYPE)
#define _BLOCKSIZE_     1
#include "gather_buff.fpp"

    
!======================================================================
!          PGSLib_Gather_Buf_1_1_Single(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_ PGSLIb_Gather_Buf_1_1_Single
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_FLOAT_C
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_    real (PGSLib_Single_TYPE)
#define _BLOCKSIZE_    1
#include "gather_buff.fpp"

!======================================================================
!          PGSLib_Gather_Buf_1_1_DOUBLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_ PGSLIb_Gather_Buf_1_1_Double
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_DOUBLE_C
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_    real (PGSLib_DOUBLE_TYPE)
#define _BLOCKSIZE_    1
#include "gather_buff.fpp"

!======================================================================
!          PGSLib_Gather_Buf_1_1_LOG(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_ PGSLIb_Gather_Buf_1_1_Log
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_LOG
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_    logical (PGSLib_Log_TYPE)
#define _BLOCKSIZE_    1
#include "gather_buff.fpp"

!========== Gather Vector from Vector =================================

!======================================================================
!          PGSLib_Gather_Buf_2_2_INT(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_2_2_Int
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "gather_buff.fpp"
    
!======================================================================
!          PGSLib_Gather_Buf_2_2_SINGLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_2_2_Single
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_FLOAT_C
#define _DATA_TYPE_     REAL (PGSLib_Single_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "gather_buff.fpp"

!======================================================================
!          PGSLib_Gather_Buf_2_2_DOUBLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_2_2_Double
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_DOUBLE_C
#define _DATA_TYPE_     REAL (PGSLib_Double_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "gather_buff.fpp"


!======================================================================
!          PGSLib_Gather_Buf_2_2_LOG(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_2_2_Log
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "gather_buff.fpp"

!========== Gather 3-d routines =================================

!======================================================================
!          PGSLib_Gather_Buf_3_3_INT(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_3_3_Int
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "gather_buff.fpp"
    
!======================================================================
!          PGSLib_Gather_Buf_3_3_Single(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_3_3_Single
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_FLOAT_C
#define _DATA_TYPE_     real (PGSLib_Single_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "gather_buff.fpp"

!======================================================================
!          PGSLib_Gather_Buf_3_3_DOUBLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_3_3_Double
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_DOUBLE_C
#define _DATA_TYPE_     real (PGSLib_Double_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "gather_buff.fpp"

!======================================================================
!          PGSLib_Gather_Buf_3_3_LOG(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Gather_Buf_3_3_Log
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "gather_buff.fpp"

!======================================================================
!                    SUBROUTINES (PGSLib_Buffer_Gather)
!======================================================================



!========== Gather Vector from Scalar =================================

!======================================================================
!          PGSLib_Buf_Gather_1_1_INT(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_  PGSLib_Buf_Gather_1_1_INT
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_INT_C
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_     integer (PGSLib_Int_TYPE)
#define _BLOCKSIZE_     1
#include "buff_gather.fpp"

    
!======================================================================
!          PGSLib_Buf_Gather_1_1_Single(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_ PGSLIb_Buf_Gather_1_1_Single
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_FLOAT_C
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_    real (PGSLib_Single_TYPE)
#define _BLOCKSIZE_    1
#include "buff_gather.fpp"

!======================================================================
!          PGSLib_Buf_Gather_1_1_DOUBLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_ PGSLIb_Buf_Gather_1_1_Double
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_DOUBLE_C
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_    real (PGSLib_DOUBLE_TYPE)
#define _BLOCKSIZE_    1
#include "buff_gather.fpp"

!======================================================================
!          PGSLib_Buf_Gather_1_1_LOG(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================

#define _ROUTINE_NAME_ PGSLIb_Buf_Gather_1_1_Log
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_LOG
#define _DUP_DIMENSION_ dimension(:) 
#define _SUP_DIMENSION_ dimension(Trace%N_Supplement) 
#define _DATA_TYPE_    logical (PGSLib_Log_TYPE)
#define _BLOCKSIZE_    1
#include "buff_gather.fpp"

!========== Gather Vector from Vector =================================

!======================================================================
!          PGSLib_Buf_Gather_2_2_INT(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_2_2_Int
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "buff_gather.fpp"
    
!======================================================================
!          PGSLib_Buf_Gather_2_2_SINGLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_2_2_Single
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_FLOAT_C
#define _DATA_TYPE_     REAL (PGSLib_Single_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "buff_gather.fpp"

!======================================================================
!          PGSLib_Buf_Gather_2_2_DOUBLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_2_2_Double
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_DOUBLE_C
#define _DATA_TYPE_     REAL (PGSLib_Double_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "buff_gather.fpp"


!======================================================================
!          PGSLib_Buf_Gather_2_2_LOG(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_2_2_Log
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_TYPE)
#define _DUP_DIMENSION_ dimension(:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)
#include "buff_gather.fpp"

!========== Gather 3-d routines =================================

!======================================================================
!          PGSLib_Buf_Gather_3_3_INT(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_3_3_Int
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "buff_gather.fpp"
    
!======================================================================
!          PGSLib_Buf_Gather_3_3_Single(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_3_3_Single
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_FLOAT_C
#define _DATA_TYPE_     real (PGSLib_Single_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "buff_gather.fpp"

!======================================================================
!          PGSLib_Buf_Gather_3_3_DOUBLE(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_3_3_Double
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_DOUBLE_C
#define _DATA_TYPE_     real (PGSLib_Double_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "buff_gather.fpp"

!======================================================================
!          PGSLib_Buf_Gather_3_3_LOG(Supplement, Duplicate, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a gather.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLIb_Buf_Gather_3_3_Log
#define _GATHER_BUF_C_  PGSLib_Gather_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _DUP_DIMENSION_ dimension(:,:,:) 
#define _SUP_DIMENSION_ dimension(SIZE(Duplicate,1),SIZE(Duplicate,2),Trace%N_Supplement) 
#define _BLOCKSIZE_     SIZE(Duplicate, 1)*SIZE(Duplicate, 2)
#include "buff_gather.fpp"

!======================================================================
!                    SCATTER ROUTINES
!======================================================================
    
!========== 1D Routines ====================

!======================================================================
!          PGSLib_Scatter_Buf_1_1_INT(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_1_1_Int
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _SUP_DIMENSION_ dimension(:)
#define _DUP_DIMENSION_ dimension(Trace%N_Duplicate)
#define _BLOCKSIZE_     1
#include "scatter_buff.fpp"

!======================================================================
!          PGSLib_Scatter_Buf_1_1_Single(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_1_1_Single
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_FLOAT_C
#define _DATA_TYPE_     Real (PGSLib_Single_Type)
#define _SUP_DIMENSION_ dimension(:)
#define _DUP_DIMENSION_ dimension(Trace%N_Duplicate)
#define _BLOCKSIZE_     1
#include "scatter_buff.fpp"

!======================================================================
!          PGSLib_Scatter_Buf_1_1_DOUBLE(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_1_1_DOUBLE
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_DOUBLE_C
#define _DATA_TYPE_     Real (PGSLib_Double_Type)
#define _SUP_DIMENSION_ dimension(:)
#define _DUP_DIMENSION_ dimension(Trace%N_Duplicate)
#define _BLOCKSIZE_     1
#include "scatter_buff.fpp"
    
!======================================================================
!          PGSLib_Scatter_Buf_1_1_LOG(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_1_1_LOG
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _SUP_DIMENSION_ dimension(:)
#define _DUP_DIMENSION_ dimension(Trace%N_Duplicate)
#define _BLOCKSIZE_     1
#include "scatter_buff.fpp"

!========== 2D Routines ====================

!======================================================================
!          PGSLib_Scatter_Buf_2_2_INT(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_2_2_INT
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _SUP_DIMENSION_ dimension(:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1)
#include "scatter_buff.fpp"

!======================================================================
!          PGSLib_Scatter_Buf_2_2_SINGLE(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_2_2_SINGLE
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_FLOAT_C
#define _DATA_TYPE_     real (PGSLib_Single_Type)
#define _SUP_DIMENSION_ dimension(:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1)
#include "scatter_buff.fpp"

!======================================================================
!          PGSLib_Scatter_Buf_2_2_DOUBLE(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_2_2_Double
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_DOUBLE_C
#define _DATA_TYPE_     real (PGSLib_Double_Type)
#define _SUP_DIMENSION_ dimension(:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1)
#include "scatter_buff.fpp"
    
!======================================================================
!          PGSLib_Scatter_Buf_2_2_LOG(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_2_2_Log
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _SUP_DIMENSION_ dimension(:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1)
#include "scatter_buff.fpp"

!========== 3D Routines ====================

!======================================================================
!          PGSLib_Scatter_Buf_3_3_INT(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_3_3_INT
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_INT_C
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _SUP_DIMENSION_ dimension(:,:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),SIZE(Supplement,2), Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1) * SIZE(Supplement,2)
#include "scatter_buff.fpp"

!======================================================================
!          PGSLib_Scatter_Buf_3_3_SINGLE(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_3_3_SINGLE
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_FLOAT_C
#define _DATA_TYPE_     real (PGSLib_Single_Type)
#define _SUP_DIMENSION_ dimension(:,:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),SIZE(Supplement,2), Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1) * SIZE(Supplement,2)
#include "scatter_buff.fpp"

!======================================================================
!          PGSLib_Scatter_Buf_3_3_DOUBLE(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_3_3_DOUBLE
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_DOUBLE_C
#define _DATA_TYPE_     real (PGSLib_Double_Type)
#define _SUP_DIMENSION_ dimension(:,:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),SIZE(Supplement,2), Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1) * SIZE(Supplement,2)
#include "scatter_buff.fpp"
    
!======================================================================
!          PGSLib_Scatter_Buf_3_3_LOG(Duplicate, Supplement, Trace)
! PURPOSE:
!          Execute the PE<->PE communication necessary for a Scatter.
!
!======================================================================
#define _ROUTINE_NAME_  PGSLib_Scatter_Buf_3_3_LOG
#define _SCATTER_BUF_C_ PGSLib_Scatter_Buf_LOG
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _SUP_DIMENSION_ dimension(:,:,:)
#define _DUP_DIMENSION_ dimension(SIZE(Supplement,1),SIZE(Supplement,2), Trace%N_Duplicate)
#define _BLOCKSIZE_     SIZE(Supplement,1) * SIZE(Supplement,2)
#include "scatter_buff.fpp"
    
END MODULE PGSLib_GS_Comm_MODULE
