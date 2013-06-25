MODULE PGSLib_SCAN_SEG_BIT_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Provide PREFIX and SUFFIX operations, segment versions.
  ! These use CSHIFT, CSHIFT depends on the no-segment versions,
  ! and that is the reason these are split from the no-segment versions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_scan_seg_bit_module.F,v 1.1.1.1 2000/10/11 22:44:30 ferrell Exp $

  use pgslib_type_module
  use pgslib_globals_module
  use pgslib_utility_module, only : PGSLib_Fatal_ERROR, PGSLib_Scope_Check

  implicit none
  PRIVATE
  PUBLIC :: PGSLib_SEG_SUM_PREFIX, PGSLib_SEG_SUM_SUFFIX

  INTERFACE PGSLib_SEG_SUM_PREFIX
     MODULE PROCEDURE SUM_PREFIX_SEG_BIT_INT_1D
     MODULE PROCEDURE SUM_PREFIX_SEG_BIT_SINGLE_1D
     MODULE PROCEDURE SUM_PREFIX_SEG_BIT_DOUBLE_1D
  END INTERFACE
     
  INTERFACE PGSLib_SEG_SUM_SUFFIX
     MODULE PROCEDURE SUM_SUFFIX_SEG_BIT_INT_1D
     MODULE PROCEDURE SUM_SUFFIX_SEG_BIT_SINGLE_1D
     MODULE PROCEDURE SUM_SUFFIX_SEG_BIT_DOUBLE_1D
  END INTERFACE
     
  ! Interface blocks for internal support routines
  INTERFACE On_Node_SUM_PREFIX
     MODULE PROCEDURE On_Node_Prfx_SEG_BIT_INT_1D
     MODULE PROCEDURE On_Node_Prfx_SEG_BIT_SINGLE_1D
     MODULE PROCEDURE On_Node_Prfx_SEG_BIT_DOUBLE_1D
  END INTERFACE

  INTERFACE On_Node_SUM_SUFFIX
     MODULE PROCEDURE On_Node_Sufx_SEG_BIT_INT_1D
     MODULE PROCEDURE On_Node_Sufx_SEG_BIT_SINGLE_1D
     MODULE PROCEDURE On_Node_Sufx_SEG_BIT_DOUBLE_1D
  END INTERFACE

  
  ! Interface blocks for internal support routines
  INTERFACE SUM_PREFIX_Fixup
     MODULE PROCEDURE PREFIX_Fixup_SEG_BIT_INT_1D
     MODULE PROCEDURE PREFIX_Fixup_SEG_BIT_SINGLE_1D
     MODULE PROCEDURE PREFIX_Fixup_SEG_BIT_DOUBLE_1D
  END INTERFACE

  INTERFACE SUM_SUFFIX_Fixup
     MODULE PROCEDURE SUFFIX_Fixup_SEG_BIT_INT_1D
     MODULE PROCEDURE SUFFIX_Fixup_SEG_BIT_SINGLE_1D
     MODULE PROCEDURE SUFFIX_Fixup_SEG_BIT_DOUBLE_1D
  END INTERFACE

  ! Interface blocks for the C routines which call MPI
  INTERFACE Off_Node_SUM_PREFIX
     subroutine Off_Node_SUM_PREFIX_INT_C(Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       use pgslib_type_module
       implicit none
       integer(PGSLib_Int_Type),  intent (INOUT) :: Dest_Data
       integer,                   intent (INOUT) :: Dest_Seg
       integer(PGSLib_Int_Type),  intent (INOUT) :: Src_Data
       integer,                   intent (INOUT) :: Src_Seg
     end subroutine Off_Node_SUM_PREFIX_INT_C
     subroutine Off_Node_SUM_PREFIX_FLOAT_C(Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       use pgslib_type_module
       implicit none
       REAL(PGSLib_Real_Type),    intent (INOUT) :: Dest_Data
       integer,                   intent (INOUT) :: Dest_Seg
       REAL(PGSLib_Real_Type),    intent (INOUT) :: Src_Data
       integer,                   intent (INOUT) :: Src_Seg
     end subroutine Off_Node_SUM_PREFIX_FLOAT_C
     subroutine Off_Node_SUM_PREFIX_DOUBLE_C(Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       use PGSLib_Type_MODULE
       implicit none
       REAL (PGSLib_DOUBLE_Type), intent (INOUT) :: Dest_Data
       integer,                   intent (INOUT) :: Dest_Seg
       REAL (PGSLib_DOUBLE_Type), intent (INOUT) :: Src_Data
       integer,                   intent (INOUT) :: Src_Seg
     end subroutine Off_Node_SUM_PREFIX_DOUBLE_C
  END INTERFACE
  
  INTERFACE Off_Node_SUM_SUFFIX
     subroutine Off_Node_SUM_SUFFIX_INT_C(Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       use pgslib_type_module
       implicit none
       integer(PGSLib_Int_Type),  intent (INOUT) :: Dest_Data
       integer,                   intent (INOUT) :: Dest_Seg
       integer(PGSLib_Int_Type),  intent (INOUT) :: Src_Data
       integer,                   intent (INOUT) :: Src_Seg
     end subroutine Off_Node_SUM_SUFFIX_INT_C
     subroutine Off_Node_SUM_SUFFIX_FLOAT_C(Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       use pgslib_type_module
       implicit none
       REAL(PGSLib_Real_Type),    intent (INOUT) :: Dest_Data
       integer,                   intent (INOUT) :: Dest_Seg
       REAL(PGSLib_Real_Type),    intent (INOUT) :: Src_Data
       integer,                   intent (INOUT) :: Src_Seg
     end subroutine Off_Node_SUM_SUFFIX_FLOAT_C
     subroutine Off_Node_SUM_SUFFIX_DOUBLE_C(Dest_Data, Dest_Seg, Src_Data, Src_Seg)
       use PGSLib_Type_MODULE
       implicit none
       REAL (PGSLib_DOUBLE_Type), intent (INOUT) :: Dest_Data
       integer,                   intent (INOUT) :: Dest_Seg
       REAL (PGSLib_DOUBLE_Type), intent (INOUT) :: Src_Data
       integer,                   intent (INOUT) :: Src_Seg
     end subroutine Off_Node_SUM_SUFFIX_DOUBLE_C
  END INTERFACE
  
  
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUM_PREFIX Top_level routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _ROUTINE_NAME_   SUM_PREFIX_SEG_BIT_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _ON_NODE_SCAN_   On_Node_SUM_PREFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_PREFIX
#define _SCAN_FIXUP_     SUM_PREFIX_Fixup
#define _LAST_           Local_N

#include "scan-seg-bit-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_PREFIX_SEG_BIT_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _ON_NODE_SCAN_   On_Node_SUM_PREFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_PREFIX
#define _SCAN_FIXUP_     SUM_PREFIX_Fixup
#define _LAST_           Local_N

#include "scan-seg-bit-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_PREFIX_SEG_BIT_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _ON_NODE_SCAN_   On_Node_SUM_PREFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_PREFIX
#define _SCAN_FIXUP_     SUM_PREFIX_Fixup
#define _LAST_           Local_N

#include "scan-seg-bit-1d.fpp"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUM_SUFFIX Top_level routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   SUM_SUFFIX_SEG_BIT_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _ON_NODE_SCAN_   On_Node_SUM_SUFFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_SUFFIX
#define _SCAN_FIXUP_     SUM_SUFFIX_Fixup
#define _LAST_           1

#include "scan-seg-bit-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_SUFFIX_SEG_BIT_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _ON_NODE_SCAN_   On_Node_SUM_SUFFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_SUFFIX
#define _SCAN_FIXUP_     SUM_SUFFIX_Fixup
#define _LAST_           1

#include "scan-seg-bit-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_SUFFIX_SEG_BIT_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _ON_NODE_SCAN_   On_Node_SUM_SUFFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_SUFFIX
#define _SCAN_FIXUP_     SUM_SUFFIX_Fixup
#define _LAST_           1

#include "scan-seg-bit-1d.fpp"
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! On node PREFIX routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Segmented versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _ROUTINE_NAME_   On_Node_Prfx_SEG_BIT_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "on-node-scan-seg-bit-1d.fpp"
       
  !======================================================================

#define  _ROUTINE_NAME_ On_Node_Prfx_SEG_BIT_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "on-node-scan-seg-bit-1d.fpp"
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   On_Node_Prfx_SEG_BIT_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "on-node-scan-seg-bit-1d.fpp"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! On node SUFFIX routines
  !

#define _ROUTINE_NAME_   On_Node_Sufx_SEG_BIT_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "on-node-scan-seg-bit-1d.fpp"
       
  !======================================================================

#define  _ROUTINE_NAME_ On_Node_Sufx_SEG_BIT_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "on-node-scan-seg-bit-1d.fpp"
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   On_Node_Sufx_SEG_BIT_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "on-node-scan-seg-bit-1d.fpp"
  
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PREFIX Fixup routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Segmented versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _ROUTINE_NAME_   Prefix_Fixup_SEG_BIT_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "scan-fixup-seg-bit-1d.fpp"
       
  !======================================================================

#define  _ROUTINE_NAME_ Prefix_Fixup_SEG_BIT_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "scan-fixup-seg-bit-1d.fpp"
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   Prefix_Fixup_SEG_BIT_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "scan-fixup-seg-bit-1d.fpp"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! On node SUFFIX routines
  !

#define _ROUTINE_NAME_   Suffix_Fixup_SEG_BIT_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "scan-fixup-seg-bit-1d.fpp"
       
  !======================================================================

#define  _ROUTINE_NAME_ Suffix_Fixup_SEG_BIT_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "scan-fixup-seg-bit-1d.fpp"
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   Suffix_Fixup_SEG_BIT_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "scan-fixup-seg-bit-1d.fpp"
  
    
END MODULE PGSLib_SCAN_SEG_BIT_MODULE
