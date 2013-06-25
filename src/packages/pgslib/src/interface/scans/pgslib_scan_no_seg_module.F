MODULE PGSLib_SCAN_NO_SEG_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Provide PREFIX and SUFFIX operations, no-segment versions
  ! These do not use CSHIFT, so they may be USEd by the SHIFT routines.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_scan_no_seg_module.F,v 1.1.1.1 2000/10/11 22:44:30 ferrell Exp $

  use pgslib_type_module
  use pgslib_globals_module
  use pgslib_utility_module, only : PGSLib_Fatal_ERROR, PGSLib_Scope_Check

  implicit none
  PRIVATE
  PUBLIC :: PGSLib_SUM_PREFIX, PGSLib_SUM_SUFFIX
  PUBLIC :: PGSLib_PARITY_PREFIX, PGSLib_PARITY_SUFFIX

  INTERFACE PGSLib_SUM_PREFIX
     MODULE PROCEDURE SUM_PREFIX_NO_SEG_INT_1D
     MODULE PROCEDURE SUM_PREFIX_NO_SEG_SINGLE_1D
     MODULE PROCEDURE SUM_PREFIX_NO_SEG_DOUBLE_1D
  END INTERFACE
     
  INTERFACE PGSLib_SUM_SUFFIX
     MODULE PROCEDURE SUM_SUFFIX_NO_SEG_INT_1D
     MODULE PROCEDURE SUM_SUFFIX_NO_SEG_SINGLE_1D
     MODULE PROCEDURE SUM_SUFFIX_NO_SEG_DOUBLE_1D
  END INTERFACE
     
  INTERFACE PGSLib_PARITY_PREFIX
     MODULE PROCEDURE PARITY_PREFIX_1D
  END INTERFACE

  INTERFACE PGSLib_PARITY_SUFFIX
     MODULE PROCEDURE PARITY_SUFFIX_1D
  END INTERFACE

  ! Interface blocks for internal support routines
  INTERFACE On_Node_SUM_PREFIX
     MODULE PROCEDURE On_Node_PREFIX_NO_SEG_INT_1D
     MODULE PROCEDURE On_Node_PREFIX_NO_SEG_SINGLE_1D
     MODULE PROCEDURE On_Node_PREFIX_NO_SEG_DOUBLE_1D
  END INTERFACE

  INTERFACE On_Node_SUM_SUFFIX
     MODULE PROCEDURE On_Node_SUFFIX_NO_SEG_INT_1D
     MODULE PROCEDURE On_Node_SUFFIX_NO_SEG_SINGLE_1D
     MODULE PROCEDURE On_Node_SUFFIX_NO_SEG_DOUBLE_1D
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


#define _ROUTINE_NAME_   SUM_PREFIX_NO_SEG_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _ON_NODE_SCAN_   On_Node_SUM_PREFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_PREFIX
#define _LAST_           Local_N

#include "scan-no-seg-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_PREFIX_NO_SEG_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _ON_NODE_SCAN_   On_Node_SUM_PREFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_PREFIX
#define _LAST_           Local_N

#include "scan-no-seg-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_PREFIX_NO_SEG_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _ON_NODE_SCAN_   On_Node_SUM_PREFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_PREFIX
#define _LAST_           Local_N

#include "scan-no-seg-1d.fpp"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUM_SUFFIX Top_level routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   SUM_SUFFIX_NO_SEG_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _ON_NODE_SCAN_   On_Node_SUM_SUFFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_SUFFIX
#define _LAST_           1

#include "scan-no-seg-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_SUFFIX_NO_SEG_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _ON_NODE_SCAN_   On_Node_SUM_SUFFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_SUFFIX
#define _LAST_           1

#include "scan-no-seg-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_   SUM_SUFFIX_NO_SEG_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _ON_NODE_SCAN_   On_Node_SUM_SUFFIX
#define _OFF_NODE_SCAN_  Off_Node_SUM_SUFFIX
#define _LAST_           1

#include "scan-no-seg-1d.fpp"
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! On node PREFIX routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! No segment versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _ROUTINE_NAME_   On_Node_PREFIX_NO_SEG_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "on-node-scan-no-seg-1d.fpp"
       
  !======================================================================

#define  _ROUTINE_NAME_ On_Node_PREFIX_NO_SEG_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "on-node-scan-no-seg-1d.fpp"
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   On_Node_PREFIX_NO_SEG_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _FIRST_          1
#define _LAST_           Local_N
#define _INDEX_INCREMENT_ 1

#include "on-node-scan-no-seg-1d.fpp"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! On node SUFFIX routines
  !

#define _ROUTINE_NAME_   On_Node_SUFFIX_NO_SEG_INT_1D
#define _SCAN_DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_          0
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "on-node-scan-no-seg-1d.fpp"
       
  !======================================================================

#define  _ROUTINE_NAME_ On_Node_SUFFIX_NO_SEG_SINGLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_          0.0
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "on-node-scan-no-seg-1d.fpp"
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   On_Node_SUFFIX_NO_SEG_DOUBLE_1D
#define _SCAN_DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_          REAL(0.0,KIND=PGSLib_Double_Type)
#define _FIRST_          Local_N
#define _LAST_           1
#define _INDEX_INCREMENT_ (-1)

#include "on-node-scan-no-seg-1d.fpp"
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PARITY_PREFIX_1D(ARRAY, DIM)
    implicit none
    ! Subroutine arguments
    logical (PGSLib_Log_Type),        &
         &   intent(IN   ),           &
         &   TARGET,                  &
         &   dimension(            :) :: ARRAY
    integer(PGSLib_INT_Type),         &
         &   intent(IN   )          , &
         &   optional                 :: DIM
    ! Function type
    logical (PGSLib_Log_TYPE),        &
         &   dimension(SIZE(ARRAY,1)) :: PARITY_PREFIX_1D
    

    ! Local variables
    
    integer (PGSLib_INT_TYPE)      ,  &
         &   dimension(SIZE(ARRAY,1)) :: Temp_Source, Temp_Dest

    ! Check that we got valid argument combinations
    if (PRESENT(DIM)) then
       if (DIM /= 1) then
          call PGSLib_Fatal_Error('In PARITY_PREFIX if DIM is present it must == 1')
       end if
    end if

    ! Note that parity is XOR, so if we put
    ! T = 1, F = 0, then use + MOD 2, that is equivalent

    Temp_Source = MERGE(1,0,ARRAY)
    Temp_Dest = PGSLib_SUM_PREFIX(Temp_Source)
    PARITY_PREFIX_1D = MOD(Temp_Dest, 2) == 1

    return
  end function PARITY_PREFIX_1D
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PARITY_SUFFIX_1D(ARRAY, DIM)
    implicit none
    ! Subroutine arguments
    logical (PGSLib_Log_Type),        &
         &   intent(IN   ),           &
         &   TARGET,                  &
         &   dimension(            :) :: ARRAY
    integer(PGSLib_INT_Type),         &
         &   intent(IN   )          , &
         &   optional                 :: DIM
    ! Function type
    logical (PGSLib_Log_TYPE),        &
         &   dimension(SIZE(ARRAY,1)) :: PARITY_SUFFIX_1D
    

    ! Local variables
    
    integer (PGSLib_INT_TYPE)      ,  &
         &   dimension(SIZE(ARRAY,1)) :: Temp_Source, Temp_Dest

    ! Check that we got valid argument combinations
    if (PRESENT(DIM)) then
       if (DIM /= 1) then
          call PGSLib_Fatal_ERROR('In PARITY_SUFFIX if DIM is present it must == 1')
       end if
    end if

    ! Note that parity is XOR, so if we put
    ! T = 1, F = 0, then use + MOD 2, that is equivalent

    Temp_Source = MERGE(1,0,ARRAY)
    Temp_Dest = PGSLib_SUM_SUFFIX(Temp_Source)
    PARITY_SUFFIX_1D = MOD(Temp_Dest, 2) == 1

    return
  end function PARITY_SUFFIX_1D
    

    
END MODULE PGSLib_SCAN_NO_SEG_MODULE
