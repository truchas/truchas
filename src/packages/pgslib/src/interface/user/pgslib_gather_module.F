MODULE PGSLib_Gather_MODULE
  !====================================================================
  ! PUPROSE
  !   Provide a user level interface to pgslib_gather.  This 
  !   module provides only the most common interfaces.  Others must
  !   be constructed as needed by the user.
  !====================================================================

  ! $Id: pgslib_gather_module.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

  use pgslib_utility_module
  implicit none

  PRIVATE
  PUBLIC :: PGSLib_Gather

  INTERFACE PGSLib_Gather
     MODULE PROCEDURE Gather_1_1_Int
     MODULE PROCEDURE Gather_1_1_Single
     MODULE PROCEDURE Gather_1_1_Double
     MODULE PROCEDURE Gather_1_1_Log
     MODULE PROCEDURE Gather_1_2_Int
     MODULE PROCEDURE Gather_1_2_Single
     MODULE PROCEDURE Gather_1_2_Double
     MODULE PROCEDURE Gather_1_2_Log
  END INTERFACE

CONTAINS

  !======================================================================
  !          Gather(Dest, Src, Index, Mask, Trace)
  ! PURPOSE:
  !          Gather from src to dest.
  !          This include part builds the specific routines:
  !          1_1_(Int, Single, Double)
  !          1_2_(Int, Single, Double)
  !
  !======================================================================
#define _INCLUDE_FILE_ "gather.fpp"
#define _OP_ 
#define _OP_NAME_       Gather
#define _GEN_OP_ID_(op_null) op_null
#define _DST_DIMENSION_ _FAT_DIM_
#define _SRC_DIMENSION_ _NARROW_DIM_
#include "arithmetic_gath-scatt.fpp"
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

  !======================================================================
  !          Gather_1_1_Log
  !======================================================================
#define _ROUTINE_NAME_  Gather_1_1_Log
#define _OP_(a,b)       
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _OP_ID_         .TRUE.
#define _DST_DIMENSION_ dimension(:)
#define _SRC_DIMENSION_ dimension(:)
#define _INDEX_SIZE_    SIZE(Index,1)
#define _ARRAY_ELEMENT_ i
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    1
#include "gather.fpp"
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_
#undef  _OP_

  !======================================================================
  !          Gather_1_2_Log
  !======================================================================
#define _ROUTINE_NAME_  Gather_1_2_Log
#define _OP_(a,b)       
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _OP_ID_         .TRUE.
#define _DST_DIMENSION_ dimension(:,:)
#define _SRC_DIMENSION_ dimension(:)
#define _INDEX_SIZE_      SIZE(Index,1), SIZE(Index,2)
#define _ARRAY_ELEMENT_ i,j
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    SIZE(Local_Index,2)
#include "gather.fpp"
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_
#undef  _OP_

END MODULE PGSLib_Gather_MODULE
