MODULE PGSLib_Scatter_LOG
  !====================================================================
  ! PUPROSE
  !   Provide a user level interface to pgslib_OP?_scatter.  This
  !   module provides only the most common interfaces.  Others must
  !   be constructed as needed by the user.
  !====================================================================

  ! $Id: pgslib_scatter_log.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $
  use pgslib_type_module
  use pgslib_globals_module
  use pgslib_gs_module
  use pgslib_timing_module
  use pgslib_utility_module
  implicit none

  PRIVATE
  PUBLIC :: PGSLib_Scatter_And
  PUBLIC :: PGSLib_Scatter_Or

  INTERFACE PGSLib_Scatter_AND
     MODULE PROCEDURE And_Scatter_1_1_Log
     MODULE PROCEDURE And_Scatter_1_2_Log
  END INTERFACE

  INTERFACE PGSLib_Scatter_OR
     MODULE PROCEDURE Or_Scatter_1_1_Log
     MODULE PROCEDURE Or_Scatter_1_2_Log
  END INTERFACE
CONTAINS

  !======================================================================
  !          Logical scatter operations.
  ! PURPOSE:
  !          Scatter from 1D or 2D source to 1D destination
  !
  !======================================================================

  !======================================================================
  !          And_Scatter_1_1_Log
  !======================================================================
#define _ROUTINE_NAME_  And_Scatter_1_1_Log
#define _OP_(a,b)       a .AND. b
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _OP_ID_         .TRUE.
#define _DST_DIMENSION_ dimension(:)
#define _SRC_DIMENSION_ dimension(:)
#define _INDEX_SIZE_      SIZE(Index,1)
#define _ARRAY_ELEMENT_ i
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    1
#include "scatter.fpp"
#undef  _OP_
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

  !======================================================================
  !          And_Scatter_1_2_Log
  !======================================================================
#define _ROUTINE_NAME_  And_Scatter_1_2_Log
#define _OP_(a,b)       a .AND. b
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _OP_ID_         .TRUE.
#define _DST_DIMENSION_ dimension(:)
#define _SRC_DIMENSION_ dimension(:,:)
#define _INDEX_SIZE_      SIZE(Index,1), SIZE(Index,2)
#define _ARRAY_ELEMENT_ i,j
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    SIZE(Local_Index,2)
#include "scatter.fpp"
#undef  _OP_
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

  !======================================================================
  !          Or_Scatter_1_1_Log
  !======================================================================
#define _ROUTINE_NAME_  Or_Scatter_1_1_Log
#define _OP_(a,b)       a .OR. b
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _OP_ID_         .FALSE.
#define _DST_DIMENSION_ dimension(:)
#define _SRC_DIMENSION_ dimension(:)
#define _INDEX_SIZE_      SIZE(Index,1)
#define _ARRAY_ELEMENT_ i
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    1
#include "scatter.fpp"
#undef  _OP_
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

  !======================================================================
  !          Or_Scatter_1_2_Log
  !======================================================================
#define _ROUTINE_NAME_  Or_Scatter_1_2_Log
#define _OP_(a,b)       a .OR. b
#define _DATA_TYPE_     logical (PGSLib_Log_Type)
#define _OP_ID_         .FALSE.
#define _DST_DIMENSION_ dimension(:)
#define _SRC_DIMENSION_ dimension(:,:)
#define _INDEX_SIZE_      SIZE(Index,1), SIZE(Index,2)
#define _ARRAY_ELEMENT_ i,j
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    SIZE(Local_Index,2)
#include "scatter.fpp"
#undef  _OP_
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

END MODULE PGSLib_Scatter_LOG
