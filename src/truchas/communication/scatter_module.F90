MODULE SCATTER_MODULE
  !=======================================================================
  ! PURPOSE - 
  !   Routines for scatter operations, used by both 
  !   Element<->Element and Element<->Node.
  !=======================================================================

  use truchas_logging_services
  use kind_module
  use gs_info_module
  use mesh_module,  only: MESH_CONNECTIVITY, CllNgbr, Vrtx
  use pgslib_module,only: PGSLib_GS_Trace,      &
                          PGSLib_Size_Of_Dup,   &
                          PGSLib_Size_Of_Sup,   &
                          PGSLib_Dup_Index,     &
                          PGSLib_Scatter_Buffer

  implicit none
  save
!  private
  public :: Sum_Scatter, Min_Scatter, Max_Scatter, Or_Scatter


  ! Interface blocks
  INTERFACE Sum_Scatter
     MODULE PROCEDURE Sum_Scatter_Scalar_INT
     MODULE PROCEDURE Sum_Scatter_Scalar_Single
     MODULE PROCEDURE Sum_Scatter_Scalar_Double
     MODULE PROCEDURE Sum_Scatter_Vector_INT
     MODULE PROCEDURE Sum_Scatter_Vector_Single
     MODULE PROCEDURE Sum_Scatter_Vector_Double
  END INTERFACE

  INTERFACE Min_Scatter
     MODULE PROCEDURE Min_Scatter_Scalar_INT
     MODULE PROCEDURE Min_Scatter_Scalar_Single
     MODULE PROCEDURE Min_Scatter_Scalar_Double
     MODULE PROCEDURE Min_Scatter_Vector_INT
     MODULE PROCEDURE Min_Scatter_Vector_Single
     MODULE PROCEDURE Min_Scatter_Vector_Double
  END INTERFACE

  INTERFACE Max_Scatter
     MODULE PROCEDURE Max_Scatter_Scalar_INT
     MODULE PROCEDURE Max_Scatter_Scalar_Single
     MODULE PROCEDURE Max_Scatter_Scalar_Double
     MODULE PROCEDURE Max_Scatter_Vector_INT
     MODULE PROCEDURE Max_Scatter_Vector_Single
     MODULE PROCEDURE Max_Scatter_Vector_Double
  END INTERFACE

  INTERFACE Or_Scatter
     MODULE PROCEDURE Or_Scatter_Scalar_Log
     MODULE PROCEDURE Or_Scatter_Vector_Log
  END INTERFACE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><>   SUM_SCATTER_SCALAR_INT       <><><><><><><><><><><><><>

#define _ROUTINE_NAME_ SUM_SCATTER_SCALAR_INT
#define _DATA_TYPE_    integer (int_kind)
#define _OP_ID_        INT(0,   int_kind)
#define _PLUS_

#include "scatter_infix_parallel_op_s_s_include.fpp"

#undef _PLUS_

  ! <><><><><><>   SUM_SCATTER_SCALAR_SINGLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ SUM_SCATTER_SCALAR_SINGLE
#define _DATA_TYPE_    real (single_kind)
#define _OP_ID_        REAL(0.0,   single_kind)
#define _PLUS_

#include "scatter_infix_parallel_op_s_s_include.fpp"

#undef _PLUS_

  ! <><><><><><>   SUM_SCATTER_SCALAR_DOUBLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ SUM_SCATTER_SCALAR_DOUBLE
#define _DATA_TYPE_    real (double_kind)
#define _OP_ID_        REAL(0.0,   double_kind)
#define _PLUS_

#include "scatter_infix_parallel_op_s_s_include.fpp"

#undef _PLUS_

  ! <><><><><><>    SUM_SCATTER_VECTOR_INT    <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ SUM_SCATTER_VECTOR_INT
#define _DATA_TYPE_    integer (int_kind)
#define _OP_ID_        INT(0,   int_kind)
#define _PLUS_

#include "scatter_infix_parallel_op_s_v_include.fpp"

#undef _PLUS_

  ! <><><><><><>    SUM_SCATTER_VECTOR_SINGLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ SUM_SCATTER_VECTOR_SINGLE
#define _DATA_TYPE_    real (single_kind)
#define _OP_ID_        REAL(0.0,   single_kind)
#define _PLUS_

#include "scatter_infix_parallel_op_s_v_include.fpp"

#undef _PLUS_

  ! <><><><><><>   SUM_SCATTER_VECTOR_DOUBLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ SUM_SCATTER_VECTOR_DOUBLE
#define _DATA_TYPE_    real (double_kind)
#define _OP_ID_        REAL(0.0,   double_kind)
#define _PLUS_

#include "scatter_infix_parallel_op_s_v_include.fpp"

#undef _PLUS_

  ! <><><><><><>   MIN_SCATTER_SCALAR_INT     <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MIN_SCATTER_SCALAR_INT
#define _DATA_TYPE_    integer (int_kind)
#define _OP_ID_        MINVAL(Src,   MASK=.FALSE.)
#define _MIN_

#include "scatter_prefix_parallel_op_s_s_include.fpp"

#undef _MIN_

  ! <><><><><><>   MIN_SCATTER_SCALAR_SINGLE     <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MIN_SCATTER_SCALAR_SINGLE
#define _DATA_TYPE_    real (single_kind)
#define _OP_ID_        MINVAL(Src, MASK=.FALSE.)
#define _MIN_

#include "scatter_prefix_parallel_op_s_s_include.fpp"

#undef _MIN_

  ! <><><><><><>   MIN_SCATTER_SCALAR_DOUBLE     <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MIN_SCATTER_SCALAR_DOUBLE
#define _DATA_TYPE_    real (double_kind)
#define _OP_ID_        MINVAL(Src, MASK=.FALSE.)
#define _MIN_

#include "scatter_prefix_parallel_op_s_s_include.fpp"

#undef _MIN_

  ! <><><><><><>   MIN_SCATTER_VECTOR_INT     <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MIN_SCATTER_VECTOR_INT
#define _DATA_TYPE_    integer (int_kind)
#define _OP_ID_        MINVAL(Src,   MASK=.FALSE.)
#define _MIN_

#include "scatter_prefix_parallel_op_s_v_include.fpp"

#undef _MIN_

  ! <><><><><><>   MIN_SCATTER_VECTOR_SINGLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MIN_SCATTER_VECTOR_SINGLE
#define _DATA_TYPE_  real (single_kind)
#define _OP_ID_      MINVAL(Src, MASK=.FALSE.)
#define _MIN_

#include "scatter_prefix_parallel_op_s_v_include.fpp"

#undef _MIN_

  ! <><><><><><>   MIN_SCATTER_VECTOR_DOUBLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MIN_SCATTER_VECTOR_DOUBLE
#define _DATA_TYPE_  real (double_kind)
#define _OP_ID_      MINVAL(Src, MASK=.FALSE.)
#define _MIN_

#include "scatter_prefix_parallel_op_s_v_include.fpp"

#undef _MIN_

  ! <><><><><><>   MAX_SCATTER_SCALAR_INT     <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MAX_SCATTER_SCALAR_INT
#define _DATA_TYPE_  integer (int_kind)
#define _OP_ID_      MAXVAL(Src, MASK=.FALSE.)
#define _MAX_

#include "scatter_prefix_parallel_op_s_s_include.fpp"

#undef _MAX_

  ! <><><><><><>   MAX_SCATTER_SCALAR_SINGLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MAX_SCATTER_SCALAR_SINGLE
#define _DATA_TYPE_  real (single_kind)
#define _OP_ID_      MAXVAL(Src, MASK=.FALSE.)
#define _MAX_

#include "scatter_prefix_parallel_op_s_s_include.fpp"

#undef _MAX_

  ! <><><><><><>   MAX_SCATTER_SCALAR_DOUBLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MAX_SCATTER_SCALAR_DOUBLE
#define _DATA_TYPE_  real (double_kind)
#define _OP_ID_      MAXVAL(Src, MASK=.FALSE.)
#define _MAX_


#include "scatter_prefix_parallel_op_s_s_include.fpp"

#undef _MAX_

  ! <><><><><><>   MAX_SCATTER_VECTOR_INT     <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MAX_SCATTER_VECTOR_INT
#define _DATA_TYPE_  integer (int_kind)
#define _OP_ID_      MAXVAL(Src, MASK=.FALSE.)
#define _MAX_

#include "scatter_prefix_parallel_op_s_v_include.fpp"

#undef _MAX_

  ! <><><><><><>   MAX_SCATTER_VECTOR_SINGLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MAX_SCATTER_VECTOR_SINGLE
#define _DATA_TYPE_  real (single_kind)
#define _OP_ID_      MAXVAL(Src, MASK=.FALSE.)
#define _MAX_

#include "scatter_prefix_parallel_op_s_v_include.fpp"

#undef _MAX_
       
  ! <><><><><><>   MAX_SCATTER_VECTOR_DOUBLE   <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ MAX_SCATTER_VECTOR_DOUBLE
#define _DATA_TYPE_  real (double_kind)
#define _OP_ID_      MAXVAL(Src, MASK=.FALSE.)
#define _MAX_

#include "scatter_prefix_parallel_op_s_v_include.fpp"

#undef _MAX_

  ! <><><><><><>   OR_SCATTER_SCALAR_LOG      <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ OR_SCATTER_SCALAR_LOG
#define _DATA_TYPE_  logical (log_kind)
#define _OP_ID_      .FALSE.
#define _OR_

#include "scatter_infix_parallel_op_s_s_include.fpp"

#undef _OR_

  ! <><><><><><>   OR_SCATTER_VECTOR_LOG      <><><><><><><><><><><><><><><>

#define _ROUTINE_NAME_ OR_SCATTER_VECTOR_LOG
#define _DATA_TYPE_  logical (log_kind)
#define _OP_ID_      .FALSE.
#define _OR_

#include "scatter_infix_parallel_op_s_v_include.fpp"

#undef _OR_

END MODULE SCATTER_MODULE
