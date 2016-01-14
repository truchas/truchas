MODULE PGSLib_Scatter_MINMAX
  !====================================================================
  ! PUPROSE
  !   Provide a user level interface to pgslib_OP?_scatter.  This 
  !   module provides only the most common interfaces.  Others must
  !   be constructed as needed by the user.
  !====================================================================

  ! $Id: pgslib_scatter_minmax.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

  use pgslib_utility_module
  implicit none

  PRIVATE
  PUBLIC :: PGSLib_Scatter_Max
  PUBLIC :: PGSLib_Scatter_Min

  INTERFACE PGSLib_Scatter_MAX
     MODULE PROCEDURE Max_Scatter_1_1_INT
     MODULE PROCEDURE Max_Scatter_1_1_Single
     MODULE PROCEDURE Max_Scatter_1_1_Double
     MODULE PROCEDURE Max_Scatter_1_2_INT
     MODULE PROCEDURE Max_Scatter_1_2_Single
     MODULE PROCEDURE Max_Scatter_1_2_Double
  END INTERFACE

  INTERFACE PGSLib_Scatter_MIN
     MODULE PROCEDURE Min_Scatter_1_1_INT
     MODULE PROCEDURE Min_Scatter_1_1_Single
     MODULE PROCEDURE Min_Scatter_1_1_Double
     MODULE PROCEDURE Min_Scatter_1_2_INT
     MODULE PROCEDURE Min_Scatter_1_2_Single
     MODULE PROCEDURE Min_Scatter_1_2_Double
  END INTERFACE

CONTAINS

  !======================================================================
  !          Max Scatter operations.
  !======================================================================

#define _INCLUDE_FILE_ "scatter.fpp"
#define _OP_(a,b) MAX(a,b)
#define _OP_NAME_ Max_Scatter
#define _GEN_OP_ID_(data_type,op_null) MAXVAL([data_type::])
#define _DST_DIMENSION_ _NARROW_DIM_
#define _SRC_DIMENSION_ _FAT_DIM_
#include "arithmetic_gath-scatt.fpp"
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

  !======================================================================
  !          Min Scatter operations.
  !======================================================================

#define _INCLUDE_FILE_ "scatter.fpp"
#define _OP_(a,b) MIN(a,b)
#define _OP_NAME_ Min_Scatter
#define _GEN_OP_ID_(data_type,op_null) MINVAL([data_type::])
#define _DST_DIMENSION_ _NARROW_DIM_
#define _SRC_DIMENSION_ _FAT_DIM_
#include "arithmetic_gath-scatt.fpp"
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

END MODULE PGSLib_Scatter_MINMAX
