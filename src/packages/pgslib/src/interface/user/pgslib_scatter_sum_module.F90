MODULE PGSLib_Scatter_SUM_MODULE
  !====================================================================
  ! PUPROSE
  !   Provide a user level interface to pgslib_OP?_scatter.  This
  !   module provides only the most common interfaces.  Others must
  !   be constructed as needed by the user.
  !====================================================================

  ! $Id: pgslib_scatter_sum_module.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $
  use pgslib_type_module
  use pgslib_globals_module
  use pgslib_gs_module
  use pgslib_timing_module
  use pgslib_utility_module
  implicit none

  PRIVATE
  PUBLIC :: PGSLib_Scatter_SUM

  INTERFACE PGSLib_Scatter_SUM
     MODULE PROCEDURE Sum_Scatter_1_1_INT
     MODULE PROCEDURE Sum_Scatter_1_1_Single
     MODULE PROCEDURE Sum_Scatter_1_1_Double
     MODULE PROCEDURE Sum_Scatter_1_2_INT
     MODULE PROCEDURE Sum_Scatter_1_2_Single
     MODULE PROCEDURE Sum_Scatter_1_2_Double
  END INTERFACE

CONTAINS

  !======================================================================
  !          Arithmetic scatter operations.
  ! PURPOSE:
  !          Scatter from 1D or 2D source to 1D destination
  !          Provides Integer, Single and Double versions
  !
  !======================================================================

  !======================================================================
  !          Sum Scatter operations.
  !======================================================================

#define _INCLUDE_FILE_ "scatter.fpp"
#define _OP_(a,b) a+ b
#define _OP_NAME_ Sum_Scatter
#define _GEN_OP_ID_(data_type,op_null) (op_null)
#define _DST_DIMENSION_ _NARROW_DIM_
#define _SRC_DIMENSION_ _FAT_DIM_
#include "arithmetic_gath-scatt.fpp"
#undef  _DST_DIMENSION_
#undef  _SRC_DIMENSION_

END MODULE PGSLib_Scatter_SUM_MODULE
