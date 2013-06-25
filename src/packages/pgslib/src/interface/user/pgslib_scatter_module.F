MODULE PGSLib_Scatter_MODULE
  !====================================================================
  ! PUPROSE
  !   Provide a user level interface to pgslib_OP?_scatter.  This 
  !   module provides only the most common interfaces.  Others must
  !   be constructed as needed by the user.
  !====================================================================

  ! $Id: pgslib_scatter_module.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

  use pgslib_scatter_sum_module
  use pgslib_scatter_minmax
  use pgslib_scatter_log
  implicit none

  PRIVATE
  PUBLIC :: PGSLib_Scatter_SUM
  PUBLIC :: PGSLib_Scatter_Max
  PUBLIC :: PGSLib_Scatter_Min
  PUBLIC :: PGSLib_Scatter_And
  PUBLIC :: PGSLib_Scatter_Or

END MODULE PGSLib_Scatter_MODULE
