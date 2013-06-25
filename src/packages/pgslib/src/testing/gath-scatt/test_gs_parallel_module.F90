MODULE Test_GS_Parallel_MODULE
  !======================================================================
  ! PURPOSE -
  !   Main module for gather/scatter test code.
  !======================================================================

  ! $Id: test_gs_parallel_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

  USE Test_Gather_MODULE,   ONLY: Test_EN_Gather
  USE Test_Scatter_MODULE,  ONLY: Test_EN_Scatter_SUM, Test_EN_Scatter_MAX, Test_EN_Scatter_OR
  USE Test_GS_Setup_MODULE, ONLY: Test_EN_GS_Setup
  USE Gather_Test_MODULE,   ONLY: Test_Gather
  USE Scatter_Test_MODULE,   ONLY: Test_Scatter_Sum, Test_Scatter_Max

END MODULE Test_GS_Parallel_MODULE
