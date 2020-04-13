MODULE PGSLib_Misc_MODULE

  ! $Id: pgslib_misc_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $
  USE PGSLIB_PERMUTE_MODULE,  ONLY : PGSLIB_PERMUTE, PGSLib_Redistribute, PGSLib_Global_PAck
  USE PGSLib_Shift_MODULE,    ONLY : PGSLib_Global_CShift, PGSLib_Global_EOShift
  USE PGSLib_Decomp_Module,   ONLY : PGSLib_Block_Decompose, PGSLib_Block_Size

  PRIVATE
  PUBLIC :: PGSLib_Permute, PGSLib_Redistribute, PGSLib_Global_Pack
  PUBLIC :: PGSLib_Global_CShift, PGSLib_Global_EOSHIFT
  PUBLIC :: PGSLib_Block_Decompose, PGSLib_Block_Size

END MODULE PGSLib_Misc_MODULE
