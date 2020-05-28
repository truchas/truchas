MODULE PGSLib_Shift_MODULE
  !======================================================================
  !  PURPOSE -
  !    Provide CSHIFT like capability
  !    Given Dest(N), Src(N), SHIFT, with Dest&Src global arrays,
  !    Roughly Dest(I-SHIFT) = Src(I) with periodic boundary conditions
  !    More exactly, Dest(Index(I)) = Src(I) where
  !                  Index(I) = MOD((I-SHIFT-1),N) + 1, I=[1,N]
  !
  !======================================================================

  ! $Id: pgslib_shift_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

  USE PGSLib_Type_MODULE
  USE PGSLib_Permute_MODULE
  USE PGSLib_Reductions_MODULE
  USE PGSLib_Scan_No_Seg_MODULE
  USE PGSLib_Utility_MODULE

  implicit none
  PRIVATE
  PUBLIC :: PGSLib_Global_CSHIFT, PGSLib_Global_EOSHIFT


  INTERFACE PGSLib_Global_CSHIFT
     MODULE PROCEDURE CSHIFT_Int_1D
     MODULE PROCEDURE CSHIFT_Sin_1D
     MODULE PROCEDURE CSHIFT_Dbl_1D
     MODULE PROCEDURE CSHIFT_Log_1D
  END INTERFACE

  INTERFACE PGSLib_Global_EOSHIFT
     MODULE PROCEDURE EOSHIFT_Int_1D
     MODULE PROCEDURE EOSHIFT_Sin_1D
     MODULE PROCEDURE EOSHIFT_Dbl_1D
     MODULE PROCEDURE EOSHIFT_Log_1D
  END INTERFACE

CONTAINS
#define _ROUTINE_NAME_ EOSHIFT_Int_1D
#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#include "eoshift.fpp"

#define _ROUTINE_NAME_ EOSHIFT_Sin_1D
#define _DATA_TYPE_ real (PGSLib_Single_Type)
#include "eoshift.fpp"

#define _ROUTINE_NAME_ EOSHIFT_Dbl_1D
#define _DATA_TYPE_ real (PGSLib_Double_Type)
#include "eoshift.fpp"

#define _ROUTINE_NAME_ EOSHIFT_Log_1D
#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#include "eoshift.fpp"

#define _ROUTINE_NAME_ CSHIFT_Int_1D
#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#include "cshift.fpp"

#define _ROUTINE_NAME_ CSHIFT_Sin_1D
#define _DATA_TYPE_ real (PGSLib_Single_Type)
#include "cshift.fpp"

#define _ROUTINE_NAME_ CSHIFT_Dbl_1D
#define _DATA_TYPE_ real (PGSLib_Double_Type)
#include "cshift.fpp"

#define _ROUTINE_NAME_ CSHIFT_Log_1D
#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#include "cshift.fpp"

end MODULE PGSLib_Shift_MODULE


