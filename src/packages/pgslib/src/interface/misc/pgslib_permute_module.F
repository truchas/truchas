MODULE PGSLib_Permute_MODULE
  !======================================================================
  ! PURPOSE -
  !   Routines to permute an array according to an index.
  !   Given Src(N), Dest(N), Index(N), global arrays (N = global size)
  !   Return Dest(Index) = Src, in the global sense.
  !======================================================================

  ! $Id: pgslib_permute_module.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

!  USE PGSLib_Type_MODULE

  implicit none
  
  PRIVATE
  PUBLIC :: PGSLib_Permute, PGSLib_Redistribute, PGSLib_Global_Pack

  INTERFACE PGSLib_Permute
     MODULE PROCEDURE PGSLib_Permute_Int_1D
     MODULE PROCEDURE PGSLib_Permute_Real_1D
     MODULE PROCEDURE PGSLib_Permute_Double_1D
     MODULE PROCEDURE PGSLib_Permute_Log_1D
  END INTERFACE

  INTERFACE PGSLib_Redistribute
     MODULE PROCEDURE PGSLib_Redist_Int_1D
     MODULE PROCEDURE PGSLib_Redist_Real_1D
     MODULE PROCEDURE PGSLib_Redist_Double_1D
     MODULE PROCEDURE PGSLib_Redist_Log_1D
  END INTERFACE

  INTERFACE PGSLib_Global_Pack
     MODULE PROCEDURE PGSLib_Pack_Int_1D
     MODULE PROCEDURE PGSLib_Pack_Real_1D
     MODULE PROCEDURE PGSLib_Pack_Double_1D
     MODULE PROCEDURE PGSLib_Pack_Log_1D
  END INTERFACE

CONTAINS
#define _ROUTINE_NAME_ PGSLib_Permute_Int_1D
#define _DATA_TYPE_    integer (PGSLib_Int_Type)
#define _OP_ID_        0
#define _SCATTER_OP_   pgslib_scatter_sum
#define _COMP_OP_      == 
#include "permute.fpp"


#define _ROUTINE_NAME_ PGSLib_Permute_Real_1D
#define _DATA_TYPE_    REAL (PGSLib_Real_Type)
#define _OP_ID_        0.0_PGSlib_Real_Type
#define _SCATTER_OP_   pgslib_scatter_sum
#define _COMP_OP_      == 
#include "permute.fpp"


#define _ROUTINE_NAME_ PGSLib_Permute_Double_1D
#define _DATA_TYPE_    REAL (PGSLib_Double_Type)
#define _OP_ID_        0.0_PGSlib_Double_Type
#define _SCATTER_OP_   pgslib_scatter_sum
#define _COMP_OP_      == 
#include "permute.fpp"


#define _ROUTINE_NAME_ PGSLib_Permute_Log_1D
#define _DATA_TYPE_    logical (PGSLib_Log_Type)
#define _OP_ID_        .FALSE.
#define _SCATTER_OP_    pgslib_scatter_or
#define _COMP_OP_      .EQV.
#include "permute.fpp"

#define _ROUTINE_NAME_ PGSLib_Redist_Int_1D
#define _DATA_TYPE_    integer (PGSLib_Int_Type)
#define _OP_ID_        0
#define _SCATTER_OP_   pgslib_scatter_sum
#include "redist.fpp"


#define _ROUTINE_NAME_ PGSLib_Redist_Real_1D
#define _DATA_TYPE_    REAL (PGSLib_Real_Type)
#define _OP_ID_        0.0_PGSlib_Real_Type
#define _SCATTER_OP_   pgslib_scatter_sum
#include "redist.fpp"


#define _ROUTINE_NAME_ PGSLib_Redist_Double_1D
#define _DATA_TYPE_    REAL (PGSLib_Double_Type)
#define _OP_ID_        0.0_PGSlib_Double_Type
#define _SCATTER_OP_   pgslib_scatter_sum
#include "redist.fpp"


#define _ROUTINE_NAME_ PGSLib_Redist_Log_1D
#define _DATA_TYPE_    logical (PGSLib_Log_Type)
#define _OP_ID_        .FALSE.
#define _SCATTER_OP_    pgslib_scatter_or
#include "redist.fpp"



  subroutine PGSLib_Pack_Int_1D(Dest, Source, MASK)
#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#define _OP_ID_     0
#define _SCATTER_OP_ pgslib_scatter_sum
#define _RESULT_ Dest
#include "pack.fpp"
  end subroutine PGSLib_Pack_Int_1D

  subroutine PGSLib_Pack_Real_1D(Dest, Source, MASK)
#define _DATA_TYPE_ real (PGSLib_Real_Type)
#define _OP_ID_     0.0_PGslib_Real_Type
#define _SCATTER_OP_ pgslib_scatter_sum
#define _RESULT_ Dest
#include "pack.fpp"
  end subroutine PGSLib_Pack_Real_1D

  subroutine PGSLib_Pack_Double_1D(Dest, Source, MASK)
#define _DATA_TYPE_ real (PGSLib_Double_Type)
#define _OP_ID_     0.0_PGSLib_Double_Type
#define _SCATTER_OP_ pgslib_scatter_sum
#define _RESULT_ Dest
#include "pack.fpp"
  end subroutine PGSLib_Pack_Double_1D

  subroutine PGSLib_Pack_Log_1D(Dest, Source, MASK)
#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#define _OP_ID_     .FALSE.
#define _SCATTER_OP_ pgslib_scatter_or
#define _RESULT_ Dest
#include "pack.fpp"
  end subroutine PGSLib_Pack_Log_1D



end MODULE PGSLib_Permute_MODULE
