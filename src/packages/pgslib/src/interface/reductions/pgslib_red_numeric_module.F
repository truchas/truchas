MODULE PGSLib_Red_Numeric_MODULE
  use PGSLib_Error_MODULE
  use pgslib_c_binding
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC:: PGSLib_Global_MINVAL
  PUBLIC:: PGSLib_Global_MAXVAL
  PUBLIC:: PGSLib_Global_SUM
  PUBLIC:: PGSLib_Global_COUNT
  PUBLIC:: PGSLib_Global_ALL
  PUBLIC:: PGSLib_Global_ANY
  PUBLIC:: PGSLib_Global_DOT_PRODUCT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines provided in this module
  !       PGSLib_Global_MINVAL
  !       PGSLib_Global_MAXVAL
  !       PGSLib_Global_SUM
  !       PGSLib_Global_ANY
  !       PGSLib_Global_ALL
  !       PGSLib_Global_DOT_PRODUCT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_red_numeric_module.F,v 1.1.1.1 2000/10/11 22:44:29 ferrell Exp $

  INTERFACE PGSLib_Global_MINVAL
     MODULE PROCEDURE PGS_Glbl_MINVAL_INT_Scalar_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_REAL_Scalar_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_DOUBLE_Scalar_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_INT_1D_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_REAL_1D_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_DOUBLE_1D_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_INT_2D_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_REAL_2D_F
     MODULE PROCEDURE PGS_Glbl_MINVAL_DOUBLE_2D_F
  END INTERFACE
  INTERFACE PGSLib_Global_MAXVAL
     MODULE PROCEDURE PGS_Glbl_MAXVAL_INT_Scalar_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_REAL_Scalar_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_DOUBLE_Scalar_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_INT_1D_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_REAL_1D_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_DOUBLE_1D_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_INT_2D_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_REAL_2D_F
     MODULE PROCEDURE PGS_Glbl_MAXVAL_DOUBLE_2D_F
  END INTERFACE
  INTERFACE PGSLib_Global_SUM
     MODULE PROCEDURE PGS_Glbl_SUM_INT_Scalar_F
     MODULE PROCEDURE PGS_Glbl_SUM_REAL_Scalar_F
     MODULE PROCEDURE PGS_Glbl_SUM_DOUBLE_Scalar_F
     MODULE PROCEDURE PGS_Glbl_SUM_INT_1D_F
     MODULE PROCEDURE PGS_Glbl_SUM_REAL_1D_F
     MODULE PROCEDURE PGS_Glbl_SUM_DOUBLE_1D_F
     MODULE PROCEDURE PGS_Glbl_SUM_INT_2D_F
     MODULE PROCEDURE PGS_Glbl_SUM_REAL_2D_F
     MODULE PROCEDURE PGS_Glbl_SUM_DOUBLE_2D_F
  END INTERFACE
  INTERFACE PGSLib_Global_COUNT
     MODULE PROCEDURE PGS_Glbl_COUNT_LOG_Scalar_F
     MODULE PROCEDURE PGS_Glbl_COUNT_Log_1D_F
     MODULE PROCEDURE PGS_Glbl_COUNT_Log_2D_F
  END INTERFACE
  INTERFACE PGSLib_Global_ALL
     MODULE PROCEDURE PGS_Glbl_ALL_LOG_Scalar_F
     MODULE PROCEDURE PGS_Glbl_ALL_Log_1D_F
     MODULE PROCEDURE PGS_Glbl_ALL_Log_2D_F
  END INTERFACE
  INTERFACE PGSLib_Global_ANY
     MODULE PROCEDURE PGS_Glbl_ANY_LOG_Scalar_F
     MODULE PROCEDURE PGS_Glbl_ANY_Log_1D_F
     MODULE PROCEDURE PGS_Glbl_ANY_Log_2D_F
  END INTERFACE
  INTERFACE PGSLib_Global_DOT_PRODUCT
     MODULE PROCEDURE PGSLib_Global_DT_PD_INT_1D
     MODULE PROCEDURE PGSLib_Global_DT_PD_REAL_1D
     MODULE PROCEDURE PGSLib_Global_DT_PD_DOUBLE_1D
     MODULE PROCEDURE PGSLib_Global_DT_PD_LOG_1D
  END INTERFACE

CONTAINS

!!CPP!! Need this strange set of macros to get proper expansion.
#define _STRING_(S)   "S"

!!!!!!!!!! MINVAL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _GENERIC_ROUTINE_NAME_ PGSLib_Global_MINVAL

!!!!!!!!!! Scalar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_INT_Scalar_F
#define _DATA_TYPE_ integer (PGSLib_INT_TYPE) 
#define _ARRAY_  
#define _ARG_       (/ A /)
#include "red_min.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_Real_Scalar_F
#define _DATA_TYPE_ real (PGSLib_Real_TYPE) 
#define _ARRAY_  
#define _ARG_       (/ A /)
#include "red_min.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_Double_Scalar_F
#define _DATA_TYPE_ real (PGSLib_Double_TYPE) 
#define _ARRAY_  
#define _ARG_       (/ A /)
#include "red_min.fpp"

!!!!!!!!!! 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_INT_1D_F
#define _DATA_TYPE_ integer (PGSLib_INT_TYPE) 
#define _ARRAY_     ,Dimension(:)
#define _ARG_       A
#include "red_min.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_Real_1D_F
#define _DATA_TYPE_ real (PGSLib_Real_TYPE) 
#define _ARRAY_     ,Dimension(:)
#define _ARG_       A
#include "red_min.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_Double_1D_F
#define _DATA_TYPE_ real (PGSLib_Double_TYPE) 
#define _ARRAY_     ,Dimension(:)
#define _ARG_       A
#include "red_min.fpp"

!!!!!!!!!! 2D_F MINVAL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_INT_2D_F
#define _DATA_TYPE_ integer (PGSLib_INT_TYPE) 
#define _ARRAY_     ,Dimension(:,:)
#define _ARG_       A
#include "red_min.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_Real_2D_F
#define _DATA_TYPE_ real (PGSLib_Real_TYPE) 
#define _ARRAY_     ,Dimension(:,:)
#define _ARG_       A
#include "red_min.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MINVAL_Double_2D_F
#define _DATA_TYPE_ real (PGSLib_Double_TYPE) 
#define _ARRAY_     ,Dimension(:,:)
#define _ARG_       A
#include "red_min.fpp"

#undef _GENERIC_ROUTINE_NAME_
!!!!!!!!!! MAXVAL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _GENERIC_ROUTINE_NAME_ PGSLib_Global_MAXVAL

!!!!!!!!!! Scalar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_INT_Scalar_F
#define _DATA_TYPE_ integer (PGSLib_INT_TYPE) 
#define _ARRAY_  
#define _ARG_       (/ A /)
#include "red_max.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_Real_Scalar_F
#define _DATA_TYPE_ real (PGSLib_Real_TYPE) 
#define _ARRAY_  
#define _ARG_       (/ A /)
#include "red_max.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_Double_Scalar_F
#define _DATA_TYPE_ real (PGSLib_Double_TYPE) 
#define _ARRAY_  
#define _ARG_       (/ A /)
#include "red_max.fpp"

!!!!!!!!!! 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_INT_1D_F
#define _DATA_TYPE_ integer (PGSLib_INT_TYPE) 
#define _ARRAY_     ,Dimension(:)
#define _ARG_       A
#include "red_max.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_Real_1D_F
#define _DATA_TYPE_ real (PGSLib_Real_TYPE) 
#define _ARRAY_     ,Dimension(:)
#define _ARG_       A
#include "red_max.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_Double_1D_F
#define _DATA_TYPE_ real (PGSLib_Double_TYPE) 
#define _ARRAY_     ,Dimension(:)
#define _ARG_       A
#include "red_max.fpp"

!!!!!!!!!! 2D_F MAXVAL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_INT_2D_F
#define _DATA_TYPE_ integer (PGSLib_INT_TYPE) 
#define _ARRAY_     ,Dimension(:,:)
#define _ARG_       A
#include "red_max.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_Real_2D_F
#define _DATA_TYPE_ real (PGSLib_Real_TYPE) 
#define _ARRAY_     ,Dimension(:,:)
#define _ARG_       A
#include "red_max.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_   PGS_Glbl_MAXVAL_Double_2D_F
#define _DATA_TYPE_ real (PGSLib_Double_TYPE) 
#define _ARRAY_     ,Dimension(:,:)
#define _ARG_       A
#include "red_max.fpp"

#undef _GENERIC_ROUTINE_NAME_
!!!!!!!!!! SUM ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _GENERIC_ROUTINE_NAME_ PGSLib_Global_SUM

!!!!!!!!!! Scalar Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#define _ARRAY_     
#define _ARG_       (/ A /)
#define _ROUTINE_   PGS_Glbl_SUM_INT_Scalar_F
#include "red_sum.fpp"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ real (PGSLib_REAL_TYPE)
#define _ARRAY_     
#define _ARG_       (/ A /)
#define _ROUTINE_   PGS_Glbl_SUM_Real_Scalar_F
#include "red_sum.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ real (PGSLib_Double_TYPE)
#define _ARRAY_     
#define _ARG_       (/ A /)
#define _ROUTINE_   PGS_Glbl_SUM_Double_Scalar_F
#include "red_sum.fpp"

!!!!!!!!!! 1D Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#define _ARRAY_     , Dimension(:)
#define _ARG_       A 
#define _ROUTINE_   PGS_Glbl_SUM_INT_1D_F
#include "red_sum.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ real (PGSLib_REAL_TYPE)
#define _ARRAY_     , Dimension(:) 
#define _ARG_       A
#define _ROUTINE_   PGS_Glbl_SUM_Real_1D_F
#include "red_sum.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ real (PGSLib_Double_TYPE)
#define _ARRAY_     , Dimension(:)
#define _ARG_       A
#define _ROUTINE_   PGS_Glbl_SUM_Double_1D_F
#include "red_sum.fpp"

!!!!!!!!!! SUM ROUTINES, 2D_F verions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#define _ARRAY_     , Dimension(:,:)
#define _ARG_       A 
#define _ROUTINE_   PGS_Glbl_SUM_INT_2D_F
#include "red_sum.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ real (PGSLib_REAL_TYPE)
#define _ARRAY_     , Dimension(:,:) 
#define _ARG_       A
#define _ROUTINE_   PGS_Glbl_SUM_Real_2D_F
#include "red_sum.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DATA_TYPE_ real (PGSLib_Double_TYPE)
#define _ARRAY_     , Dimension(:,:)
#define _ARG_       A
#define _ROUTINE_   PGS_Glbl_SUM_Double_2D_F
#include "red_sum.fpp"

#undef _GENERIC_ROUTINE_NAME_
!!!!!!!!!! COUNT ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function PGS_Glbl_COUNT_LOG_Scalar_F(MASK)
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    integer (PGSLib_Int_Type) PGS_Glbl_COUNT_LOG_Scalar_F
    logical (PGSLib_Log_Type), INTENT(in) :: MASK

    ! Local temporaries
    ! Notice that this must be integer type to pass through to C properly
    integer (PGSLib_Int_TYPE) LocalCount

    ! On-PE reduction
    LocalCount = MERGE(1,0,Mask)

    ! Global reduction, with broadcast of result
    PGS_Glbl_COUNT_LOG_Scalar_F = PGSLib_Global_SUM(LocalCount)

    RETURN
  END FUNCTION PGS_Glbl_COUNT_LOG_Scalar_F

!!!!!!!!!! 1D COUNT ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function PGS_Glbl_COUNT_LOG_1D_F(MASK)
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    integer (PGSLib_Int_Type) PGS_Glbl_COUNT_LOG_1D_F
    logical (PGSLib_Log_Type), INTENT(in), DIMENSION(:):: MASK

    ! Local temporaries
    ! Notice that this must be integer type to pass through to C properly
    integer (PGSLib_Int_TYPE) LocalCount

    ! On-PE reduction
    LocalCount = COUNT(Mask)

    ! Global reduction, with broadcast of result
    PGS_Glbl_COUNT_LOG_1D_F = PGSLib_Global_SUM(LocalCount)

    RETURN
  END FUNCTION PGS_Glbl_COUNT_Log_1D_F


!!!!!!!!!! 2D COUNT ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function PGS_Glbl_COUNT_LOG_2D_F(MASK)
    USE PGSLib_Type_MODULE
    IMPLICIT NONE
    integer (PGSLib_Int_Type) PGS_Glbl_COUNT_LOG_2D_F
    logical (PGSLib_Log_Type), INTENT(in), DIMENSION(:,:):: MASK

    ! Local temporaries
    ! Notice that this must be integer type to pass through to C properly
    integer (PGSLib_Int_TYPE) LocalCount

    ! On-PE reduction
    LocalCount = COUNT(Mask)

    ! Global reduction, with broadcast of result
    PGS_Glbl_COUNT_Log_2D_F = PGSLib_Global_SUM(LocalCount)

    RETURN
  END FUNCTION PGS_Glbl_COUNT_Log_2D_F


!!!!!!!!!! ALL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _GENERIC_ROUTINE_NAME_ PGSLib_Global_ALL


#define _ROUTINE_     PGS_Glbl_ALL_LOG_Scalar_F
#define _MASK_SHAPE_ 
#define _MASK_        (/ MASK /)
#include "red_all.fpp"

!!!!!!!!!! 1D ALL ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_     PGS_Glbl_ALL_LOG_1D_F
#define _MASK_SHAPE_  , DIMENSION(:)
#define _MASK_        MASK
#include "red_all.fpp"

!!!!!!!!!! 2D ALL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_     PGS_Glbl_ALL_LOG_2D_F
#define _MASK_SHAPE_  , DIMENSION(:,:)
#define _MASK_        MASK
#include "red_all.fpp"

#undef _GENERIC_ROUTINE_NAME_
!!!!!!!!!! ANY ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _GENERIC_ROUTINE_NAME_ PGSLib_Global_ANY

#define _ROUTINE_     PGS_Glbl_ANY_LOG_Scalar_F
#define _MASK_SHAPE_ 
#define _MASK_        (/ MASK /)
#include "red_any.fpp"

!!!!!!!!!! 1D_F Log Routine

#define _ROUTINE_     PGS_Glbl_ANY_LOG_1D_F
#define _MASK_SHAPE_  , DIMENSION(:)
#define _MASK_        MASK
#include "red_any.fpp"

!!!!!!!!!! 2D_F ANY ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_     PGS_Glbl_ANY_LOG_2D_F
#define _MASK_SHAPE_  , DIMENSION(:,:)
#define _MASK_        MASK
#include "red_any.fpp"


#undef _GENERIC_ROUTINE_NAME_
!!!!!!!!!! DOT_PRODUCT ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _GENERIC_ROUTINE_NAME_ PGSLib_Global_DOT_PRODUCT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_    PGSLib_Global_DT_PD_INT_1D
#define _DATA_TYPE_  integer (PGSLib_Int_Type)
#include "red_dtp.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_    PGSLib_Global_DT_PD_Real_1D
#define _DATA_TYPE_  real (PGSLib_Real_Type)
#include "red_dtp.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _ROUTINE_    PGSLib_Global_DT_PD_Double_1D
#define _DATA_TYPE_  real (PGSLib_Double_Type)
#include "red_dtp.fpp"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Global_DT_PD_LOG_1D(Vector_A, Vector_B, SCOPE)
    USE PGSLib_Type_MODULE
    use pgslib_globals_module
    IMPLICIT NONE
    LOGICAL (PGSLib_LOG_TYPE) PGSLib_Global_DT_PD_LOG_1D
    LOGICAL (PGSLib_LOG_TYPE), INTENT(in), DIMENSION(:):: Vector_A
    LOGICAL (PGSLib_LOG_TYPE), INTENT(in), DIMENSION(:):: Vector_B
    type (PGSLib_SCOPE),                       &
         intent(in   ),                &
         optional       :: SCOPE


    ! Local temporaries
    logical (PGSLib_LOG_TYPE) :: AdotB
    logical (PGSLib_Log_TYPE) :: Global

    ! Local or Global operation
#include "red_global_test.fpp"

    ! On-PE reduction
    AdotB = DOT_PRODUCT(Vector_A, Vector_B)

    If (Global) then
       PGSLib_Global_DT_PD_LOG_1D = PGSLib_Global_ALL(AdotB)
    ELSE
       PGSLib_Global_DT_PD_LOG_1D = AdotB
    END If

    RETURN
  END FUNCTION PGSLib_Global_DT_PD_LOG_1D

#undef _GENERIC_ROUTINE_NAME_

END MODULE PGSLib_Red_Numeric_MODULE
    
