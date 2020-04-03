!!CPP!! This file is included by pgslib_scatter_module.F and pgslib_gather_module.F
!!CPP!! It provides the arithemetic/numerical scatter and gather interfaces.

!!CPP!! $Id: arithmetic_gath-scatt.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

!!CPP!! The contortions in this macro are to get op and specific
!!CPP!! to expand before catenation.
#define IDENTITY(a) a
#define _GEN_NAME_(op,specific) IDENTITY(op)IDENTITY(specific)

!======================================================================
!          Gather/Scatter_1_1_Int(Dest, Src, Index, Mask, Trace)
! PURPOSE:
!          Scatter from 1D source to 1D destination
!
!======================================================================
#define _SPECIFIC_SUFFIX_ _1_1_Int
#define _ROUTINE_NAME_  _GEN_NAME_(_OP_NAME_,_SPECIFIC_SUFFIX_)
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _OP_NULL_       0
#define _OP_ID_         _GEN_OP_ID_( _DATA_TYPE_, _OP_NULL_ )
#define _NARROW_DIM_ dimension(:)
#define _FAT_DIM_ dimension(:)
#define _INDEX_SIZE_      SIZE(Index,1)
#define _ARRAY_ELEMENT_ i
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    1
#include _INCLUDE_FILE_
#undef  _SPECIFIC_SUFFIX_
#undef  _OP_NULL_
#undef  _FAT_DIM_
#undef  _NARROW_DIM_

!======================================================================
!          Gather/Scatter_1_1_Single(Dest, Src, Index, Mask, Trace)
! PURPOSE:
!          Scatter from 1D source to 1D destination
!
!======================================================================
#define _SPECIFIC_SUFFIX_ _1_1_Single
#define _ROUTINE_NAME_  _GEN_NAME_(_OP_NAME_,_SPECIFIC_SUFFIX_)
#define _DATA_TYPE_     real (PGSLib_Single_Type)
#define _OP_NULL_       REAL(0.0, PGSLib_Single_Type)
#define _OP_ID_         _GEN_OP_ID_( _DATA_TYPE_, _OP_NULL_ )
#define _NARROW_DIM_ dimension(:)
#define _FAT_DIM_ dimension(:)
#define _INDEX_SIZE_      SIZE(Index,1)
#define _ARRAY_ELEMENT_ i
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    1
#include _INCLUDE_FILE_
#undef  _SPECIFIC_SUFFIX_
#undef  _OP_NULL_
#undef  _FAT_DIM_
#undef  _NARROW_DIM_

!======================================================================
!          Gather/Scatter_1_1_Double(Dest, Src, Index, Mask, Trace)
! PURPOSE:
!          Scatter from 1D source to 1D destination
!
!======================================================================
#define _SPECIFIC_SUFFIX_ _1_1_Double
#define _ROUTINE_NAME_  _GEN_NAME_(_OP_NAME_,_SPECIFIC_SUFFIX_)
#define _DATA_TYPE_     real (PGSLib_Double_Type)
#define _OP_NULL_       REAL(0.0, PGSLib_Double_Type)
#define _OP_ID_         _GEN_OP_ID_( _DATA_TYPE_, _OP_NULL_ )
#define _NARROW_DIM_ dimension(:)
#define _FAT_DIM_ dimension(:)
#define _INDEX_SIZE_      SIZE(Index,1)
#define _ARRAY_ELEMENT_ i
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    1
#include _INCLUDE_FILE_
#undef  _SPECIFIC_SUFFIX_
#undef  _OP_NULL_
#undef  _FAT_DIM_
#undef  _NARROW_DIM_

!======================================================================
!          Gather/Scatter_1_2_Int(Dest, Src, Index, Mask, Trace)
! PURPOSE:
!          Scatter from 2D source to 1D destination
!
!======================================================================
#define _SPECIFIC_SUFFIX_ _1_2_Int
#define _ROUTINE_NAME_  _GEN_NAME_(_OP_NAME_,_SPECIFIC_SUFFIX_)
#define _DATA_TYPE_     integer (PGSLib_Int_Type)
#define _OP_NULL_       0_PGSLib_Int_Type
#define _OP_ID_         _GEN_OP_ID_( _DATA_TYPE_, _OP_NULL_ )
#define _NARROW_DIM_    dimension(:)
#define _FAT_DIM_       dimension(:,:)
#define _INDEX_SIZE_    SIZE(Index,1), SIZE(Index,2)
#define _ARRAY_ELEMENT_ i,j
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    SIZE(Local_Index,2)
#include _INCLUDE_FILE_
#undef  _SPECIFIC_SUFFIX_
#undef  _OP_NULL_
#undef  _FAT_DIM_
#undef  _NARROW_DIM_

!======================================================================
!          Gather/Scatter_1_2_Single(Dest, Src, Index, Mask, Trace)
! PURPOSE:
!          Scatter from 2D source to 1D destination
!
!======================================================================
#define _SPECIFIC_SUFFIX_ _1_2_Single
#define _ROUTINE_NAME_  _GEN_NAME_(_OP_NAME_,_SPECIFIC_SUFFIX_)
#define _DATA_TYPE_     real (PGSLib_Single_Type)
#define _OP_NULL_       REAL(0.0, PGSLib_Single_Type)
#define _OP_ID_         _GEN_OP_ID_( _DATA_TYPE_, _OP_NULL_ )
#define _NARROW_DIM_ dimension(:)
#define _FAT_DIM_ dimension(:,:)
#define _INDEX_SIZE_      SIZE(Index,1), SIZE(Index,2)
#define _ARRAY_ELEMENT_ i,j
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    SIZE(Local_Index,2)
#include _INCLUDE_FILE_
#undef  _SPECIFIC_SUFFIX_
#undef  _OP_NULL_
#undef  _FAT_DIM_
#undef  _NARROW_DIM_

!======================================================================
!          Gather/Scatter_1_2_Double(Dest, Src, Index, Mask, Trace)
! PURPOSE:
!          Scatter from 2D source to 1D destination
!
!======================================================================
#define _SPECIFIC_SUFFIX_ _1_2_Double
#define _ROUTINE_NAME_  _GEN_NAME_(_OP_NAME_,_SPECIFIC_SUFFIX_)
#define _DATA_TYPE_     real (PGSLib_Double_Type)
#define _OP_NULL_       REAL(0.0, PGSLib_Double_Type)
#define _OP_ID_         _GEN_OP_ID_( _DATA_TYPE_, _OP_NULL_ )
#define _NARROW_DIM_ dimension(:)
#define _FAT_DIM_ dimension(:,:)
#define _INDEX_SIZE_      SIZE(Index,1), SIZE(Index,2)
#define _ARRAY_ELEMENT_ i,j
#define _I_LOOP_MAX_    SIZE(Local_Index,1)
#define _J_LOOP_MAX_    SIZE(Local_Index,2)
#include _INCLUDE_FILE_
#undef  _SPECIFIC_SUFFIX_
#undef  _OP_NULL_
#undef  _FAT_DIM_
#undef  _NARROW_DIM_

#undef  _GEN_NAME_
#undef  _OP_
#undef  _OP_NAME_
#undef  _GEN_OP_ID_
#undef _INCLUDE_FILE_
