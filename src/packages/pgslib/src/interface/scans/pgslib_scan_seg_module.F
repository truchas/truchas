MODULE PGSLib_SCAN_SEG_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Provide PREFIX and SUFFIX operations, segment versions.
  ! These use CSHIFT, CSHIFT depends on the no-segment versions,
  ! and that is the reason these are split from the no-segment versions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_scan_seg_module.F,v 1.1.1.1 2000/10/11 22:44:30 ferrell Exp $

  use pgslib_type_module
  use pgslib_globals_module
  use pgslib_utility_module, only : PGSLib_Fatal_ERROR, PGSLib_Scope_Check

  implicit none
  PRIVATE
  PUBLIC :: PGSLib_SUM_PREFIX, PGSLib_SUM_SUFFIX

  INTERFACE PGSLib_SUM_PREFIX
     MODULE PROCEDURE SUM_PREFIX_SEG_INT_1D
     MODULE PROCEDURE SUM_PREFIX_SEG_SINGLE_1D
     MODULE PROCEDURE SUM_PREFIX_SEG_DOUBLE_1D
  END INTERFACE
     
  INTERFACE PGSLib_SUM_SUFFIX
     MODULE PROCEDURE SUM_SUFFIX_SEG_INT_1D
     MODULE PROCEDURE SUM_SUFFIX_SEG_SINGLE_1D
     MODULE PROCEDURE SUM_SUFFIX_SEG_DOUBLE_1D
  END INTERFACE
     
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUM_PREFIX Top_level routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _ROUTINE_NAME_          SUM_PREFIX_SEG_INT_1D
#define _SCAN_DATA_TYPE_        integer (PGSLib_Int_Type)
#define _SHIFT_DIRECTION_       -1
#define _GLOBAL_MAX_OR_MIN_VAL_ PGSLib_Global_MINVAL
#define _SEG_BIT_SCAN_          PGSLib_SEG_SUM_PREFIX

#include "scan-segment-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_          SUM_PREFIX_SEG_SINGLE_1D
#define _SCAN_DATA_TYPE_        real (PGSLib_Real_Type)
#define _SHIFT_DIRECTION_       -1
#define _GLOBAL_MAX_OR_MIN_VAL_ PGSLib_Global_MINVAL
#define _SEG_BIT_SCAN_          PGSLib_SEG_SUM_PREFIX

#include "scan-segment-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_          SUM_PREFIX_SEG_DOUBLE_1D
#define _SCAN_DATA_TYPE_        real (PGSLib_Double_Type)
#define _SHIFT_DIRECTION_       -1
#define _GLOBAL_MAX_OR_MIN_VAL_ PGSLib_Global_MINVAL
#define _SEG_BIT_SCAN_          PGSLib_SEG_SUM_PREFIX

#include "scan-segment-1d.fpp"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUM_SUFFIX Top_level routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_          SUM_SUFFIX_SEG_INT_1D
#define _SCAN_DATA_TYPE_        integer (PGSLib_Int_Type)
#define _SHIFT_DIRECTION_       1
#define _GLOBAL_MAX_OR_MIN_VAL_ PGSLib_Global_MAXVAL
#define _SEG_BIT_SCAN_          PGSLib_SEG_SUM_SUFFIX

#include "scan-segment-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_          SUM_SUFFIX_SEG_SINGLE_1D
#define _SCAN_DATA_TYPE_        real (PGSLib_Real_Type)
#define _SHIFT_DIRECTION_       1
#define _GLOBAL_MAX_OR_MIN_VAL_ PGSLib_Global_MAXVAL
#define _SEG_BIT_SCAN_          PGSLib_SEG_SUM_SUFFIX

#include "scan-segment-1d.fpp"

  !======================================================================

#define _ROUTINE_NAME_          SUM_SUFFIX_SEG_DOUBLE_1D
#define _SCAN_DATA_TYPE_        real (PGSLib_Double_Type)
#define _SHIFT_DIRECTION_       1
#define _GLOBAL_MAX_OR_MIN_VAL_ PGSLib_Global_MAXVAL
#define _SEG_BIT_SCAN_          PGSLib_SEG_SUM_SUFFIX

#include "scan-segment-1d.fpp"
  

END MODULE PGSLib_SCAN_SEG_MODULE
