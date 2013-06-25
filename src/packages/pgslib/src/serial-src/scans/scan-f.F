!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Scan Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! This file contains the DUMMY, serial versions.  This file
!     is replaced by a .c file containing MPI calls in the
!     actual library.
!
!


! $Id: scan-f.F,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PREFIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define _ROUTINE_NAME_   Off_Node_SUM_PREFIX_INT_C
#define PGSLIB_DATA_TYPE integer

#include "scan-f.fpp"

#define _ROUTINE_NAME_   Off_Node_SUM_PREFIX_FLOAT_C
#define PGSLIB_DATA_TYPE real

#include "scan-f.fpp"

#define _ROUTINE_NAME_   Off_Node_SUM_PREFIX_DOUBLE_C
#define PGSLIB_DATA_TYPE real(KIND(1.0d0))

#include "scan-f.fpp"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! SUFFIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _ROUTINE_NAME_   Off_Node_SUM_SUFFIX_INT_C
#define PGSLIB_DATA_TYPE integer

#include "scan-f.fpp"

#define _ROUTINE_NAME_   Off_Node_SUM_SUFFIX_FLOAT_C
#define PGSLIB_DATA_TYPE real

#include "scan-f.fpp"

#define _ROUTINE_NAME_   Off_Node_SUM_SUFFIX_DOUBLE_C
#define PGSLIB_DATA_TYPE real(KIND(1.0d0))

#include "scan-f.fpp"
