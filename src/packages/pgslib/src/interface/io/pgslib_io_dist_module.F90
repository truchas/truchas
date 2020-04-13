MODULE PGSLib_IO_Dist_MODULE
  use pgslib_c_binding
  use PGSLib_Type_MODULE
  use pgslib_io_collate_module
  USE pgslib_globals_module
  use PGSLib_Utility_MODULE,     only : pgslib_check_error,      &
      PGSLib_Inquire_IO_P,     &
      PGSLib_Inquire_nPE
  use,intrinsic :: iso_fortran_env, only: int8
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_Dist

  !  The routines supported in this module are
  !  PGSLib_Dist
  !
  !
  ! For each of these generic routines, scalar and vector versions,
  !          as well as four or five data types are supported.
  ! The data types are
  !
  !          PGSLib_INT_TYPE    (typically the default integer on each system)
  !          PGSLib_REAL_TYPE   (typically the default single precision)
  !          PGSLib_DOUBLE_TYPE (typically the default double precision)
  !          PGSLib_Log_Type    (typeically the default logical type)
  !          CHARACTER          (for some of the routines)
  ! The types are set in the module PGSLib_Types_MODULE, except for CHARACTER, which is the
  !          default character type.
  !

  ! $Id: pgslib_io_dist_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PGSLib_Dist Interfaces  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE PGSLib_Dist

     ! Broadcast Scalars
     ! These routines provide F90 interfaces to the C routines
     MODULE PROCEDURE PGSLib_Dist_Int8_Scalar_F
     MODULE PROCEDURE PGSLib_Dist_Int_Scalar_F
     MODULE PROCEDURE PGSLib_Dist_REAL_Scalar_F
     MODULE PROCEDURE PGSLib_Dist_Double_Scalar_F
     MODULE PROCEDURE PGSLib_Dist_Log_Scalar_F

     ! Broadcast Vectors
     ! Interfaces for the F90 interface routines, which take a single vector argument
     MODULE Procedure PGSLib_Dist_Int8_Vector_F
     MODULE Procedure PGSLib_Dist_Int_Vector_F
     MODULE Procedure PGSLib_Dist_Real_Vector_F
     MODULE Procedure PGSLib_Dist_Double_Vector_F
     MODULE Procedure PGSLib_Dist_Log_Vector_F

  END INTERFACE!For the user callable PGSLib_Dist

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PGSLib_Dist SUBROUTINE Bodies !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Distribute scalars, one per PE.
  ! These are the F90 versions, which do some checking, and then
  !    call C working routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_dist_int8_scalar_f(scalar_out, scalarv_in)

#define _DATA_TYPE_ integer (int8)
#include "dist_scalar.fpp"

  END SUBROUTINE PGSLib_dist_int8_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_dist_int_scalar_f(scalar_out, scalarv_in)

#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#include "dist_scalar.fpp"

  END SUBROUTINE PGSLib_dist_int_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE pgslib_dist_real_scalar_f(scalar_out, scalarv_in)

#define _DATA_TYPE_ real (PGSLib_Real_Type)
#include "dist_scalar.fpp"

    return
  end SUBROUTINE pgslib_dist_real_scalar_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE pgslib_dist_double_scalar_f(scalar_out, scalarv_in)

#define _DATA_TYPE_ real (PGSLib_Double_Type)
#include "dist_scalar.fpp"

  end SUBROUTINE pgslib_dist_double_scalar_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE pgslib_dist_log_scalar_f(scalar_out, scalarv_in)

#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#include "dist_scalar.fpp"

  end SUBROUTINE pgslib_dist_log_scalar_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Distribute vectors.
  ! These routines do some checking, and then call C routines
  !          to do the work.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Dist_Int8_Vector_F(Vector_out, Vector_In, Lengths)

#define _DATA_TYPE_ integer (int8)
#include "dist_vector.fpp"

  END SUBROUTINE PGSLib_Dist_Int8_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Dist_Int_Vector_F(Vector_out, Vector_In, Lengths)

#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#include "dist_vector.fpp"

  END SUBROUTINE PGSLib_Dist_Int_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Dist_Real_Vector_F(Vector_out, Vector_In, Lengths)

#define _DATA_TYPE_ real (PGSLib_Real_Type)
#include "dist_vector.fpp"

  END SUBROUTINE PGSLib_Dist_Real_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Dist_Double_Vector_F(Vector_out, Vector_In, Lengths)

#define _DATA_TYPE_ real (PGSLib_Double_Type)
#include "dist_vector.fpp"

  END SUBROUTINE PGSLib_Dist_Double_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Dist_Log_Vector_F(Vector_out, Vector_In, Lengths)

#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#include "dist_vector.fpp"

  END SUBROUTINE PGSLib_Dist_Log_Vector_F

END MODULE PGSLib_IO_Dist_MODULE
