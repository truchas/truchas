MODULE PGSLib_IO_Collate_MODULE
  use pgslib_c_binding
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_Collate

  !  The routines supported in this module are
  !  PGSLib_Collate
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

  ! $Id: pgslib_io_collate_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PGSLib_Collate Interfaces  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE PGSLib_Collate

     ! Broadcast Scalars
     ! These routines provide F90 interfaces to the C routines

     MODULE PROCEDURE PGSLib_Collate_Int_Scalar_F
     MODULE PROCEDURE PGSLib_Collate_REAL_Scalar_F
     MODULE PROCEDURE PGSLib_Collate_Double_Scalar_F
     MODULE PROCEDURE PGSLib_Collate_Log_Scalar_F
     MODULE PROCEDURE PGSLib_Collate_Char_Scalar_F

     ! Broadcast Vectors
     ! Interfaces for the F90 interface routines, which take a single vector argument
     MODULE Procedure PGSLib_Collate_Int_Vector_F
     MODULE Procedure PGSLib_Collate_Real_Vector_F
     MODULE Procedure PGSLib_Collate_Double_Vector_F
     MODULE Procedure PGSLib_Collate_Log_Vector_F
     MODULE Procedure PGSLib_Collate_Char_Vector_F

  END INTERFACE!For the user callable PGSLib_Collate

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PGSLib_Collate SUBROUTINE Bodies !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Collate scalars, one per PE.
  ! These are the F90 versions, which do some checking, and then
  !    call C working routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Collate_int_scalar_f(scalarv_out, scalar_in)

#define _DATA_TYPE_ integer (PGSLib_Int_Type)
#include "collate_scalar.fpp"

  END SUBROUTINE PGSLib_Collate_int_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE pgslib_Collate_real_scalar_f(scalarv_out, scalar_in)

#define _DATA_TYPE_ real (PGSLib_Real_Type)
#include "collate_scalar.fpp"

  end SUBROUTINE pgslib_Collate_real_scalar_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE pgslib_Collate_double_scalar_f(scalarv_out, scalar_in)

#define _DATA_TYPE_ real (PGSLib_Double_Type)
#include "collate_scalar.fpp"

  end SUBROUTINE pgslib_Collate_double_scalar_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE pgslib_Collate_log_scalar_f(scalarv_out, scalar_in)

#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#include "collate_scalar.fpp"

  end SUBROUTINE pgslib_Collate_log_scalar_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE pgslib_Collate_char_scalar_f(scalarv_out, scalar_in)
    USE PGSLib_Type_MODULE
    USE PGSLib_Utility_MODULE,  ONLY : pgslib_check_error,     &
         &                             PGSLib_Inquire_thisPE,  &
         &                             PGSLib_Inquire_nPE,     &
         &                             PGSLib_Inquire_IO_P

    USE PGSLib_Red_Numeric_Module, ONLY:PGSLib_Global_MAXVAL,   &
         &                             PGSLib_Global_MINVAL
    implicit none
    CHARACTER (LEN = *), intent(  OUT), dimension(:):: scalarv_out
    CHARACTER (LEN = *), intent(IN   )              :: scalar_in
    ! Local variables 
    integer in_len, max_in_len, min_in_len, out_len
    integer i, start, c
    integer vec_len
    integer (PGSLib_Int_Type), POINTER, dimension(:) :: lengths
    logical error_flag
    CHARACTER (LEN=1), dimension(LEN(scalar_in)) :: vector_in
    CHARACTER (LEN=1), POINTER, dimension(:) :: vector_out_total

    ! Check that all strings, in and out, have same length.
    ! Check in lengths
    in_len  = LEN(scalar_in)
    max_in_len = PGSLib_Global_MAXVAL(in_len)
    min_in_len = PGSLib_Global_MINVAL(in_len)

    error_flag = (in_len /= max_in_len) .OR. (in_len /= min_in_len)
    call pgslib_check_error(error_flag, "String lengths of all sources must be the same in PGSLib_Collate.")

    ! Check out length, but only on IO PE
    error_flag = .FALSE.
    if (PGSLib_Inquire_IO_P()) then
       out_len = LEN(scalarv_out)
       error_flag = (out_len /= in_len)
    end if
    call pgslib_check_error(error_flag, "String lengths of source and destination must be the same in PGSLib_Collate.")

    ! Check that destination is large enough
    error_flag = .FALSE.
    if (PGSLib_Inquire_IO_P()) then
       vec_len = SIZE(scalarv_out)
       error_flag = (vec_len .LT. PGSLib_Inquire_nPE())
    endif
    call pgslib_check_error(error_flag, "Desintation array too small in pgslib_Collate for log scalar")

    ! Setup lengths for receiving
    ALLOCATE(lengths(PGSLib_Inquire_nPE()))
    call pgslib_collate(lengths, in_len)

    ! Now do "vector" collate, each string is a vector of characters
    do i = 1, SIZE(vector_in,1)
       vector_in(i) = scalar_in(i:i)
    enddo
    if (PGSLib_Inquire_IO_P()) then
       allocate(vector_out_total(vec_len*out_len))
    else
       allocate(vector_out_total(1))
    endif

    call pgslib_collate_vector_c(vector_out_total, lengths, vector_in, in_len)

    if (PGSLib_Inquire_IO_P()) then
       scalarv_out = ''
       do i=1,SIZE(scalarv_out,1)
          start = (i-1) * out_len + 1
          do c = 1, out_len
             scalarv_out(i)(c:c) = vector_out_total(start + c - 1)
          end do
       end do
    end if
    deallocate(vector_out_total)
    
    return
    end SUBROUTINE pgslib_Collate_char_scalar_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Collatee vectors.
  ! These routines do some checking, and then call C routines
  !          to do the work.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Collate_Int_Vector_F(Vector_out, Vector_In)

#define _DATA_TYPE_ INTEGER (PGSLib_Int_Type)
#include "collate_vector.fpp"

  END SUBROUTINE PGSLib_Collate_Int_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Collate_Real_Vector_F(Vector_out, Vector_In)

#define _DATA_TYPE_ Real (PGSLib_Real_Type)
#include "collate_vector.fpp"

  END SUBROUTINE PGSLib_Collate_Real_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Collate_Double_Vector_F(Vector_out, Vector_In)

#define _DATA_TYPE_ Real (PGSLib_Double_Type)
#include "collate_vector.fpp"

  END SUBROUTINE PGSLib_Collate_Double_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Collate_Log_Vector_F(Vector_out, Vector_In)

#define _DATA_TYPE_ logical (PGSLib_Log_Type)
#include "collate_vector.fpp"

  END SUBROUTINE PGSLib_Collate_Log_Vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE pgslib_Collate_char_vector_f(vector_out, vector_in)
    USE PGSLib_TYPE_MODULE
    USE PGSLib_Utility_MODULE,  ONLY : pgslib_check_error,     &
         &                             PGSLib_Inquire_thisPE,  &
         &                             PGSLib_Inquire_nPE,     &
         &                             PGSLib_Inquire_IO_P

    USE PGSLib_Red_Numeric_Module, ONLY:PGSLib_Global_MAXVAL,   &
         &                             PGSLib_Global_MINVAL
    implicit none
    CHARACTER (LEN = *), intent(  OUT), dimension(:):: vector_out
    CHARACTER (LEN = *), intent(IN   ), dimension(:):: vector_in
    ! Local variables
    integer string_in_len, max_string_in_len, min_string_in_len
    integer i, start, lf, c
    integer string_out_len
    integer vec_in_len
    integer (PGSLib_Int_Type), POINTER, dimension(:) :: vec_out_lengths
    logical error_flag
    CHARACTER (LEN=1), POINTER, dimension(:) :: vector_in_total, vector_out_total

    ! Check that all strings, in and out, have same length.
    ! Check in lengths
    string_in_len  = LEN(vector_in)
    max_string_in_len = PGSLib_Global_MAXVAL(string_in_len)
    min_string_in_len = PGSLib_Global_MINVAL(string_in_len)

    error_flag = (string_in_len /= max_string_in_len) .OR. &
         &       (String_in_len /= min_string_in_len)
    call pgslib_check_error(error_flag, "String lengths of all sources must be the same in PGSLib_Collate.")

    ! Check out length, but only on IO PE
    error_flag = .FALSE.
    if (PGSLib_Inquire_IO_P()) then
      string_out_len = LEN(vector_out)
       error_flag = (string_out_len /= string_in_len)
    end if
    call pgslib_check_error(error_flag, "String lengths of source and destination must be the same in PGSLib_Collate.")

    ! Setup lengths for receiving
    ALLOCATE(vec_out_lengths(PGSLib_Inquire_nPE()))
    vec_in_len = SIZE(vector_in,1) 
    call pgslib_collate(vec_out_lengths, vec_in_len)

    ! Check that destination is large enough
    error_flag = .FALSE.
    if (PGSLib_Inquire_IO_P()) then
       error_flag = SUM(vec_out_lengths) > SIZE(vector_out, 1)
    endif
    call pgslib_check_error(error_flag, "Desintation array too small in pgslib_Collate for character vector")

    ! Now do "vector" collate, each string is a vector of characters
    ! The lengths have to be adjusted to include the length of the strings
    vec_in_len  = vec_in_len  * string_in_len
    vec_out_lengths = vec_out_lengths * string_in_len
    allocate(vector_in_total(vec_in_len))
    vector_in_total = ''
    do i = 1,SIZE(vector_in,1)
       start = (i-1)*string_in_len + 1
       lf = LEN_TRIM(vector_in(i))
       do c = 1, lf
          vector_in_total(start + c -1) = vector_in(i)(c:c)
       end do
    enddo
    if (PGSLib_Inquire_IO_P()) then
       allocate(vector_out_total(SUM(vec_out_lengths)))
       vector_out_total = ''
    else
       allocate(vector_out_total(1))
    end if

    call pgslib_collate_vector_c(vector_out_total, vec_out_lengths, vector_in_total, vec_in_len)

    deallocate(vector_in_total)
    if (PGSLib_Inquire_IO_P()) then
       vector_out = ''
       do i=1,SIZE(vector_out,1)
          start = (i-1) * string_out_len + 1
          do c = 1, string_out_len
             vector_out(i)(c:c) = vector_out_total(start + c -1)
          end do
       enddo
    end if
    deallocate(vector_out_total)

    return
  end SUBROUTINE pgslib_Collate_char_vector_f

END MODULE PGSLib_IO_Collate_MODULE
