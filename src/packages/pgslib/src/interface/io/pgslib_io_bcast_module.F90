MODULE PGSLib_IO_BCast_MODULE
  USE PGSLib_Type_MODULE
  use pgslib_c_binding
  use,intrinsic :: iso_fortran_env, only: int8
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_BCast

  !  The routines supported in this module are
  !  PGSLib_BCast
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

  ! $Id: pgslib_io_bcast_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PGSLib_BCast Interfaces  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE PGSLib_bcast

     ! Broadcast Scalars
     ! These routines provide F90 interfaces to the C routines

     MODULE PROCEDURE PGSLib_BCast_Int8_Scalar_F
     MODULE PROCEDURE PGSLib_BCast_Int_Scalar_F
     MODULE PROCEDURE PGSLib_BCast_REAL_Scalar_F
     MODULE PROCEDURE PGSLib_BCast_Double_Scalar_F
     MODULE PROCEDURE PGSLib_BCast_Log_Scalar_F

     ! Broadcast Vectors
     ! Interfaces for the F90 interface routines, which take a single vector argument
     MODULE Procedure PGSLib_BCast_Int8_Vector_F
     MODULE Procedure PGSLib_BCast_Int_Vector_F
     MODULE Procedure PGSLib_BCast_Real_Vector_F
     MODULE Procedure PGSLib_BCast_Double_Vector_F
     MODULE Procedure PGSLib_BCast_Log_Vector_F
     MODULE Procedure PGSLib_BCast_Char_Vector_F

     ! Broadcast 2D Vectors
     ! Interfaces for the F90 interface routines, which take a single vector argument
     MODULE Procedure PGS_BCast_Int8_vector_2d_F
     MODULE Procedure PGS_BCast_Int_vector_2d_F
     MODULE Procedure PGS_BCast_Real_vector_2d_F
     MODULE Procedure PGS_BCast_Double_vector_2d_F
     MODULE Procedure PGS_BCast_Log_vector_2d_F
     MODULE Procedure PGS_BCast_Char_vector_2d_F

     ! Broadcast 3D Vectors
     ! Interfaces for the F90 interface routines, which take a single vector argument
     MODULE Procedure PGS_BCast_Int8_vector_3d_F
     MODULE Procedure PGS_BCast_Int_vector_3d_F
     MODULE Procedure PGS_BCast_Real_vector_3d_F
     MODULE Procedure PGS_BCast_Double_vector_3d_F
     MODULE Procedure PGS_BCast_Log_vector_3d_F
     MODULE Procedure PGS_BCast_Char_vector_3d_F

  END INTERFACE!For the user callable PGSLib_BCast

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! PGSLib_BCast Subroutine Bodies !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Broadcast of Scalars
  ! These routines only call C routines, they exist only for symmetry
  ! with the vector versions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_int8_scalar_F(scalar)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (int8) :: scalar
    Call PGSLib_BCast_int8_Scalar_c(scalar)
    RETURN
  END SUBROUTINE PGSLib_bcast_int8_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_int_scalar_F(scalar)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (KIND=PGSLib_INT_TYPE) :: scalar
    Call PGSLib_BCast_int_Scalar_c(scalar)
    RETURN
  END SUBROUTINE PGSLib_bcast_int_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_real_scalar_F(scalar)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    REAL (KIND=PGSLib_REAL_TYPE) ::scalar
    Call PGSLib_BCast_float_Scalar_c(scalar)
    RETURN
  END SUBROUTINE PGSLib_bcast_real_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_double_scalar_F(scalar)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    REAL (KIND=PGSLib_DOUBLE_TYPE) :: scalar
    Call PGSLib_BCast_double_Scalar_c(scalar)
    RETURN
  END SUBROUTINE PGSLib_bcast_double_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_log_scalar_F(scalar)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    LOGICAL (PGSLib_Log_Type) :: scalar
    Call PGSLib_BCast_log_Scalar_c(scalar)
    RETURN
  END SUBROUTINE PGSLib_bcast_log_scalar_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Broadcast of Vectors
  ! These routines do some checking, and then call C routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_int8_vector_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (int8) ,dimension(:):: vector

    ! Local variables
    integer vec_len
    vec_len = SIZE(vector)
    call pgslib_bcast_int8_vector_c(vector, vec_len)
    return

  END SUBROUTINE PGSLib_bcast_int8_vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_int_vector_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (KIND=PGSLib_INT_TYPE) ,dimension(:):: vector

    ! Local variables
    integer vec_len
    vec_len = SIZE(vector)
    call pgslib_bcast_int_vector_c(vector, vec_len)
    return

  END SUBROUTINE PGSLib_bcast_int_vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_real_vector_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    REAL (KIND=PGSLib_REAL_TYPE) ,dimension(:) ::vector

    ! Local variables
    integer vec_len
    vec_len = SIZE(vector)
    call pgslib_bcast_float_vector_c (vector, vec_len)
    return

  END SUBROUTINE PGSLib_bcast_real_vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_double_vector_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    REAL (KIND=PGSLib_DOUBLE_TYPE) ,dimension(:):: vector

    ! Local variables
    integer vec_len
    vec_len = SIZE(vector)
    call pgslib_bcast_double_vector_c (vector, vec_len)
    return

  END SUBROUTINE PGSLib_bcast_double_vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_log_vector_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    LOGICAL (PGSLib_Log_Type) ,dimension(:):: vector

    ! Local variables
    integer vec_len
    vec_len = SIZE(vector)
    call pgslib_bcast_log_vector_c (vector, vec_len)
    return

  END SUBROUTINE PGSLib_bcast_log_vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_bcast_char_vector_F(vector)
    USE PGSLib_TYPE_MODULE
    USE PGSLib_Utility_MODULE, ONLY : PGSLib_Inquire_IO_P
    IMPLICIT NONE
    CHARACTER (LEN=*) :: vector

    ! Local variables
    integer vec_len

    vec_len = len_trim(vector)
    call pgslib_bcast(vec_len)

    ! We need to overwrite anything that might be in the character string 
    ! on non-root pe''s.
    if (.NOT.PGSLib_Inquire_IO_P()) then
       vector = " "
    endif

    call pgslib_bcast_char_vector_c (vector, vec_len)

    return
  END SUBROUTINE PGSLib_bcast_char_vector_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Broadcast of 2D Vectors
  ! These routines do some checking, and then call C routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_int8_vector_2d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (int8), intent(INOUT), dimension(:,:):: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:))
    enddo

    return

  END SUBROUTINE PGS_bcast_int8_vector_2d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_int_vector_2d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (KIND=PGSLib_INT_TYPE), intent(INOUT), dimension(:,:):: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:))
    enddo

    return

  END SUBROUTINE PGS_bcast_int_vector_2d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_real_vector_2d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    real (PGSLib_Real_Type), intent(INOUT), dimension(:,:) :: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:))
    enddo

    RETURN
  END SUBROUTINE PGS_bcast_real_vector_2d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_double_vector_2d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    real (PGSLib_Double_Type), intent(INOUT), dimension(:,:) :: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:))
    enddo

    RETURN
  END SUBROUTINE PGS_bcast_double_vector_2d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_log_vector_2d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    logical (PGSLib_Log_Type), intent(INOUT), dimension(:,:):: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:))
    enddo

    return

  END SUBROUTINE PGS_bcast_log_vector_2d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_char_vector_2d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    CHARACTER (LEN=*), dimension(:) :: vector

    ! Local variables
    integer dim1, idim

    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim))
    enddo

    RETURN
  END SUBROUTINE PGS_bcast_char_vector_2d_F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Broadcast of 3D Vectors
  ! These routines do some checking, and then call C routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_int8_vector_3d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (int8), intent(INOUT), dimension(:,:,:):: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:,:))
    enddo

    return

  END SUBROUTINE PGS_bcast_int8_vector_3d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_int_vector_3d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    INTEGER (KIND=PGSLib_INT_TYPE), intent(INOUT), dimension(:,:,:):: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:,:))
    enddo

    return

  END SUBROUTINE PGS_bcast_int_vector_3d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_real_vector_3d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    real (PGSLib_Real_Type), intent(INOUT), dimension(:,:,:) :: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:,:))
    enddo

    RETURN
  END SUBROUTINE PGS_bcast_real_vector_3d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_double_vector_3d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    real (PGSLib_Double_Type), intent(INOUT), dimension(:,:,:) :: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:,:))
    enddo

    RETURN
  END SUBROUTINE PGS_bcast_double_vector_3d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_log_vector_3d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    logical (PGSLib_Log_Type), intent(INOUT), dimension(:,:,:):: vector

    ! Local variables
    integer dim1, idim
    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:,:))
    enddo

    return

  END SUBROUTINE PGS_bcast_log_vector_3d_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGS_bcast_char_vector_3d_F(vector)
    USE PGSLib_TYPE_MODULE
    IMPLICIT NONE
    CHARACTER (LEN=*), dimension(:,:) :: vector

    ! Local variables
    integer dim1, idim, vec_len

    dim1 = SIZE(vector,1)

    do idim = 1, dim1
       call pgslib_bcast(vector(idim,:))
    enddo

    RETURN
  END SUBROUTINE PGS_bcast_char_vector_3d_F

END MODULE PGSLib_IO_BCast_MODULE
