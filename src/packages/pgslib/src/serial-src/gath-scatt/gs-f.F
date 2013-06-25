! Support routines for gath-scat directory
! In the parallel code these are C routines.  In the
! serial code these are fortran routines.

! $Id: gs-f.F,v 1.2 2001/03/22 00:26:13 ferrell Exp $

subroutine PGSLib_GS_Init_Trace_C(GS_Trace)
  implicit none
  integer, dimension(*):: GS_Trace

  GS_Trace(1) =  1
  GS_Trace(2) = -1
  RETURN
END subroutine PGSLib_GS_Init_Trace_C

subroutine PGSLib_GS_Release_Trace_C(GS_Trace)
  implicit none
  integer, dimension(*):: GS_Trace

  RETURN
END subroutine PGSLib_GS_Release_Trace_C

subroutine PGSLib_Setup_N_Duplicate_C(N_Duplicate, GS_Trace)
  implicit none
  integer :: N_Duplicate
  integer, dimension(*):: GS_Trace

  N_Duplicate = 0
  RETURN
END subroutine PGSLib_Setup_N_Duplicate_C

subroutine PGSLib_Prep_Supplement_C(N_Supplement, PEs, GS_Trace)
  implicit none
  integer :: N_Supplement
  integer, dimension(*) :: PEs
  integer, dimension(*):: GS_Trace

  GS_Trace(1) = 0
  RETURN
END subroutine PGSLib_Prep_Supplement_C

subroutine PGSLib_Setup_Duplicate_Buffer_C(GS_Trace)
  implicit none
  integer, dimension(*):: GS_Trace
  
  RETURN
END subroutine PGSLib_Setup_Duplicate_Buffer_C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gather routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PGSLib_Gather_Buf_INT_C( Supplement, &
         &                     Duplicate, BlockSize, GS_Trace)
  implicit none 
  integer, dimension(*):: Supplement
  integer, dimension(*):: Duplicate
  integer :: BlockSize
  integer, dimension(*):: GS_Trace

  ! For the serial code, there shouldn't be any off-PE moves, so don't need to do anything.
  RETURN
End subroutine PGSLib_Gather_Buf_INT_C

subroutine PGSLib_Gather_Buf_FLOAT_C( Supplement, &
         &                     Duplicate, BlockSize, GS_Trace)
  implicit none 
  REAL, dimension(*):: Supplement
  REAL, dimension(*):: Duplicate
  integer :: BlockSize
  integer, dimension(*):: GS_Trace

  ! For the serial code, there shouldn't be any off-PE moves, so don't need to do anything.
  RETURN
End subroutine PGSLib_Gather_Buf_FLOAT_C

subroutine PGSLib_Gather_Buf_DOUBLE_C( Supplement, &
         &                     Duplicate, BlockSize, GS_Trace)
  implicit none 
  REAL , dimension(*):: Supplement
  REAL , dimension(*):: Duplicate
  integer :: BlockSize
  integer, dimension(*):: GS_Trace

  ! For the serial code, there shouldn't be any off-PE moves, so don't need to do anything.
  RETURN
End subroutine PGSLib_Gather_Buf_DOUBLE_C

subroutine PGSLib_Gather_Buf_LOG_C( Supplement, &
         &                     Duplicate, BlockSize, GS_Trace)
  implicit none 
  logical , dimension(*):: Supplement
  logical , dimension(*):: Duplicate
  integer :: BlockSize
  integer, dimension(*):: GS_Trace

  ! For the serial code, there shouldn't be any off-PE moves, so don't need to do anything.
  RETURN
End subroutine PGSLib_Gather_Buf_LOG_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scatter Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PGSLib_Scatter_Buf_INT_C(  Duplicate, &
         &                 Supplement, GS_Trace)
  implicit none
  INTEGER, dimension(*):: Duplicate
  INTEGER, dimension(*):: Supplement
  integer, dimension(*):: GS_Trace

  ! For serial code there shouldn't be any off-PE moves, so don't need to do anything
  RETURN
END subroutine PGSLib_Scatter_Buf_INT_C

subroutine PGSLib_Scatter_Buf_FLOAT_C(  Duplicate, &
         &                 Supplement, GS_Trace)
  implicit none
  REAL, dimension(*):: Duplicate
  REAL, dimension(*):: Supplement
  integer, dimension(*):: GS_Trace

  ! For serial code there shouldn't be any off-PE moves, so don't need to do anything
  RETURN
END subroutine PGSLib_Scatter_Buf_FLOAT_C

subroutine PGSLib_Scatter_Buf_DOUBLE_C(  Duplicate, &
         &                 Supplement, GS_Trace)
  implicit none
  REAL , dimension(*):: Duplicate
  REAL , dimension(*):: Supplement
  integer, dimension(*):: GS_Trace

  ! For serial code there shouldn't be any off-PE moves, so don't need to do anything
  RETURN
END subroutine PGSLib_Scatter_Buf_DOUBLE_C


subroutine PGSLib_Scatter_Buf_LOG_C(  Duplicate, &
         &                 Supplement, GS_Trace)
  implicit none
  logical , dimension(*):: Duplicate
  logical , dimension(*):: Supplement
  integer, dimension(*):: GS_Trace

  ! For serial code there shouldn't be any off-PE moves, so don't need to do anything
  RETURN
END subroutine PGSLib_Scatter_Buf_LOG_C

subroutine PGSLib_Trace_Degree_C(Scatter_Degree, Gather_Degree, GS_Trace)
  implicit none
  integer :: Scatter_Degree, Gather_Degree
  integer, dimension(*):: GS_Trace

  ! For serial code there should not be any off-PE moves, so  return 0
  Scatter_Degree = 0
  Gather_Degree = 0
  RETURN
END subroutine PGSLib_Trace_Degree_C
