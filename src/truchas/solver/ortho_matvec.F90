!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ORTHO_MATVEC
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures necessary to perform a matrix-vector multiply
  !   (Ax), returning it in vector y
  !
  !   Public Interface:
  !
  !     * call ORTHO_MATVEC_Y_EQ_AX (nunk, nfc, Matrix, X, X_Neighbors,
  !                                  Y, status)
  !
  ! Contains: ORTHO_MATVEC_Y_EQ_AX
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_timers
  Implicit None
  Private

  Public :: ORTHO_MATVEC_Y_EQ_AX

CONTAINS

  SUBROUTINE ORTHO_MATVEC_Y_EQ_AX (nunk, nfc, Matrix, X, X_Neighbors, Y, status)
    !=======================================================================
    ! Purpose:
    !
    !   Compute y = Ax where A is an ortho matrix stored in ELL format
    !
    !=======================================================================

    ! Arguments
    integer, intent(IN) :: nunk, nfc
    real(r8), dimension(nunk), intent(IN) :: X
    real(r8), dimension(0:nfc,nunk), intent(IN) :: Matrix  
    real(r8), dimension(nfc,nunk), intent(IN) :: X_Neighbors
    integer, intent(OUT) :: status
    real(r8), dimension(nunk), intent(INOUT) :: Y

    ! Local Variables
    integer :: f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the timer.
    call start_timer("Solver TMP2")

    ! Perform the matrix-vector multiply
    Y = Matrix(0,:) * X
    do f = 1,nfc
       Y = Y + Matrix(f,:) * X_Neighbors(f,:)
    end do

    ! Stop the timer.
    call stop_timer ("Solver TMP2")

    status = 0

  END SUBROUTINE ORTHO_MATVEC_Y_EQ_AX

END MODULE ORTHO_MATVEC
