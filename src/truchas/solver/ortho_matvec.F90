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
  use kind_module, only: int_kind, real_kind

  Implicit None

  Private

  ! public procedures
  Public :: ORTHO_MATVEC_Y_EQ_AX

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE ORTHO_MATVEC_Y_EQ_AX (nunk, nfc, Matrix, X, X_Neighbors, Y, status)
    !=======================================================================
    ! Purpose:
    !
    !   Compute y = Ax where A is an ortho matrix stored in ELL format
    !
    !=======================================================================
    use timing_tree

    implicit none

    ! Arguments
    integer(KIND = int_kind),                      intent(IN)    :: nunk, nfc
    real(KIND = real_kind), dimension(nunk),       intent(IN)    :: X
    real(KIND = real_kind), dimension(0:nfc,nunk), intent(IN)    :: Matrix  
    real(KIND = real_kind), dimension(nfc,nunk),   intent(IN)    :: X_Neighbors
    integer(KIND = int_kind),                      intent(OUT)   :: status
    real(KIND = real_kind), dimension(nunk),       intent(INOUT) :: Y

    ! Local Variables
    integer(KIND = int_kind) :: f

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
    return

  END SUBROUTINE ORTHO_MATVEC_Y_EQ_AX

END MODULE ORTHO_MATVEC
