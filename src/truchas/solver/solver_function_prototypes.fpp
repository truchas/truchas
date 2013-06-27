  ! -*- Mode: f90 -*-

  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  ! define interfaces for preconditioners, matrix-vector multiplies,
  ! NK dampers, residuals, etc., that are visible in multiple places
  !
  ! this file should be included everywhere that the prototyped
  ! procedures are used
  !
  ! these prototypes must be kept consistent with
  !
  !    1) the precondition and y_eq_ax routines
  !    2) the interfaces used by the linear solvers (particularly UbikSolve)
  !
  ! contains: interface prototypes for
  !    MATVEC
  !    PRECONDITIONER
  !    PRECONDITIONER_UPDATE
  !    RESIDUAL
  !
  ! Author(s): ferrell@cpca.com
  !-----------------------------------------------------------------------------

    interface
       subroutine MATVEC (X, Y, status)
          use kinds, only: r8
          use UbikSolve_module, only: Ubik_vector_type
          type(Ubik_vector_type), intent(INOUT) :: X
          real(r8), target, intent(INOUT) :: Y(:)
          integer, intent(OUT) :: status
       end subroutine MATVEC
    end interface

    interface
       subroutine PRECONDITIONER (B, X, status)
          use kinds, only: r8
          use UbikSolve_module, only: Ubik_vector_type
          real(r8), target, intent(IN) :: B(:)
          type(Ubik_vector_type), intent(INOUT) :: X
          integer, intent(OUT) :: status
       end subroutine PRECONDITIONER
    end interface

    interface
       subroutine PRECONDITIONER_UPDATE (x)
          use kinds, only: r8
          real(r8), intent(IN) :: x(:)
       end subroutine PRECONDITIONER_UPDATE
    end interface

    interface
       subroutine RESIDUAL (x_old, x, r)
          use kinds, only: r8
          real(r8), intent(IN)  :: x_old(:)
          real(r8), intent(IN)  :: x(:)
          real(r8), intent(OUT) :: r(:)
       end subroutine RESIDUAL
    end interface
