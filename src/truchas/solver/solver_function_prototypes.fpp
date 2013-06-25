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
          use kind_module, only: int_kind, real_kind
          use UbikSolve_module
          type (Ubik_vector_type),                intent(INOUT) :: X
          real (real_kind), dimension(:), target, intent(INOUT) :: Y
          integer (int_kind),                     intent(OUT)   :: status
       end subroutine MATVEC
    end interface

    interface
       subroutine PRECONDITIONER (B, X, status)
          use kind_module, only: int_kind, real_kind
          use UbikSolve_module
          real (real_kind), dimension(:), target, intent(IN)    :: B
          type (Ubik_vector_type),                intent(INOUT) :: X
          integer (int_kind),                     intent(OUT)   :: status
       end subroutine PRECONDITIONER
    end interface

    interface
       subroutine PRECONDITIONER_UPDATE (x)
          use kind_module, only: real_kind
          real (real_kind), dimension(:), intent(IN) :: x
       end subroutine PRECONDITIONER_UPDATE
    end interface

    interface
       subroutine RESIDUAL (x_old, x, r)
          use kind_module,                only: real_kind
          real (real_kind), dimension(:), intent(IN)  :: x_old
          real (real_kind), dimension(:), intent(IN)  :: x
          real (real_kind), dimension(:), intent(OUT) :: r
       end subroutine RESIDUAL
    end interface
