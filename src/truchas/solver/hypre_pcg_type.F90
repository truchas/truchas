!!
!! HYPRE_PCG_TYPE
!!
!! This module defines a preconditioned CG solver class built on the
!! ParCSR PCG solver from Hypre that uses BoomerAMG as the preconditioner.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module hypre_pcg_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fhypre
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  type, public :: hypre_pcg
    private
    type(pcsr_matrix), pointer :: A => null()
    integer :: nrows = 0, ilower = 0, iupper = 0
    integer, allocatable :: rows(:)
    type(hypre_obj) :: solver = hypre_null_obj ! HYPRE_Solver object handle
    type(hypre_obj) :: precon = hypre_null_obj ! HYPRE_Solver object handle
    type(hypre_obj) :: Ah = hypre_null_obj     ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: bh = hypre_null_obj     ! HYPRE_IJVector object handles
    type(hypre_obj) :: xh = hypre_null_obj     ! HYPRE_IJVector object handles
    type(parameter_list), pointer :: params    ! reference only -- do not own
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: matrix
    procedure :: get_metrics
    final :: hypre_pcg_delete
  end type hypre_pcg

contains

  !! Final subroutine for HYPRE_PCG objects
  subroutine hypre_pcg_delete (this)
    type(hypre_pcg), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy (this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy (this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy (this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_PCGDestroy (this%solver, ierr)
    if (hypre_associated(this%precon)) call fHYPRE_BoomerAMGDestroy (this%precon, ierr)
    INSIST(ierr == 0)
  end subroutine hypre_pcg_delete

  subroutine init (this, A, params)

    class(hypre_pcg), intent(out) :: this
    type(pcsr_matrix), intent(in), target :: A
    type(parameter_list), pointer, intent(in) :: params

    integer :: ierr, i

    this%A => A
    this%params => params

    this%nrows  = A%graph%row_imap%onp_size
    this%ilower = A%graph%row_imap%first_gid
    this%iupper = A%graph%row_imap%last_gid
    allocate(this%rows(this%nrows))
    this%rows = [(i, i = this%ilower, this%iupper)]

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%bh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%xh, 0, ierr)
    INSIST(ierr == 0)

  end subroutine init

  function matrix (this)
    class(hypre_pcg), intent(in) :: this
    type(pcsr_matrix), pointer :: matrix
    matrix => this%A
  end function matrix

  subroutine get_metrics (this, num_itr)
    class(hypre_pcg), intent(in) :: this
    integer, intent(out), optional :: num_itr
    integer :: ierr
    INSIST(hypre_associated(this%solver))
    if (present(num_itr)) then
      call fHYPRE_PCGGetNumIterations (this%solver, num_itr, ierr)
      INSIST(ierr == 0)
    end if
  end subroutine get_metrics

  subroutine setup (this)

    class(hypre_pcg), intent(inout) :: this

    integer :: ierr, ipar
    real(r8) :: dummy(this%nrows), rpar

    ASSERT(associated(this%A))
    call this%A%copy_to_ijmatrix(this%Ah)

    !! Create the Hypre BoomerAMG solver object (preconditioner).  Note that
    !! once the solver has been setup, it is not possible to change the matrix
    !! values without completely destroying the solver and recreating it.
    if (hypre_associated(this%precon)) call fHYPRE_BoomerAMGDestroy (this%precon, ierr)
    call fHYPRE_BoomerAMGCreate (this%precon, ierr)
    call fHYPRE_BoomerAMGSetOldDefault(this%precon, ierr)
    INSIST(ierr == 0)

    !! Create the Hypre PCG solver object. This supposes that once the solver
    !! has been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it.  VERIFY THIS.
    if (hypre_associated(this%solver)) call fHYPRE_PCGDestroy (this%solver, ierr)
    call fHYPRE_PCGCreate (this%solver, ierr)
    INSIST(ierr == 0)

    !! Print level (optional, default is no output)
    call this%params%get ('print-level', ipar, default=0)
    INSIST(ipar >= 0 .and. ipar <= 3) !TODO: replace with proper error handling
    call fHYPRE_PCGSetPrintLevel (this%solver, ipar, ierr)
    INSIST(ierr == 0)
    call this%params%get ('amg-print-level', ipar, default=0)
    INSIST(ipar >= 0) !TODO: replace with proper error handling
    call fHYPRE_BoomerAMGSetPrintLevel  (this%precon, ipar, ierr)
    INSIST(ierr == 0)

    !! Logging level (optional, default is none)
    call this%params%get ('logging-level', ipar, default=0)
    INSIST(ipar >= 0) !TODO: replace with proper error handling
    call fHYPRE_BoomerAMGSetLogging (this%precon, ipar, ierr)
    INSIST(ierr == 0)

    !! PCG solver relative tolerance (required)
    call this%params%get ('rel-tol', rpar)
    INSIST(rpar >= 0.0_r8) !TODO: replace with proper error handling
    call fHYPRE_PCGSetTol (this%solver, rpar, ierr)
    INSIST(ierr == 0)

    !! Absolute error tolerance for PCG solver (optional, use Hypre default -- none)
    if (this%params%is_parameter ('abs-tol')) then
      call this%params%get ('abs-tol', rpar)
      INSIST(rpar > 0.0_r8) !TODO: replace with proper error handling
      call fHYPRE_PCGSetAbsoluteTol (this%solver, rpar, ierr)
      INSIST(ierr == 0)
    end if

    !! Use the two norm instead of the energy norm (hardwired).  Why?
    call fHYPRE_PCGSetTwoNorm (this%solver, 1, ierr)
    INSIST(ierr == 0)

    !! Maximum number of PCG iterations (required)
    call this%params%get ('max-iter', ipar)
    INSIST(ipar > 0) !TODO: replace with proper error handling
    call fHYPRE_PCGSetTol (this%solver, rpar, ierr)
    call fHYPRE_PCGSetMaxIter (this%solver, ipar, ierr)
    INSIST(ierr == 0)

    !! Number of AMG cycles (using as a preconditioner, optional, default=1)
    call this%params%get ('amg-num-cycles', ipar, default=1)
    INSIST(ipar >= 1) !TODO: replace with proper error handling
    call fHYPRE_BoomerAMGSetMaxIter (this%precon, ipar, ierr)
    call fHYPRE_BoomerAMGSetTol (this%precon, 0.0_r8, ierr)
    INSIST(ierr == 0)

    !! AMG strength threshold (default 0.5 (3D Laplace))
    call this%params%get ('amg-strong-threshold', rpar, default=0.5_r8)
    INSIST(rpar > 0.0_r8 .and. rpar < 1.0_r8) !TODO: replace with proper error handling
    call fHYPRE_BoomerAMGSetStrongThreshold  (this%precon, rpar, ierr)
    INSIST(ierr == 0)

    !! Maximum number of AMG levels (optional, use Hypre default)
    if (this%params%is_parameter('amg-max-levels')) then
      call this%params%get ('amg-max-levels', ipar)
      INSIST(ipar > 1) !TODO: replace with proper error handling
      call fHYPRE_BoomerAMGSetMaxLevels   (this%precon, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! The AMG coarsening method (optional, use Hypre default)
    if (this%params%is_parameter('amg-coarsen-method')) then
      call this%params%get ('amg-coarsen-type', ipar)
      INSIST(ipar >= 0) !TODO: replace with proper error handling
      call fHYPRE_BoomerAMGSetCoarsenType (this%precon, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! Number of AMG smoother sweeps (optional, default 1)
    call this%params%get ('amg-smoothing-sweeps', ipar, default=1)
    INSIST(ipar >= 0) !TODO: replace with proper error handling
    call fHYPRE_BoomerAMGSetNumSweeps (this%precon, ipar, ierr)
    INSIST(ierr == 0)

    !! The AMG smoothing method (optional, use Hypre default)
    if (this%params%is_parameter('amg-smoothing-method')) then
      call this%params%get ('amg-smoothing-method', ipar)
      INSIST(ipar >= 0) !TODO: replace with proper error handling
      call fHYPRE_BoomerAMGSetRelaxType (this%precon, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! The AMG interpolation method (optional, use Hypre default)
    if (this%params%is_parameter('amg-interp-method')) then
      call this%params%get ('amg-interp-method', ipar)
      INSIST(ipar >= 0) !TODO: replace with proper error handling
      call fHYPRE_BoomerAMGSetInterpType (this%solver, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! Hardwire a few AMG parameters.
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 3, 1, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 4, 2, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 9, 3, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre solution and RHS vectors;
    !! structure (not values) needed to setup PCG solver.
    dummy = 0.0_r8
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, this%rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, this%rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.
    call fHYPRE_PCGSetPrecond (this%solver, this%precon, ierr)
    call fHYPRE_PCGSetup (this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

  end subroutine setup

  subroutine solve (this, b, x, stat)

    class(hypre_pcg), intent(inout) :: this
    real(r8), intent(in) :: b(:)
    real(r8), intent(inout) :: x(:)
    integer, intent(out) :: stat

    integer :: ierr

    call fHYPRE_ClearAllErrors

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, this%rows, b, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, this%rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! Solve the system.
    call fHYPRE_PCGSolve (this%solver, this%Ah, this%bh, this%xh, ierr)
    if (ierr /= 0) then
      INSIST(ierr == HYPRE_ERROR_CONV)
      stat = 1
      call fHYPRE_ClearAllErrors
    else
      stat = 0
    end if

    !! Retrieve the solution vector from HYPRE.
    call fHYPRE_IJVectorGetValues (this%xh, this%nrows, this%rows, x, ierr)
    INSIST(ierr == 0)

  end subroutine solve

end module hypre_pcg_type
