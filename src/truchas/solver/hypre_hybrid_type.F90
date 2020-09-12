!!
!! HYPRE_HYBRID_TYPE
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

module hypre_hybrid_type

  use kinds, only: r8
  use fhypre
  use pcsr_matrix_type
  use index_partitioning
  use parameter_list_type
  implicit none
  private

  type, public :: hypre_hybrid
    private
    type(pcsr_matrix), pointer :: A => null()
    integer :: nrows = 0, ilower = 0, iupper = 0
    type(hypre_obj) :: solver = hypre_null_obj ! HYPRE_Solver object handle
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
    procedure :: metrics_string
    final :: hypre_hybrid_delete
  end type hypre_hybrid

contains

  !! Final subroutine for HYPRE_HYBRID objects.
  subroutine hypre_hybrid_delete (this)
    type(hypre_hybrid) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy (this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy (this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy (this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_ParCSRHybridDestroy (this%solver, ierr)
    INSIST(ierr == 0)
    call fHYPRE_ClearAllErrors
  end subroutine hypre_hybrid_delete

  subroutine init (this, A, params)

    class(hypre_hybrid), intent(out) :: this
    type(pcsr_matrix), intent(in), target :: A
    type(parameter_list), pointer, intent(in) :: params

    integer :: ierr

    this%A => A
    this%params => params

    this%nrows  = A%graph%row_ip%onP_size()    ! number of on-process rows (if parallel)
    this%ilower = A%graph%row_ip%first_index() ! global index of first on-process row (if parallel)
    this%iupper = A%graph%row_ip%last_index()  ! global index of last on-process row (if parallel)

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%bh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%xh, 0, ierr)
    INSIST(ierr == 0)

  end subroutine init

  function matrix (this)
    class(hypre_hybrid), intent(in) :: this
    type(pcsr_matrix), pointer :: matrix
    matrix => this%A
  end function matrix

  subroutine get_metrics (this, num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
    class(hypre_hybrid), intent(in) :: this
    integer, intent(out), optional :: num_itr, num_dscg_itr, num_pcg_itr
    real(r8), intent(out), optional :: rel_res_norm
    integer :: ierr
    INSIST(hypre_associated(this%solver))
    if (present(num_itr)) then
      call fHYPRE_ParCSRHybridGetNumIterations (this%solver, num_itr, ierr)
      INSIST(ierr == 0)
    end if
    if (present(num_dscg_itr)) then
      call fHYPRE_ParCSRHybridGetDSCGNumIterations (this%solver, num_dscg_itr, ierr)
      INSIST(ierr == 0)
    end if
    if (present(num_pcg_itr)) then
      call fHYPRE_ParCSRHybridGetPCGNumIterations (this%solver, num_pcg_itr, ierr)
      INSIST(ierr == 0)
    end if
    if (present(rel_res_norm)) then
      call fHYPRE_ParCSRHybridGetFinalRelativeResidualNorm (this%solver, rel_res_norm, ierr)
      INSIST(ierr == 0)
    end if
  end subroutine get_metrics

  function metrics_string(this) result(string)
    class(hypre_hybrid), intent(in) :: this
    character(:), allocatable :: string
    character(80) :: buffer
    integer :: num_itr, num_dscg_itr, num_pcg_itr
    real(r8) :: rel_res_norm
    call get_metrics(this, num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
    write(buffer,'(i4," (DS), ",i4," (AMG), ",es10.4," (|r|/|b|)")') &
        num_dscg_itr, num_pcg_itr, rel_res_norm
    string = trim(buffer)
  end function metrics_string

  subroutine setup (this)

    class(hypre_hybrid), intent(inout) :: this

    integer :: ierr, i, ipar
    integer :: rows(this%nrows)
    real(r8) :: dummy(this%nrows), rpar
    character(:), allocatable :: cpar
    logical :: lpar

    !! Create the Hypre ParCSR Hybrid solver object. This supposes that once the
    !! solver has been setup, it is not possible to change the matrix values
    !! without completely destroying the solver and recreating it.  VERIFY THIS.
    if (hypre_associated(this%solver)) call fHYPRE_ParCSRHybridDestroy (this%solver, ierr)
    call fHYPRE_ParCSRHybridCreate (this%solver, ierr)
    INSIST(ierr == 0)

    !! Krylov solver relative tolerance for Krylov solver (required)
    call this%params%get ('rel-tol', rpar)
    INSIST(rpar >= 0.0_r8)  !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetTol (this%solver, rpar, ierr)
    INSIST(ierr == 0)

    !! Krylov solver absolute tolerance (optional, default none (=0))
    call this%params%get ('abs-tol', rpar, default=0.0_r8)
    INSIST(rpar >= 0.0_r8)  !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetAbsoluteTol (this%solver, rpar, ierr)
    INSIST(ierr == 0)

    !! Convergence rate tolerance (optional, use Hypre default)
    if (this%params%is_parameter('conv-rate-tol')) then
      call this%params%get ('conv-rate-tol', rpar)
      INSIST(rpar > 0.0_r8 .and. rpar < 1.0_r8) !TODO: replace with proper error handling
      call fHYPRE_ParCSRHybridSetConvergenceTol (this%solver, rpar, ierr)
      INSIST(ierr == 0)
    end if

    !! Max number of diagonally-scaled Krylov iterations (required)
    call this%params%get ('max-ds-iter', ipar)
    INSIST(ipar > 0)  !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetDSCGMaxIter (this%solver, ipar, ierr)
    INSIST(ierr == 0)

    !! Max number of AMG-preconditioned Krylov iterations (required)
    call this%params%get ('max-amg-iter', ipar)
    INSIST(ipar > 0)  !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetPCGMaxIter (this%solver, ipar, ierr)
    INSIST(ierr == 0)

    !! Kyrlov method (required)
    call this%params%get ('krylov-method', cpar, default='cg')
    select case (cpar)
    case ('cg')
      call fHYPRE_ParCSRHybridSetSolverType (this%solver, 1, ierr)
    case ('gmres')
      call fHYPRE_ParCSRHybridSetSolverType (this%solver, 2, ierr)
    case ('bicgstab')
      call fHYPRE_ParCSRHybridSetSolverType (this%solver, 3, ierr)
    case default
      INSIST(.false.) !TODO: replace with proper error handling
    end select
    INSIST(ierr == 0)

    !! Krylov space dimension for restarted GMRES (required for GMRES)
    if (cpar == 'gmres') then
      call this%params%get ('gmres-krylov-dim', ipar)
      INSIST(ipar > 0)  !TODO: replace with proper error handling
      call fHYPRE_ParCSRHybridSetKDim (this%solver, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! Two norm for CG (optional, use Hypre default)
    if (cpar == 'cg') then ! TODO: should this apply to bicgstab too?
      if (this%params%is_parameter('cg-use-two-norm')) then
        call this%params%get ('cg-use-two-norm', lpar)
        if (lpar) then
          call fHYPRE_ParCSRHybridSetTwoNorm (this%solver, 1, ierr)
        else
          call fHYPRE_ParCSRHybridSetTwoNorm (this%solver, 0, ierr)
        end if
        INSIST(ierr == 0)
      end if
    end if

    !! Logging level (optional, default is none)
    call this%params%get ('logging-level', ipar, default=1)
    INSIST(ipar >= 0) !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetLogging (this%solver, ipar, ierr)
    INSIST(ierr == 0)

    !! Print level (optional, default is no output)
    call this%params%get ('print-level', ipar, default=0)
    INSIST(ipar >= 0) !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetPrintLevel (this%solver, ipar, ierr)
    INSIST(ierr == 0)

    !! AMG strength threshold (default 0.5 (3D Laplace))
    call this%params%get ('amg-strong-threshold', rpar, default=0.5_r8)
    INSIST(rpar > 0.0_r8 .and. rpar < 1.0_r8) !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetStrongThreshold (this%solver, rpar, ierr)
    INSIST(ierr == 0)

    !! Maximum number of AMG levels (optional, use Hypre default)
    if (this%params%is_parameter('amg-max-levels')) then
      call this%params%get ('amg-max-levels', ipar)
      INSIST(ipar > 1) !TODO: replace with proper error handling
      call fHYPRE_ParCSRHybridSetMaxLevels (this%solver, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! The AMG coarsening method (optional, use Hypre default)
    if (this%params%is_parameter('amg-coarsen-method')) then
      call this%params%get ('amg-coarsen-type', ipar)
      INSIST(ipar >= 0) !TODO: replace with proper error handling
      call fHYPRE_ParCSRHybridSetCoarsenType (this%solver, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! Number of AMG smoother sweeps (optional, default 1)
    call this%params%get ('amg-smoothing-sweeps', ipar, default=1)
    INSIST(ipar >= 0) !TODO: replace with proper error handling
    call fHYPRE_ParCSRHybridSetNumSweeps (this%solver, ipar, ierr)
    INSIST(ierr == 0)

    !! The AMG smoothing method (optional, use Hypre default)
    if (this%params%is_parameter('amg-smoothing-method')) then
      call this%params%get ('amg-smoothing-method', ipar)
      INSIST(ipar >= 0) !TODO: replace with proper error handling
      call fHYPRE_ParCSRHybridSetRelaxType (this%solver, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! The AMG interpolation method (optional, use Hypre default)
    if (this%params%is_parameter('amg-interp-method')) then
      call this%params%get ('amg-interp-method', ipar)
      INSIST(ipar >= 0) !TODO: replace with proper error handling
      call fHYPRE_ParCSRHybridSetInterpType (this%solver, ipar, ierr)
      INSIST(ierr == 0)
    end if

    !! Initialize the Hypre matrix.
    ASSERT(associated(this%A))
    call copy_to_ijmatrix (this%A, this%Ah)

    !! Initialize the Hypre solution and RHS vectors.  These are passed to the
    !! solver setup but are apparently ignored for the hybrid solver.
    !! TODO: can we omit this initialization?
    rows = [ (i, i = this%ilower, this%iupper) ]  ! global row indices for this process.
    dummy = 0.0_r8
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.
    call fHYPRE_ParCSRHybridSetup (this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0) !TODO: replace with proper error handling

  end subroutine setup

  subroutine solve (this, b, x, stat)

    class(hypre_hybrid), intent(inout) :: this
    real(r8), intent(in) :: b(:)
    real(r8), intent(inout) :: x(:)
    integer, intent(out) :: stat

    integer :: i, ierr, rows(this%nrows)
    real(r8) :: norm

    call fHYPRE_ClearAllErrors

    !! Global row indices for this process.
    rows = [ (i, i = this%ilower, this%iupper) ]

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, rows, b, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! Solve the system.
    call fHYPRE_ParCSRHybridSolve (this%solver, this%Ah, this%bh, this%xh, ierr)
    if (ierr /= 0) then
      INSIST(ierr == HYPRE_ERROR_CONV)
      stat = 1
      call fHYPRE_ClearAllErrors
    else
      stat = 0
    end if

    !! Retrieve the solution vector from HYPRE.
    call fHYPRE_IJVectorGetValues (this%xh, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)

  end subroutine solve

 !!
 !! This auxillary routine copies a PCSR_MATRIX object SRC to an equivalent
 !! HYPRE_IJMatrix object.  The HYPRE matrix is created if it does not exist.
 !! Otherwise the elements of the existing HYPRE matrix are overwritten with
 !! the values from SRC.  In the latter case the sparsity pattern of the two
 !! matrices must be identical.
 !!
 !! The parallel CSR matrix is defined over local indices, both on-process
 !! and off-process.  By nature of our particular construction, the on-process
 !! rows of this matrix are complete and describe a partitioning of the global
 !! matrix by rows.  The off-process rows, however, are partial and extraneous
 !! and should be ignored.
 !!

  subroutine copy_to_ijmatrix (src, matrix)

    type(pcsr_matrix), intent(in) :: src
    type(hypre_obj), intent(inout) :: matrix

    integer :: j, ierr, ilower, iupper, nrows, nnz
    integer, allocatable :: ncols_onP(:), ncols_offP(:), ncols(:), rows(:), cols(:)

    nrows  = src%graph%row_ip%onP_size()
    ilower = src%graph%row_ip%first_index()
    iupper = src%graph%row_ip%last_index()

    call fHYPRE_ClearAllErrors

    if (.not.hypre_associated(matrix)) then
      call fHYPRE_IJMatrixCreate (ilower, iupper, ilower, iupper, matrix, ierr)
      !! For each row we know how many column entries are on-process and how many
      !! are off-process.  HYPRE is allegedly much faster at forming its CSR matrix
      !! if it knows this info up front.
      allocate(ncols_onP(nrows), ncols_offP(nrows))
      do j = 1, nrows
        ncols_offP(j) = count(src%graph%adjncy(src%graph%xadj(j):src%graph%xadj(j+1)-1) > nrows)
        ncols_onP(j)  = src%graph%xadj(j+1) - src%graph%xadj(j) - ncols_offP(j)
      end do
      call fHYPRE_IJMatrixSetDiagOffdSizes (matrix, ncols_onP, ncols_offP, ierr)
      deallocate(ncols_onP, ncols_offP)
      !! Let HYPRE know that we won't be setting any off-process matrix values.
      call fHYPRE_IJMatrixSetMaxOffProcElmts (matrix, 0, ierr)
      INSIST(ierr == 0)
    end if

    !! After initialization the HYPRE matrix elements can be set.
    call fHYPRE_IJMatrixInitialize (matrix, ierr)
    INSIST(ierr == 0)

    !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    !! nonzero structure of the matrix and the values of those elements. HYPRE
    !! expects global row and column indices.
    nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    allocate(ncols(nrows), rows(nrows), cols(nnz))
    rows = [ (j, j = ilower, iupper) ]
    ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    cols = src%graph%row_ip%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    call fHYPRE_IJMatrixSetValues (matrix, nrows, ncols, rows, cols, src%values, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble (matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module hypre_hybrid_type
