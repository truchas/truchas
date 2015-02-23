!!
!! HYPRE_PCG_TYPE
!!
!! This module defines a preconditioned CG solver class built on the
!! ParCSR PCG solver from Hypre that uses BoomerAMG as the preconditioner.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!

#include "f90_assert.fpp"

module hypre_pcg_type

  use kinds
  use pcsr_matrix_type
  use index_partitioning
  use fhypre
  implicit none
  private

  type, public :: hypre_pcg
    private
    type(pcsr_matrix), pointer :: Asrc => null()
    integer :: nrows = 0, ilower = 0, iupper = 0
    integer, pointer :: rows(:) => null()
    type(hypre_obj) :: solver = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: precon = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: A = hypre_null_obj         ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: b = hypre_null_obj, x = hypre_null_obj  ! HYPRE_IJVector object handles
    !! PCG parameters -- these are set at initialization
    integer  :: print_level       ! ??? OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: max_iter          ! maximum number of iterations
    real(r8) :: err_tol           ! error tolerance (||r||_2 / ||b||_2)
    !! BoomerAMG preconditioner parameters -- these are set at initialization
    integer  :: num_cycles        ! number of cycles
    integer  :: debug_level       ! OFF=0, ON=1
    integer  :: logging_level     ! OFF=0, ON=1, >1=residual available from hypre
    !! BoomerAMG preconditioner parameters -- these are currently hardwired
    integer  :: coarsen_type = 6  ! Falgout coarsening
    integer  :: relax_type = 3    ! hybrid Gauss-Seidel smoothing
    integer  :: num_sweeps = 1    ! number of smoother sweeps
    integer  :: max_levels = 25   ! max number of multigrid levels
    real(r8) :: strong_threshold = 0.5_r8 ! should be 0.5 for 3D problems and 0.25 for 2D
  end type hypre_pcg

  public :: hypre_pcg_init, hypre_pcg_compute, hypre_pcg_solve, hypre_pcg_delete
  public :: hypre_pcg_matrix, hypre_pcg_get_metrics

  type, public :: hypre_pcg_params
    real(r8) :: err_tol           ! error tolerance (||r||_2 / ||b||_2)
    integer  :: max_iter = 100    ! maximum number of iterations
    integer  :: num_cycles = 1    ! number of BoomerAMG preconditioner cycles
    integer  :: print_level = 0   ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: debug_level = 0   ! OFF=0, ON=1
    integer  :: logging_level = 0 ! OFF=0, ON=1, >1=residual available from hypre
  end type hypre_pcg_params

contains

  subroutine hypre_pcg_init (this, A, params)

    type(hypre_pcg), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(hypre_pcg_params), intent(in) :: params

    integer :: ierr, i

    this%Asrc => A

    this%nrows  = onP_size(A%graph%row_ip)
    this%ilower = first_index(A%graph%row_ip)
    this%iupper = last_index(A%graph%row_ip)
    allocate(this%rows(this%nrows))
    this%rows = (/ (i, i = this%ilower, this%iupper) /)

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%b, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%b, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%x, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%x, 0, ierr)
    INSIST(ierr == 0)

    !! Process the parameters.
    INSIST(params%err_tol > 0.0_r8)
    this%err_tol = params%err_tol
    INSIST(params%max_iter > 0)
    this%max_iter = params%max_iter
    INSIST(params%num_cycles > 0)
    this%num_cycles = params%num_cycles
    INSIST(params%print_level >= 0 .and. params%print_level <= 3)
    this%print_level = params%print_level
    INSIST(params%debug_level >= 0)
    this%debug_level = params%debug_level
    INSIST(params%logging_level >= 0)
    this%logging_level = params%logging_level

  end subroutine hypre_pcg_init

  subroutine hypre_pcg_delete (this)
    type(hypre_pcg), intent(inout) :: this
    type(hypre_pcg) :: default
    integer :: ierr
    ierr = 0
    if (associated(this%rows)) deallocate(this%rows)
    if (hypre_associated(this%A)) call fHYPRE_IJMatrixDestroy (this%A, ierr)
    if (hypre_associated(this%b)) call fHYPRE_IJVectorDestroy (this%b, ierr)
    if (hypre_associated(this%x)) call fHYPRE_IJVectorDestroy (this%x, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_PCGDestroy (this%solver, ierr)
    if (hypre_associated(this%precon)) call fHYPRE_BoomerAMGDestroy (this%precon, ierr)
    INSIST(ierr == 0)
    this = default  ! assign default initialization values
  end subroutine hypre_pcg_delete

  function hypre_pcg_matrix (this) result (matrix)
    type(hypre_pcg), intent(in) :: this
    type(pcsr_matrix), pointer :: matrix
    matrix => this%Asrc
  end function hypre_pcg_matrix

  subroutine hypre_pcg_get_metrics (this, num_itr)
    type(hypre_pcg), intent(in) :: this
    integer, intent(out), optional :: num_itr
    integer :: ierr
    INSIST(hypre_associated(this%solver))
    if (present(num_itr)) then
      call fHYPRE_PCGGetNumIterations (this%solver, num_itr, ierr)
      INSIST(ierr == 0)
    end if
  end subroutine hypre_pcg_get_metrics

  subroutine hypre_pcg_compute (this)

    type(hypre_pcg), intent(inout) :: this

    integer :: ierr
    real(r8) :: dummy(this%nrows)

    ASSERT(associated(this%Asrc))
    call copy_to_ijmatrix (this%Asrc, this%A)

    !! Create the Hypre BoomerAMG solver object (preconditioner).  Note that
    !! once the solver has been setup, it is not possible to change the matrix
    !! values without completely destroying the solver and recreating it.
    if (hypre_associated(this%precon)) call fHYPRE_BoomerAMGDestroy (this%precon, ierr)
    call fHYPRE_BoomerAMGCreate (this%precon, ierr)
    INSIST(ierr == 0)

    !! Set some debugging/diagnostic output options.
    call fHYPRE_BoomerAMGSetDebugFlag (this%precon, this%debug_level, ierr)
    call fHYPRE_BoomerAMGSetLogging (this%precon, this%logging_level, ierr)
    INSIST(ierr == 0)

    !! Set the Boomer AMG parameters.
    call fHYPRE_BoomerAMGSetPrintLevel  (this%precon, this%print_level, ierr)
    call fHYPRE_BoomerAMGSetCoarsenType (this%precon, this%coarsen_type, ierr)
    call fHYPRE_BoomerAMGSetRelaxType   (this%precon, this%relax_type, ierr)
    call fHYPRE_BoomerAMGSetNumSweeps   (this%precon, this%num_sweeps, ierr)
    call fHYPRE_BoomerAMGSetMaxLevels   (this%precon, this%max_levels, ierr)
    call fHYPRE_BoomerAMGSetMaxIter     (this%precon, this%num_cycles, ierr)
    call fHYPRE_BoomerAMGSetTol         (this%precon, 0.0_r8, ierr)
    call fHYPRE_BoomerAMGSetStrongThreshold  (this%precon, this%strong_threshold, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 3, 1, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 4, 2, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 9, 3, ierr)
    INSIST(ierr == 0)

    !! Create the Hypre PCG solver object. This supposes that once the solver
    !! has been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it.  VERIFY THIS.
    if (hypre_associated(this%solver)) call fHYPRE_PCGDestroy (this%solver, ierr)
    call fHYPRE_PCGCreate (this%solver, ierr)
    INSIST(ierr == 0)

    !! Set the PCG parameters.
    call fHYPRE_PCGSetPrintLevel (this%solver, this%print_level, ierr)
    call fHYPRE_PCGSetTwoNorm (this%solver, 1, ierr)
    call fHYPRE_PCGSetTol (this%solver, this%err_tol, ierr)
    call fHYPRE_PCGSetAbsoluteTol (this%solver, 0.0_r8, ierr)
    call fHYPRE_PCGSetMaxIter (this%solver, this%max_iter, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre solution and RHS vectors;
    !! structure (not values) needed to setup PCG solver.
    dummy = 0.0_r8
    call fHYPRE_IJVectorInitialize (this%b, ierr)
    call fHYPRE_IJVectorSetValues  (this%b, this%nrows, this%rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%b, ierr)
    INSIST(ierr == 0)
    call fHYPRE_IJVectorInitialize (this%x, ierr)
    call fHYPRE_IJVectorSetValues  (this%x, this%nrows, this%rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%x, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.
    call fHYPRE_PCGSetPrecond (this%solver, this%precon, ierr)
    call fHYPRE_PCGSetup (this%solver, this%A, this%b, this%x, ierr)
    INSIST(ierr == 0)

  end subroutine hypre_pcg_compute

  subroutine hypre_pcg_solve (this, b, x)

    type(hypre_pcg), intent(inout) :: this
    real(r8), intent(in) :: b(:)
    real(r8), intent(inout) :: x(:)

    integer :: ierr

    call fHYPRE_ClearAllErrors

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize (this%b, ierr)
    call fHYPRE_IJVectorSetValues  (this%b, this%nrows, this%rows, b, ierr)
    call fHYPRE_IJVectorAssemble   (this%b, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    call fHYPRE_IJVectorInitialize (this%x, ierr)
    call fHYPRE_IJVectorSetValues  (this%x, this%nrows, this%rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%x, ierr)
    INSIST(ierr == 0)

    !! Solve the system.
    call fHYPRE_PCGSolve (this%solver, this%A, this%b, this%x, ierr)
    INSIST(ierr == 0)

    !! Retrieve the solution vector from HYPRE.
    call fHYPRE_IJVectorGetValues (this%x, this%nrows, this%rows, x, ierr)
    INSIST(ierr == 0)

  end subroutine hypre_pcg_solve

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

    nrows  = onP_size(src%graph%row_ip)
    ilower = first_index(src%graph%row_ip)
    iupper = last_index(src%graph%row_ip)

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
    rows = (/ (j, j = ilower, iupper) /)
    ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    cols = global_index(src%graph%row_ip, src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    call fHYPRE_IJMatrixSetValues (matrix, nrows, ncols, rows, cols, src%values, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble (matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module hypre_pcg_type
