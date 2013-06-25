!!
!! BOOMER_AMG_PRECON_TYPE
!!
!! This module defines a derived type and associated procedures that describe
!! a Hypre BoomerAMG preconditioner for a parallel CSR matrix.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived data type BOOMER_AMG_PRECON (private
!!  components) and the following procedures that operate on instances of the
!!  type passed as the THIS argument.  An instance describes a BoomerAMG
!!  preconditioner for a parallel CSR matrix of type PCSR_MATRIX.
!!
!!  CALL BAMG_PRECON_INIT (THIS, A, PARAMS) initializes the object to be a
!!    preconditioner for the parallel CSR matrix A of type PCSR_MATRIX.  The
!!    object holds a reference to the matrix A, and so the matrix must never
!!    go out of scope during the lifetime of the object.  Morover the actual
!!    argument must be a pointer or have the target attribute.  Only the
!!    structure of the matrix A needs to be defined at this point; the matrix
!!    values are not referenced.  Note that the matrix will not be modified
!!    in any way by this, or the other procedures.
!!
!!    PARAMS is an intent-in argument of type BOOMER_AMG_PRECON_PARAMS which
!!    has the following components:
!!      MAX_ITER        The number of multigrid cycles to apply; default 1.
!!      PRINT_LEVEL     Verbosity level: 0, silent (default), 1, info from
!!                      set-up step; 2, info from solve step; 3, info from
!!                      both set-up and solve steps.
!!      DEBUG_LEVEL     0, disable debugging (default); 1 enable debugging.
!!      LOGGING_LEVEL   0, off (default); 1, on; >1, residual info available.
!!
!!  CALL BAMG_PRECON_COMPUTE (THIS) performs the final setup and configuration
!!    of the preconditioner.  It is at this point the values of the matrix A
!!    are referenced.  It must be called before using the APPLY procedure and
!!    after the matrix values are defined, and must be called again whenever
!!    the matrix values are modified.
!!
!!  CALL BAMG_PRECON_APPLY (THIS, X) applies the preconditioner to the vector X.
!!    The size of X must be at least A%NROW_ONP; only the the initial A%NROW_ONP
!!    elements will be referenced or modified.
!!
!!  BAMG_PRECON_MATRIX(THIS) returns a pointer to the parallel CSR matrix A
!!    with which the preconditioner was initialized.
!!
!!  CALL BAMG_PRECON_DELETE (THIS) frees resources allocated by the object.
!!
!! IMPLEMENTATION NOTES
!!
!!  The params derived type is intended as a temporary stand-in for the use of
!!  a generic parameter list capability similiar to the Teuchos::ParameterList
!!  class from Trilinos.  A partial F2003 implementation exists.
!!
!!  Transitioning to F2003. The design is intended make the change to an OO
!!  implementation straightforward:
!!  * All procedures become type-bound; delete is a final procedure.
!!  * Drop the BAMG_PRECON_ prefix from the method names and don't export
!!    them as loose procedures.
!!  * The interface is identical in form to that for SSOR_PRECON.
!!    Both should extend the same abstract base type, which implements
!!    the MATRIX function and defines the deferred subroutines COMPUTE
!!    and APPLY.  It could also define a deferred subroutine INIT if it
!!    were modified to accept a generic 'parameter list' argument as
!!    described above.  Otherwise the INIT method would need to be
!!    specific to a particular implementation of the base type.
!!

#include "f90_assert.fpp"

module boomer_amg_precon_type

  use kinds
  use index_partitioning
  use parallel_csr_matrix
  use fhypre
  implicit none
  private

  type, public :: boomer_amg_precon
    private
    type(pcsr_matrix), pointer :: Asrc => null()
    integer :: nrows = 0, ilower = 0, iupper = 0
    integer(hypre_obj) :: solver = 0    ! HYPRE_Solver object handle
    integer(hypre_obj) :: A = 0         ! HYPRE_IJMatrix object handle
    integer(hypre_obj) :: b = 0, x = 0  ! HYPRE_IJVector handles
    !! BoomerAMG parameters -- these are set at initialization
    integer  :: max_iter          ! number of cycles -- using as a preconditioner
    integer  :: print_level       ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: debug_level       ! OFF=0, ON=1
    integer  :: logging_level     ! OFF=0, ON=1, >1=residual available from hypre
    !! BoomerAMG parameters -- these are currently hardwired
    integer  :: coarsen_type = 6  ! Falgout coarsening
    integer  :: relax_type = 3    ! hybrid Gauss-Seidel smoothing
    integer  :: num_sweeps = 1    ! number of smoother sweeps
    integer  :: max_levels = 25   ! max number of multigrid levels
    real(r8) :: tol = 0.0d0       ! no tolerance -- using as a preconditioner
    real(r8) :: strong_threshold = 0.5_r8 ! should be 0.5 for 3D problems and 0.25 for 2D
  end type boomer_amg_precon

  public :: bamg_precon_init, bamg_precon_compute, bamg_precon_delete
  public :: bamg_precon_apply, bamg_precon_matrix

  type, public :: boomer_amg_precon_params
    integer  :: max_iter = 1      ! number of cycles -- using as a preconditioner
    integer  :: print_level = 0   ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: debug_level = 0   ! OFF=0, ON=1
    integer  :: logging_level = 0 ! OFF=0, ON=1, >1=residual available from hypre
  end type boomer_amg_precon_params

contains

  subroutine bamg_precon_init (this, A, params)

    type(boomer_amg_precon), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(boomer_amg_precon_params), intent(in) :: params

    integer :: ierr

    this%Asrc => A

    this%nrows  = onP_size(A%graph%row_ip)
    this%ilower = first_index(A%graph%row_ip)
    this%iupper = last_index(A%graph%row_ip)

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%b, ierr)
    call fHYPRE_IJVectorSetMaxOffPValues (this%b, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%x, ierr)
    call fHYPRE_IJVectorSetMaxOffPValues (this%x, 0, ierr)
    INSIST(ierr == 0)

    !! Process the parameters.
    INSIST(params%max_iter > 0)
    this%max_iter = params%max_iter
    INSIST(params%print_level >= 0 .and. params%print_level <= 3)
    this%print_level = params%print_level
    INSIST(params%debug_level >= 0)
    this%debug_level = params%debug_level
    INSIST(params%logging_level >= 0)
    this%logging_level = params%logging_level

  end subroutine bamg_precon_init

  subroutine bamg_precon_delete (this)
    type(boomer_amg_precon), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (this%A /= 0) call fHYPRE_IJMatrixDestroy (this%A, ierr)
    if (this%b /= 0) call fHYPRE_IJVectorDestroy (this%b, ierr)
    if (this%x /= 0) call fHYPRE_IJVectorDestroy (this%x, ierr)
    if (this%solver /= 0) call fHYPRE_BoomerAMGDestroy (this%solver, ierr)
    INSIST(ierr == 0)
  end subroutine bamg_precon_delete

  function bamg_precon_matrix (this) result (matrix)
    type(boomer_amg_precon), intent(in) :: this
    type(pcsr_matrix), pointer :: matrix
    matrix => this%Asrc
  end function bamg_precon_matrix

  subroutine bamg_precon_compute (this)

    type(boomer_amg_precon), intent(inout) :: this

    integer :: ierr

    call copy_to_ijmatrix (this%Asrc, this%A)

    !! Create the Hypre solver object.  Note that once the solver has
    !! been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it from scratch.
    if (this%solver /= 0) then
      call fHYPRE_BoomerAMGDestroy (this%solver, ierr)
      this%solver = 0
    end if
    call fHYPRE_BoomerAMGCreate (this%solver, ierr)
    INSIST(ierr == 0)

    !! Set some debugging/diagnostic output options.
    call fHYPRE_BoomerAMGSetDebugFlag (this%solver, this%debug_level, ierr)
    call fHYPRE_BoomerAMGSetLogging (this%solver, this%logging_level, ierr)
    INSIST(ierr == 0)

    !! Set the Boomer AMG parameters.
    call fHYPRE_BoomerAMGSetPrintLevel  (this%solver, this%print_level, ierr)
    call fHYPRE_BoomerAMGSetCoarsenType (this%solver, this%coarsen_type, ierr)
    call fHYPRE_BoomerAMGSetRelaxType   (this%solver, this%relax_type, ierr)
    call fHYPRE_BoomerAMGSetNumSweeps   (this%solver, this%num_sweeps, ierr)
    call fHYPRE_BoomerAMGSetMaxLevels   (this%solver, this%max_levels, ierr)
    call fHYPRE_BoomerAMGSetMaxIter     (this%solver, this%max_iter, ierr)
    call fHYPRE_BoomerAMGSetTol         (this%solver, this%tol, ierr)
    call fHYPRE_BoomerAMGSetStrongThld  (this%solver, this%strong_threshold, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.  Note that B and X are ignored here.
    call fHYPRE_BoomerAMGSetup (this%solver, this%A, this%b, this%x, ierr)
    INSIST(ierr == 0)

  end subroutine bamg_precon_compute

  subroutine bamg_precon_apply (this, x)

    type(boomer_amg_precon), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: i, ierr, rows(this%nrows)

    ASSERT(size(x) >= this%nrows)

    call fHYPRE_ClearAllErrors

    !! Global row indices for this process.
    rows = (/ (i, i = this%ilower, this%iupper) /)

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize (this%b, ierr)
    call fHYPRE_IJVectorSetValues  (this%b, this%nrows, rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%b, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    x(:this%nrows) = 0.0_r8
    call fHYPRE_IJVectorInitialize (this%x, ierr)
    call fHYPRE_IJVectorSetValues  (this%x, this%nrows, rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%x, ierr)
    INSIST(ierr == 0)

    !! Call the BoomerAMG solver.
    call fHYPRE_BoomerAMGSolve (this%solver, this%A, this%b, this%x, ierr)
    INSIST(ierr == 0)

    !! Retrieve the solution vector from HYPRE
    call fHYPRE_IJVectorGetValues (this%x, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)

  end subroutine bamg_precon_apply

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
    integer(hypre_obj), intent(inout) :: matrix

    integer :: j, ierr, ilower, iupper, nrows, nnz
    integer, allocatable :: ncols_onP(:), ncols_offP(:), ncols(:), rows(:), cols(:)

    nrows  = onP_size(src%graph%row_ip)
    ilower = first_index(src%graph%row_ip)
    iupper = last_index(src%graph%row_ip)

    call fHYPRE_ClearAllErrors

    if (matrix == 0) then
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
      call fHYPRE_IJMatrixSetMaxOffPValues (matrix, 0, ierr)
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
    call fHYPRE_IJMatrixSetValues (matrix, nrows, ncols, rows, cols, src%data, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble (matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module boomer_amg_precon_type
