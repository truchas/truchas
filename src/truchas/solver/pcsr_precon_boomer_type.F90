!!
!! PCSR_PRECON_BOOMER_TYPE
!!
!! A concrete implementation of the abstract base class PCSR_PRECON that
!! uses Hypre's BoomerAMG preconditioning.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, February 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type PCSR_PRECON_BOOMER which is an
!!  extension of the abstract base class PCSR_PRECON that uses the BoomerAMG
!!  algebraic multigrid preconditioner from LLNL's Hypre library.  See the
!!  base class comments for a description of the common type bound procedures.
!!
!!  The INIT procedure expects to find the following parameters in the
!!  TYPE(PARAMETER_LIST) argument PARAMS.  Parameters with a default value
!!  are optional; the others are required.
!!
!!    'num-cycles'    - The number of multigrid cycles to apply; > 0.
!!    'print-level'   - Verbosity level: 0, silent (default), 1, info from
!!                    - set-up step; 2, info from solve step; 3, info from
!!                      both set-up and solve steps.
!!    'debug-level'   - 0, disable debugging (default); 1 enable debugging.
!!    'logging-level' - 0, off (default); 1, on; >1, residual info available.
!!
!! IMPLEMENTATION NOTES
!!
!!  The INIT procedure includes some commented-out code from the Pececillo
!!  mini-app for error checking the parameter list.  This code, especially
!!  its error messages, make little sense in the current context, but are
!!  included for future reference.  The provided parameter list is created
!!  by the client code using namelist data read from the input file, and it
!!  is expected that the parameter list is fully vetted at that point.  In
!!  the future I expect the parameter list to be read directly from the input
!!  file, and when that happens this included code may become apropos.
!!  The params derived type is intended as a temporary stand-in for the use of
!!  a generic parameter list capability similiar to the Teuchos::ParameterList
!!  class from Trilinos.  A partial F2003 implementation exists.
!!

#include "f90_assert.fpp"

module pcsr_precon_boomer_type

  use kinds, only: r8
  use fhypre
  use index_partitioning
  use pcsr_matrix_type
  use pcsr_precon_class
  use parameter_list_type
  use truchas_timers
  implicit none
  private

  type, extends(pcsr_precon), public :: pcsr_precon_boomer
    private
    integer :: nrows = 0, ilower = 0, iupper = 0
    type(hypre_obj) :: solver = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: Ah = hypre_null_obj        ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: bh = hypre_null_obj, xh = hypre_null_obj  ! HYPRE_IJVector handles
    !! BoomerAMG parameters -- these are set at initialization
    integer  :: max_iter          ! number of cycles -- using as a preconditioner
    integer  :: print_level       ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: debug_level       ! OFF=0, ON=1
    integer  :: logging_level     ! OFF=0, ON=1, >1=residual available from hypre
    !! BoomerAMG parameters -- these are currently hardwired
    integer  :: coarsen_type = 10 ! HMIS coarsening
    integer  :: relax_type = 18   ! l-scaled Jacobi smoothing
    integer  :: num_sweeps = 1    ! number of smoother sweeps
    integer  :: max_levels = 25   ! max number of multigrid levels
    integer :: keep_transpose = 1 ! use matvecs and avoid transpose matvecs
    real(r8) :: tol = 0.0d0       ! no tolerance -- using as a preconditioner
    real(r8) :: strong_threshold = 0.5_r8 ! should be 0.5 for 3D problems and 0.25 for 2D
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    final :: pcsr_precon_boomer_delete
  end type pcsr_precon_boomer

contains

  !! Final subroutine for PCSR_PRECON_BOOMER objects.
  subroutine pcsr_precon_boomer_delete (this)
    type(pcsr_precon_boomer), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy (this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy (this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy (this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_BoomerAMGDestroy (this%solver, ierr)
    INSIST(ierr == 0)
  end subroutine pcsr_precon_boomer_delete

  subroutine init (this, A, params)

    class(pcsr_precon_boomer), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(parameter_list) :: params

    integer :: ierr

    this%A => A

    this%nrows  = A%graph%row_ip%onP_size()
    this%ilower = A%graph%row_ip%first_index()
    this%iupper = A%graph%row_ip%last_index()

    call fHYPRE_ClearAllErrors
    call fHYPRE_EnableGPU

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%bh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%xh, 0, ierr)
    INSIST(ierr == 0)

    !! Process the parameters.
    call params%get('num-cycles', this%max_iter)
    INSIST(this%max_iter > 0)
    call params%get('print-level', this%print_level, default=0)
    INSIST(this%print_level >= 0 .and. this%print_level <= 3)
    call params%get('debug-level', this%debug_level, default=0)
    INSIST(this%debug_level >= 0)
    call params%get('logging-level', this%logging_level, default=0)
    INSIST(this%logging_level >= 0)

!NNC    !! Process the parameters.
!NNC    use truchas_logging_services
!NNC    integer :: stat
!NNC    character(:), allocatable :: context, errmsg
!NNC    context = 'processing ' // params%name() // ': '
!NNC    call params%get('num-cycles', this%max_iter, stat=stat, errmsg=errmsg)
!NNC    if (stat /= 0) call TLS_fatal (context//errmsg)
!NNC    if (this%max_iter <= 0) call TLS_fatal (context//'"num-cycles" must be > 0')
!NNC    call params%get('print-level', this%print_level, default=0, stat=stat, errmsg=errmsg)
!NNC    if (stat /= 0) call TLS_fatal (context//errmsg)
!NNC    if (this%print_level < 0 .or. this%print_level > 3) call TLS_fatal (context//'"print-level" must be >= 0 and <= 3')
!NNC    call params%get('debug-level', this%debug_level, default=0, stat=stat, errmsg=errmsg)
!NNC    if (stat /= 0) call TLS_fatal (context//errmsg)
!NNC    if (this%debug_level < 0) call TLS_fatal (context//'"debug-level" must be >= 0')
!NNC    call params%get('logging-level', this%logging_level, default=0, stat=stat, errmsg=errmsg)
!NNC    if (stat /= 0) call TLS_fatal (context//errmsg)
!NNC    if (this%logging_level < 0) call TLS_fatal (context//'"logging-level" must be >= 0')

  end subroutine init

  subroutine compute (this)

    class(pcsr_precon_boomer), intent(inout) :: this

    integer :: ierr

    call start_timer ('hypre-matrix-copy')
    call copy_to_ijmatrix (this%A, this%Ah)
    call stop_timer ('hypre-matrix-copy')

    !! Create the Hypre solver object.  Note that once the solver has
    !! been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it from scratch.
    call start_timer ('boomer-setup')
    if (hypre_associated(this%solver)) call fHYPRE_BoomerAMGDestroy (this%solver, ierr)
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
    call fHYPRE_BoomerAMGSetStrongThreshold  (this%solver, this%strong_threshold, ierr)
    call fHYPRE_BoomerAMGSetKeepTranspose (this%solver, this%keep_transpose, ierr)
    call fHYPRE_BoomerAMGSetRAP2 (this%solver, 1, ierr)
    call fHYPRE_BoomerAMGSetModuleRAP2 (this%solver, 1, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.  Note that B and X are ignored here.
    call fHYPRE_BoomerAMGSetup (this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call stop_timer ('boomer-setup')

  end subroutine compute

  subroutine apply (this, x)

    class(pcsr_precon_boomer), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: i, ierr, rows(this%nrows)

    ASSERT(size(x) >= this%nrows)

    call start_timer ('boomer-solve')

    call fHYPRE_ClearAllErrors

    !! Global row indices for this process.
    rows = [ (i, i = this%ilower, this%iupper) ]

    ! !! Initialize the Hypre RHS & initial guess vectors.
    ! x(:this%nrows) = 0.0_r8
    ! call start_timer('transfer')
    ! call start_timer('initialize')
    ! call fHYPRE_IJVectorInitialize_v2 (this%bh, HYPRE_MEMORY_DEVICE, ierr)
    ! call fHYPRE_IJVectorInitialize_v2 (this%xh, HYPRE_MEMORY_DEVICE, ierr)
    ! call stop_timer('initialize')
    ! call fHYPRE_IJVectorSetValues_v2 (this%bh, this%nrows, rows, x, ierr)
    ! call fHYPRE_IJVectorSetValues_v2 (this%xh, this%nrows, rows, x, ierr)
    ! call start_timer('assemble')
    ! call fHYPRE_IJVectorAssemble (this%bh, ierr)
    ! call fHYPRE_IJVectorAssemble (this%xh, ierr)
    ! call stop_timer('assemble')
    ! call stop_timer('transfer')
    ! INSIST(ierr == 0)

    !! Initialize the Hypre RHS vector.
    call start_timer('transfer')
    call start_timer('initialize')
    call fHYPRE_IJVectorInitialize_v2 (this%bh, HYPRE_MEMORY_DEVICE, ierr)
    call stop_timer('initialize')
    call fHYPRE_IJVectorSetValues_v2 (this%bh, this%nrows, rows, x, ierr)
    call start_timer('assemble')
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    call stop_timer('assemble')
    call stop_timer('transfer')
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    x(:this%nrows) = 0.0_r8
    call start_timer('transfer')
    call start_timer('initialize')
    call fHYPRE_IJVectorInitialize_v2 (this%xh, HYPRE_MEMORY_DEVICE, ierr)
    call stop_timer('initialize')
    call fHYPRE_IJVectorSetValues_v2 (this%xh, this%nrows, rows, x, ierr)
    call start_timer('assemble')
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    call stop_timer('assemble')
    call stop_timer('transfer')
    INSIST(ierr == 0)

    !! Call the BoomerAMG solver.
    call start_timer('solve')
    call fHYPRE_BoomerAMGSolve (this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call stop_timer('solve')

    !! Retrieve the solution vector from HYPRE
    call start_timer('transfer')
    call fHYPRE_IJVectorGetValues_v2 (this%xh, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)
    call stop_timer('transfer')

    call stop_timer ('boomer-solve')

  end subroutine apply

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
      call start_timer('transfer')
      call fHYPRE_IJMatrixSetDiagOffdSizes (matrix, ncols_onP, ncols_offP, ierr)
      deallocate(ncols_onP, ncols_offP)
      !! Let HYPRE know that we won't be setting any off-process matrix values.
      call fHYPRE_IJMatrixSetMaxOffProcElmts (matrix, 0, ierr)
      INSIST(ierr == 0)
      call stop_timer('transfer')
    end if

    !! After initialization the HYPRE matrix elements can be set.
    call start_timer('transfer')
    call start_timer('initialize')
    call fHYPRE_IJMatrixInitialize_v2 (matrix, HYPRE_MEMORY_DEVICE, ierr)
    INSIST(ierr == 0)
    call stop_timer('initialize')
    call stop_timer('transfer')

    !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    !! nonzero structure of the matrix and the values of those elements. HYPRE
    !! expects global row and column indices.
    nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    allocate(ncols(nrows), rows(nrows), cols(nnz))
    rows = (/ (j, j = ilower, iupper) /)
    ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    cols = src%graph%row_ip%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    call start_timer('transfer')
    call fHYPRE_IJMatrixSetValues_v2 (matrix, nrows, ncols, rows, cols, src%values, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call start_timer('assemble')
    call fHYPRE_IJMatrixAssemble (matrix, ierr)
    INSIST(ierr == 0)
    call stop_timer('assemble')
    call stop_timer('transfer')

  end subroutine copy_to_ijmatrix

end module pcsr_precon_boomer_type
