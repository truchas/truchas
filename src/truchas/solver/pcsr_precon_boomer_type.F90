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

#include "f90_assert.fpp"

module pcsr_precon_boomer_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fhypre
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
    !! BoomerAMG parameters
    integer  :: coarsen_type, interp_type, relax_down_type, relax_up_type
    integer  :: num_sweeps = 1    ! number of smoother sweeps
    integer  :: max_levels = 25   ! max number of multigrid levels
    real(r8) :: tol = 0.0d0       ! no tolerance -- using as a preconditioner
    real(r8) :: strong_threshold
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    final :: pcsr_precon_boomer_delete
  end type pcsr_precon_boomer

contains

  !! Final subroutine for PCSR_PRECON_BOOMER objects.
  subroutine pcsr_precon_boomer_delete(this)
    type(pcsr_precon_boomer), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy(this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy(this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy(this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_BoomerAMGDestroy(this%solver, ierr)
    INSIST(ierr == 0)
  end subroutine


  subroutine init(this, A, params, stat, errmsg)

    class(pcsr_precon_boomer), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: ierr
    character(:), allocatable :: context

    this%A => A

    this%nrows  = A%graph%row_imap%onp_size
    this%ilower = A%graph%row_imap%first_gid
    this%iupper = A%graph%row_imap%last_gid

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%bh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%xh, 0, ierr)
    INSIST(ierr == 0)

    !! Process the parameters.
    context = 'processing ' // params%path() // ': '
    call params%get('num-cycles', this%max_iter, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%max_iter <= 0) then
      stat = 1
      errmsg = context // '"num-cycles" must be > 0'
      return
    end if
    call params%get('print-level', this%print_level, stat, errmsg, default=0)
    ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%print_level < 0 .or. this%print_level > 3) then
      stat = 1
      errmsg = context // '"print-level" must be >= 0 and <= 3'
      return
    end if
    call params%get('debug-level', this%debug_level, stat, errmsg, default=0)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%debug_level < 0) then
      stat = 1
      errmsg = context // '"debug-level" must be >= 0'
      return
    end if
    call params%get('logging-level', this%logging_level, stat, errmsg, default=0)
    ! OFF=0, ON=1, >1=residual available from hypre
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%logging_level < 0) then
      stat = 1
      errmsg = context // '"logging-level" must be >= 0'
      return
    end if

    !! Hypre's default coarsening technique is type 10 (v2.23)
    call params%get('coarsen-type', this%coarsen_type, stat, errmsg, default=10)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    !! Default is 0.5, Hypre's recommendation for 3D (heat transfer?)
    call params%get('strong-threshold', this%strong_threshold, stat, errmsg, default=0.5_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    !! Hypre's default interpolation technique is type 6 (v2.23)
    call params%get('interp-type', this%interp_type, stat, errmsg, default=6)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    !! Hypre's default relaxation technique on down cycles is type 13 (v2.23)
    call params%get('relax-down-type', this%relax_down_type, stat, errmsg, default=13)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    !! Hypre's default relaxation technique on up cycles is type 14 (v2.23)
    call params%get('relax-up-type', this%relax_up_type, stat, errmsg, default=14)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

  end subroutine init

  subroutine compute(this)

    class(pcsr_precon_boomer), intent(inout) :: this

    integer :: ierr

    call start_timer('hypre-matrix-copy')
    call copy_to_ijmatrix(this%A, this%Ah)
    call stop_timer('hypre-matrix-copy')

    !! Create the Hypre solver object.  Note that once the solver has
    !! been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it from scratch.
    call start_timer('boomer-setup')
    if (hypre_associated(this%solver)) call fHYPRE_BoomerAMGDestroy(this%solver, ierr)
    call fHYPRE_BoomerAMGCreate(this%solver, ierr)
    INSIST(ierr == 0)

    !! Set some debugging/diagnostic output options.
    call fHYPRE_BoomerAMGSetDebugFlag(this%solver, this%debug_level, ierr)
    call fHYPRE_BoomerAMGSetLogging(this%solver, this%logging_level, ierr)
    INSIST(ierr == 0)

    !! Set the Boomer AMG parameters.
    call fHYPRE_BoomerAMGSetPrintLevel  (this%solver, this%print_level, ierr)
    call fHYPRE_BoomerAMGSetCoarsenType (this%solver, this%coarsen_type, ierr)
    call fHYPRE_BoomerAMGSetInterpType  (this%solver, this%interp_type, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%solver, this%relax_down_type, 1, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%solver, this%relax_up_type, 2, ierr)
    call fHYPRE_BoomerAMGSetNumSweeps   (this%solver, this%num_sweeps, ierr)
    call fHYPRE_BoomerAMGSetMaxLevels   (this%solver, this%max_levels, ierr)
    call fHYPRE_BoomerAMGSetMaxIter     (this%solver, this%max_iter, ierr)
    call fHYPRE_BoomerAMGSetTol         (this%solver, this%tol, ierr)
    call fHYPRE_BoomerAMGSetStrongThreshold (this%solver, this%strong_threshold, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.  Note that B and X are ignored here.
    call fHYPRE_BoomerAMGSetup(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call stop_timer('boomer-setup')

  end subroutine compute

  subroutine apply(this, x)

    class(pcsr_precon_boomer), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: i, ierr, rows(this%nrows)

    ASSERT(size(x) >= this%nrows)

    call start_timer('boomer-solve')

    call fHYPRE_ClearAllErrors

    !! Global row indices for this process.
    rows = [ (i, i = this%ilower, this%iupper) ]

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    x(:this%nrows) = 0.0_r8
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! Call the BoomerAMG solver.
    call fHYPRE_BoomerAMGSolve(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

    !! Retrieve the solution vector from HYPRE
    call fHYPRE_IJVectorGetValues(this%xh, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)

    call stop_timer('boomer-solve')

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

  subroutine copy_to_ijmatrix(src, matrix)

    type(pcsr_matrix), intent(in) :: src
    type(hypre_obj), intent(inout) :: matrix

    integer :: j, ierr, ilower, iupper, nrows, nnz
    integer, allocatable :: ncols_onP(:), ncols_offP(:), ncols(:), rows(:), cols(:)

    nrows  = src%graph%row_imap%onp_size
    ilower = src%graph%row_imap%first_gid
    iupper = src%graph%row_imap%last_gid

    call fHYPRE_ClearAllErrors

    if (.not.hypre_associated(matrix)) then
      call fHYPRE_IJMatrixCreate(ilower, iupper, ilower, iupper, matrix, ierr)
      !! For each row we know how many column entries are on-process and how many
      !! are off-process.  HYPRE is allegedly much faster at forming its CSR matrix
      !! if it knows this info up front.
      allocate(ncols_onP(nrows), ncols_offP(nrows))
      do j = 1, nrows
        ncols_offP(j) = count(src%graph%adjncy(src%graph%xadj(j):src%graph%xadj(j+1)-1) > nrows)
        ncols_onP(j)  = src%graph%xadj(j+1) - src%graph%xadj(j) - ncols_offP(j)
      end do
      call fHYPRE_IJMatrixSetDiagOffdSizes(matrix, ncols_onP, ncols_offP, ierr)
      deallocate(ncols_onP, ncols_offP)
      !! Let HYPRE know that we won't be setting any off-process matrix values.
      call fHYPRE_IJMatrixSetMaxOffProcElmts(matrix, 0, ierr)
      INSIST(ierr == 0)
    end if

    !! After initialization the HYPRE matrix elements can be set.
    call fHYPRE_IJMatrixInitialize(matrix, ierr)
    INSIST(ierr == 0)

    !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    !! nonzero structure of the matrix and the values of those elements. HYPRE
    !! expects global row and column indices.
    nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    allocate(ncols(nrows), rows(nrows), cols(nnz))
    rows = [ (j, j = ilower, iupper) ]
    ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    cols = src%graph%row_imap%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    call fHYPRE_IJMatrixSetValues(matrix, nrows, ncols, rows, cols, src%values, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble(matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module pcsr_precon_boomer_type
