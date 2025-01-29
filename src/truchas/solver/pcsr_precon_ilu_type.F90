!!
!! PCSR_PRECON_ILU_TYPE
!!
!! A concrete implementation of the abstract base class PCSR_PRECON that
!! uses Hypre's ILU preconditioning.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!! January 2025
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module pcsr_precon_ilu_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fhypre
  use pcsr_matrix_type
  use pcsr_precon_class
  use parameter_list_type
  use truchas_timers
  implicit none
  private

  type, extends(pcsr_precon), public :: pcsr_precon_ilu
    private
    integer :: nrows = 0, ilower = 0, iupper = 0
    type(hypre_obj) :: solver = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: Ah = hypre_null_obj        ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: bh = hypre_null_obj, xh = hypre_null_obj  ! HYPRE_IJVector handles
    !! ILU parameters -- these are set at initialization
    integer :: fill_level
    integer :: print_level ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer :: logging_level ! NONE=0, RESIDUAL>=1
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    final :: pcsr_precon_ilu_delete
  end type

contains

  !! Final subroutine for PCSR_PRECON_ILU objects.
  subroutine pcsr_precon_ilu_delete(this)
    type(pcsr_precon_ilu), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy(this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy(this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy(this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_ILUDestroy(this%solver, ierr)
    INSIST(ierr == 0)
  end subroutine


  subroutine init(this, A, params, stat, errmsg)

    class(pcsr_precon_ilu), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: ierr, comm
    character(:), allocatable :: context

    this%A => A

    comm = A%graph%row_imap%comm
    this%nrows  = A%graph%row_imap%onp_size
    this%ilower = A%graph%row_imap%first_gid
    this%iupper = A%graph%row_imap%last_gid

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%bh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%xh, 0, ierr)
    INSIST(ierr == 0)

    !! Process the parameters.
    context = 'processing ' // params%path() // ': '

    !! Level of fill for ILU(k); default is 0 for ILU(0)
    call params%get('fill-level', this%fill_level, stat, errmsg, default=0)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%fill_level < 0) then
      stat = 1
      errmsg = context // '"fill-level" must be >= 0'
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

  end subroutine init

  subroutine compute(this)

    class(pcsr_precon_ilu), intent(inout) :: this

    integer :: ierr

    call start_timer('hypre-matrix-copy')
    call this%A%copy_to_ijmatrix(this%Ah)
    call stop_timer('hypre-matrix-copy')

    !! Create the Hypre solver object.  Note that once the solver has
    !! been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it from scratch.
    call start_timer('ilu-setup')
    if (hypre_associated(this%solver)) call fHYPRE_ILUDestroy(this%solver, ierr)
    call fHYPRE_ILUCreate(this%solver, ierr)
    INSIST(ierr == 0)

    !! Set some debugging/diagnostic output options.
    call fHYPRE_ILUSetLogging(this%solver, this%logging_level, ierr)
    call fHYPRE_ILUSetPrintLevel(this%solver, this%print_level, ierr)
    INSIST(ierr == 0)

    !! Set the ILU parameters for use as a preconditioner.
    call fHYPRE_ILUSetTol(this%solver, 0.0_r8, ierr)
    call fHYPRE_ILUSetMaxIter(this%solver, 1, ierr)
    INSIST(ierr == 0)

    !! Hardwire for ILU(k) with user-specified k.
    call fHYPRE_ILUSetType(this%solver, 0, ierr) ! default, block Jacobi ILU(k)
    call fHYPRE_ILUSetLevelOfFill(this%solver, this%fill_level, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go. Note that B and X are ignored here.
    call fHYPRE_ILUSetup(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call stop_timer('ilu-setup')

  end subroutine compute

  subroutine apply(this, x)

    class(pcsr_precon_ilu), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: i, ierr, rows(this%nrows)

    ASSERT(size(x) >= this%nrows)

    call start_timer('ilu-solve')

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

    !! Call the ILU solver.
    call fHYPRE_ILUSolve(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

    !! Retrieve the solution vector from HYPRE
    call fHYPRE_IJVectorGetValues(this%xh, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)

    call stop_timer('ilu-solve')

  end subroutine apply

end module pcsr_precon_ilu_type
