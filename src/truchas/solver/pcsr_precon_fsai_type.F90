!!
!! PCSR_PRECON_FSAI_TYPE
!!
!! A concrete implementation of the abstract base class PCSR_PRECON that
!! uses Hypre's FSAI preconditioning.
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

module pcsr_precon_fsai_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fhypre
  use pcsr_matrix_type
  use pcsr_precon_class
  use parameter_list_type
  use truchas_timers
  implicit none
  private

  type, extends(pcsr_precon), public :: pcsr_precon_fsai
    private
    integer :: nrows = 0, ilower = 0, iupper = 0
    type(hypre_obj) :: solver = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: Ah = hypre_null_obj        ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: bh = hypre_null_obj, xh = hypre_null_obj  ! HYPRE_IJVector handles
    !! FSAI parameters -- these are set at initialization
    integer  :: max_steps, max_step_size
    real(r8) :: kap_tol
    integer  :: print_level ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    final :: pcsr_precon_fsai_delete
  end type pcsr_precon_fsai

contains

  !! Final subroutine for PCSR_PRECON_fsai objects.
  subroutine pcsr_precon_fsai_delete(this)
    type(pcsr_precon_fsai), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy(this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy(this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy(this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_FSAIDestroy(this%solver, ierr)
    INSIST(ierr == 0)
  end subroutine


  subroutine init(this, A, params, stat, errmsg)

    class(pcsr_precon_fsai), intent(out) :: this
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

    call params%get('max-steps', this%max_steps, stat, errmsg, default=5)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%max_steps <= 0) then
      stat = 1
      errmsg = context // '"max-steps" must be > 0'
      return
    end if

    call params%get('max-step-size', this%max_step_size, stat, errmsg, default=3)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%max_step_size <= 0) then
      stat = 1
      errmsg = context // '"max-step-size" must be > 0'
      return
    end if

    call params%get('kap-tol', this%kap_tol, stat, errmsg, default=1e-3_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%kap_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"kap-tol" must be >= 0.0'
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

  end subroutine init

  subroutine compute(this)

    class(pcsr_precon_fsai), intent(inout) :: this

    integer :: ierr

    call start_timer('hypre-matrix-copy')
    call this%A%copy_to_ijmatrix(this%Ah)
    call stop_timer('hypre-matrix-copy')

    !! Create the Hypre solver object.  Note that once the solver has
    !! been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it from scratch.
    call start_timer('fsai-setup')
    if (hypre_associated(this%solver)) call fHYPRE_FSAIDestroy(this%solver, ierr)
    call fHYPRE_FSAICreate(this%solver, ierr)
    INSIST(ierr == 0)

    !! Set some debugging/diagnostic output options.
    call fHYPRE_FSAISetPrintLevel  (this%solver, this%print_level, ierr)
    INSIST(ierr == 0)

    !! Set the FSAI parameters for use as a preconditioner.
    call fHYPRE_FSAISetTolerance(this%solver, 0.0_r8, ierr)
    call fHYPRE_FSAISetMaxIterations(this%solver, 1, ierr)
    call fHYPRE_FSAISetZeroGuess(this%solver, 1, ierr)
    INSIST(ierr == 0)

    !! Set the FSAI parameters.
    call fHYPRE_FSAISetAlgoType(this%solver, 1, ierr)
    call fHYPRE_FSAISetMaxSteps(this%solver, this%max_steps, ierr)
    call fHYPRE_FSAISetMaxStepSize(this%solver, this%max_step_size, ierr)
    call fHYPRE_FSAISetKapTolerance(this%solver, this%kap_tol, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.  Note that B and X are ignored here.
    call fHYPRE_FSAISetup(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call stop_timer('fsai-setup')

  end subroutine compute

  subroutine apply(this, x)

    class(pcsr_precon_fsai), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: i, ierr, rows(this%nrows)

    ASSERT(size(x) >= this%nrows)

    call start_timer('fsai-solve')

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

    !! Call the FSAI solver.
    call fHYPRE_FSAISolve(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

    !! Retrieve the solution vector from HYPRE
    call fHYPRE_IJVectorGetValues(this%xh, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)

    call stop_timer('fsai-solve')

  end subroutine apply

end module pcsr_precon_fsai_type
