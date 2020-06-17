!!
!! PCSR_PRECON_AMGX_TYPE
!!
!! A concrete implementation of the abstract base class PCSR_PRECON that
!! uses AmgX's preconditioning.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type PCSR_PRECON_AMGX which is an
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

module pcsr_precon_amgx_type

  use kinds, only: r8
  use famgx
  use index_partitioning
  use pcsr_matrix_type
  use pcsr_precon_class
  use parameter_list_type
  use truchas_timers
  implicit none
  private

  type, extends(pcsr_precon), public :: pcsr_precon_amgx
    private
    integer :: nrows = 0, ilower = 0, iupper = 0
    type(amgx_obj) :: solver = amgx_null_obj    ! HYPRE_Solver object handle
    type(amgx_obj) :: Ah = amgx_null_obj        ! HYPRE_IJMatrix object handle
    type(amgx_obj) :: bh = amgx_null_obj, xh = amgx_null_obj  ! HYPRE_IJVector handles
    type(amgx_obj) :: amgx_config = amgx_null_obj, amgx_resources = amgx_null_obj
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
    final :: pcsr_precon_amgx_delete
  end type pcsr_precon_amgx

contains

  !! Final subroutine for PCSR_PRECON_AMGX objects.
  subroutine pcsr_precon_amgx_delete (this)
    type(pcsr_precon_amgx), intent(inout) :: this
    integer :: ierr
    ierr = 0
    if (amgx_associated(this%Ah)) call famgx_matrix_destroy(this%Ah, ierr)
    if (amgx_associated(this%bh)) call famgx_vector_destroy(this%bh, ierr)
    if (amgx_associated(this%xh)) call famgx_vector_destroy(this%xh, ierr)
    if (amgx_associated(this%solver)) call famgx_solver_destroy(this%solver, ierr)
    if (amgx_associated(this%amgx_resources)) call famgx_resources_destroy(this%amgx_resources, ierr)
    if (amgx_associated(this%amgx_config)) call famgx_config_destroy(this%amgx_config, ierr)
    call famgx_finalize_plugins(ierr)
    call famgx_finalize(ierr)
    INSIST(ierr == 0)
  end subroutine pcsr_precon_amgx_delete

  subroutine init (this, A, params)

    class(pcsr_precon_amgx), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(parameter_list) :: params

    integer :: ierr

    this%A => A

    this%nrows  = A%graph%row_ip%onP_size()
    this%ilower = A%graph%row_ip%first_index()
    this%iupper = A%graph%row_ip%last_index()

    !! Process the parameters.
    call params%get('num-cycles', this%max_iter)
    INSIST(this%max_iter > 0)
    call params%get('print-level', this%print_level, default=0)
    INSIST(this%print_level >= 0 .and. this%print_level <= 3)
    call params%get('debug-level', this%debug_level, default=0)
    INSIST(this%debug_level >= 0)
    call params%get('logging-level', this%logging_level, default=0)
    INSIST(this%logging_level >= 0)

    call famgx_initialize(ierr)
    call famgx_initialize_plugins(ierr)
    INSIST(ierr == 0)

    call famgx_config_create(this%amgx_config, "", ierr)
    call famgx_resources_create(this%amgx_resources, this%amgx_config, 1, [0], ierr)
    INSIST(ierr)

    call famgx_vector_create(this%bh, this%amgx_resources, ierr)
    INSIST(ierr == 0)

    call famgx_vector_create(this%xh, this%amgx_resources, ierr)
    INSIST(ierr == 0)

  end subroutine init

  subroutine compute(this)

    class(pcsr_precon_amgx), intent(inout) :: this

    integer :: ierr

    call start_timer('precon-compute')
    call copy_to_ijmatrix(this%amgx_resources, this%A, this%Ah)

    !! Create & setup solver.  Note that B and X are ignored here.
    !! TODO: Is it possible to modify AmgX matrices without destroying the solver
    !!       and rebuilding from scratch?
    call start_timer('setup')
    if (amgx_associated(this%solver)) call famgx_solver_destroy(this%solver, ierr)
    call famgx_solver_create(this%solver, this%amgx_resources, this%amgx_config, ierr)
    INSIST(ierr == 0)
    call famgx_solver_setup(this%solver, this%Ah, ierr)
    INSIST(ierr == 0)
    call stop_timer('setup')

    call stop_timer('precon-compute')

  end subroutine compute

  subroutine apply(this, x)

    class(pcsr_precon_amgx), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: ierr, stat !, rows(this%nrows)

    ASSERT(size(x) >= this%nrows)

    call start_timer('precon-apply')

    ! !! Global row indices for this process.
    ! rows = [(i, i = this%ilower, this%iupper)]

    !! Initialize the AmgX RHS & initial guess vector.
    call start_timer('transfer')
    call famgx_pin_memory(x, ierr); INSIST(ierr == 0)
    call famgx_vector_upload(this%bh, this%nrows, 1, x, ierr); INSIST(ierr == 0)
    call famgx_vector_set_zero(this%xh, this%nrows, 1, ierr); INSIST(ierr == 0)
    call stop_timer('transfer')

    !! Call the AmgX solver.
    call start_timer('solve')
    call famgx_solver_solve_with_0_initial_guess(this%solver, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call famgx_solver_get_status(this%solver, stat, ierr)
    INSIST(ierr == 0)
    INSIST(stat == 0)
    call stop_timer('solve')

    !! Retrieve the solution vector from AmgX
    call start_timer('transfer')
    call famgx_vector_download(this%xh, x, ierr)
    INSIST(ierr == 0)
    call stop_timer('transfer')

    call stop_timer('precon-apply')

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
  subroutine copy_to_ijmatrix(amgx_resources, src, matrix)

    type(amgx_obj), intent(in) :: amgx_resources
    type(pcsr_matrix), intent(in) :: src
    type(amgx_obj), intent(inout) :: matrix

    integer :: ierr, nrow_onP, nrow_global, len_onP
    integer, allocatable :: row_offset(:), col_global(:), global_row_rank(:)
    !integer, allocatable :: ncols_onP(:), ncols_offP(:), ncols(:), rows(:), cols(:)

    call start_timer('matrix-copy')

    if (.not.amgx_associated(matrix)) then
      call famgx_matrix_create(matrix, amgx_resources, ierr)
      INSIST(ierr == 0)
    end if

    !! TODO-WARN: Currently hardwired for serial runs.
    !! Indices are converted to C 0-indexing here.
    ! ilower = src%graph%row_ip%first_index()
    ! iupper = src%graph%row_ip%last_index()
    nrow_onP  = src%graph%row_ip%onP_size()
    len_onP = nrow_onP
    nrow_global = nrow_onP
    allocate(global_row_rank(nrow_global), col_global(len_onP), row_offset(nrow_onP))
    global_row_rank = 0

    call start_timer('transfer')
    call famgx_matrix_upload_all_global(matrix, nrow_global, nrow_onP, len_onP, 1, 1, &
        row_offset, col_global, src%values, &
        0, 0, global_row_rank, ierr) ! here set halo depth & num import rings to 1, 1 in parallel
    INSIST(ierr == 0)
    call stop_timer('transfer')

    call stop_timer('matrix-copy')

    ! call fHYPRE_IJMatrixCreate(ilower, iupper, ilower, iupper, matrix, ierr)
    ! !! For each row we know how many column entries are on-process and how many
    ! !! are off-process.  HYPRE is allegedly much faster at forming its CSR matrix
    ! !! if it knows this info up front.
    ! allocate(ncols_onP(nrows), ncols_offP(nrows))
    ! do j = 1, nrows
    !   ncols_offP(j) = count(src%graph%adjncy(src%graph%xadj(j):src%graph%xadj(j+1)-1) > nrows)
    !   ncols_onP(j)  = src%graph%xadj(j+1) - src%graph%xadj(j) - ncols_offP(j)
    ! end do
    ! call start_timer('transfer')
    ! call fHYPRE_IJMatrixSetDiagOffdSizes(matrix, ncols_onP, ncols_offP, ierr)
    ! deallocate(ncols_onP, ncols_offP)
    ! !! Let HYPRE know that we won't be setting any off-process matrix values.
    ! call fHYPRE_IJMatrixSetMaxOffProcElmts(matrix, 0, ierr)
    ! INSIST(ierr == 0)
    ! call stop_timer('transfer')

    ! !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    ! !! nonzero structure of the matrix and the values of those elements. HYPRE
    ! !! expects global row and column indices.
    ! nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    ! allocate(ncols(nrows), rows(nrows), cols(nnz))
    ! rows = [(j, j = ilower, iupper)]
    ! ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    ! cols = src%graph%row_ip%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    ! call start_timer('transfer')
    ! call fHYPRE_IJMatrixSetValues_v2(matrix, nrows, ncols, rows, cols, src%values, ierr)
    ! deallocate(ncols, rows, cols)
    ! INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module pcsr_precon_amgx_type
