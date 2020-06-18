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
    integer :: nrows, nnz !, ilower = 0, iupper = 0
    logical :: amgx_initialized = .false., matrix_initialized = .false.
    type(amgx_obj) :: solver = amgx_null_obj    ! amgx_solver object handle
    type(amgx_obj) :: Ah = amgx_null_obj        ! amgx_matrix object handle
    type(amgx_obj) :: bh = amgx_null_obj, xh = amgx_null_obj  ! amgx_vector handles
    type(amgx_obj) :: amgx_config = amgx_null_obj, amgx_resources = amgx_null_obj

    !! AmgX parameters -- these are set at initialization
    integer  :: max_iter      ! number of cycles -- using as a preconditioner
    integer  :: print_level   ! OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: debug_level   ! OFF=0, ON=1
    integer  :: logging_level ! OFF=0, ON=1, >1=residual available from hypre
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    final :: pcsr_precon_amgx_delete
  end type pcsr_precon_amgx

contains

  !! Final subroutine for PCSR_PRECON_AMGX objects.
  subroutine pcsr_precon_amgx_delete(this)
    type(pcsr_precon_amgx), intent(inout) :: this
    integer :: ierr
    ierr = 0
    if (amgx_associated(this%Ah)) call famgx_matrix_destroy(this%Ah, ierr)
    if (amgx_associated(this%bh)) call famgx_vector_destroy(this%bh, ierr)
    if (amgx_associated(this%xh)) call famgx_vector_destroy(this%xh, ierr)
    if (amgx_associated(this%solver)) call famgx_solver_destroy(this%solver, ierr)
    if (amgx_associated(this%amgx_resources)) call famgx_resources_destroy(this%amgx_resources, ierr)
    if (amgx_associated(this%amgx_config)) call famgx_config_destroy(this%amgx_config, ierr)
    if (this%amgx_initialized) then
      call famgx_finalize_plugins(ierr)
      call famgx_finalize(ierr)
    end if
    INSIST(ierr == 0)
    this%amgx_initialized = .false. ! TODO: do these need to be here?
    this%matrix_initialized = .false.
  end subroutine pcsr_precon_amgx_delete


  subroutine init(this, A, params)

    class(pcsr_precon_amgx), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(parameter_list) :: params

    integer :: ierr

    this%A => A

    this%nrows  = A%graph%row_ip%onP_size()
    this%nnz = A%graph%xadj(this%nrows+1) - A%graph%xadj(1)

    !! Process the parameters.
    call params%get('num-cycles', this%max_iter)
    INSIST(this%max_iter > 0)
    call params%get('print-level', this%print_level, default=0)
    INSIST(this%print_level >= 0 .and. this%print_level <= 2)
    call params%get('debug-level', this%debug_level, default=0)
    INSIST(this%debug_level >= 0)
    call params%get('logging-level', this%logging_level, default=0)
    INSIST(this%logging_level >= 0)

    print *, "Initializing AmgX ..."
    this%print_level = 2
    !this%logging_level = 1

    call famgx_initialize(ierr)
    INSIST(ierr == 0)
    call famgx_initialize_plugins(ierr)
    INSIST(ierr == 0)
    this%amgx_initialized = .true.

    call famgx_config_create(this%amgx_config, config_string(this), ierr)
    INSIST(ierr == 0)
    call famgx_resources_create(this%amgx_resources, this%amgx_config, 1, [0], ierr)
    INSIST(ierr == 0)
    call famgx_solver_create(this%solver, this%amgx_resources, this%amgx_config, ierr)
    INSIST(ierr == 0)

    call famgx_matrix_create(this%Ah, this%amgx_resources, ierr); INSIST(ierr == 0)
    call famgx_vector_create(this%bh, this%amgx_resources, ierr); INSIST(ierr == 0)
    call famgx_vector_create(this%xh, this%amgx_resources, ierr); INSIST(ierr == 0)
    !! TODO-WARN: famgx_vector_bind

    print *, "AmgX initialized."

  end subroutine init


  !! Copy matrix & setup solver. Note that B and X are ignored here.
  subroutine compute(this)

    class(pcsr_precon_amgx), intent(inout) :: this

    integer :: ierr

    call start_timer('precon-compute')
    call start_timer('memcopy-matrix')
    if (.not.this%matrix_initialized) then
      call copy_to_ijmatrix(this%A, this%Ah)
      this%matrix_initialized = .true.
    else
      call famgx_pin_memory(this%A%values, ierr); INSIST(ierr == 0)
      call famgx_matrix_replace_coefficients(this%Ah, this%nrows, this%nnz, this%A%values, ierr)
      INSIST(ierr == 0)
      call famgx_unpin_memory(this%A%values, ierr); INSIST(ierr == 0)
    end if
    call stop_timer('memcopy-matrix')

    print *, "AmgX: setup solver"
    call start_timer('setup')
    call famgx_solver_setup(this%solver, this%Ah, ierr)
    INSIST(ierr == 0)
    call stop_timer('setup')

    call stop_timer('precon-compute')

  end subroutine compute


  subroutine apply(this, x)

    class(pcsr_precon_amgx), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: ierr, stat

    ASSERT(size(x) >= this%nrows)

    call start_timer('precon-apply')

    !! Initialize the AmgX RHS & initial guess vector.
    call start_timer('memcpy-vector')
    call famgx_pin_memory(x, ierr); INSIST(ierr == 0)
    call famgx_vector_upload(this%bh, this%nrows, 1, x, ierr); INSIST(ierr == 0)
    call famgx_vector_set_zero(this%xh, this%nrows, 1, ierr); INSIST(ierr == 0)
    call stop_timer('memcpy-vector')

    !! Call the AmgX solver.
    call start_timer('solve')
    call famgx_solver_solve_with_0_initial_guess(this%solver, this%bh, this%xh, ierr)
    INSIST(ierr == 0)
    call famgx_solver_get_status(this%solver, stat, ierr)
    INSIST(ierr == 0)
    INSIST(stat == 0)
    call stop_timer('solve')

    !! Retrieve the solution vector from AmgX
    call start_timer('memcpy-vector')
    call famgx_vector_download(this%xh, x, ierr)
    INSIST(ierr == 0)
    call famgx_unpin_memory(x, ierr); INSIST(ierr == 0)
    call stop_timer('memcpy-vector')

    call stop_timer('precon-apply')

  end subroutine apply


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
  subroutine copy_to_ijmatrix(src, matrix)

    use,intrinsic :: iso_fortran_env, only: int64

    type(pcsr_matrix), intent(in) :: src
    type(amgx_obj), intent(inout) :: matrix

    integer :: ierr, nrows, nrows_global, nnz
    integer, allocatable :: row_offset(:), global_row_rank(:)
    !integer(int64), allocatable :: col_global(:)
    integer, allocatable :: col_global(:)

    !! TODO-WARN: Currently hardwired for serial runs.
    !! Indices are converted to C 0-indexing here.
    nrows  = src%graph%row_ip%onP_size()
    nrows_global = src%graph%row_ip%global_size()
    nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    allocate(global_row_rank(nrows_global), col_global(nnz), row_offset(nrows+1))
    global_row_rank = 0
    col_global = src%graph%row_ip%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1)) - 1
    row_offset = src%graph%xadj(:nrows+1) - 1

    call famgx_pin_memory(row_offset, ierr); INSIST(ierr == 0)
    call famgx_pin_memory(col_global, ierr); INSIST(ierr == 0)
    call famgx_pin_memory(global_row_rank, ierr); INSIST(ierr == 0)
    call famgx_pin_memory(src%values, ierr); INSIST(ierr == 0)

    print *, "AmgX: upload matrix"
    ! call famgx_matrix_upload_all_global(matrix, nrows_global, nrows, nnz, 1, 1, &
    !     row_offset, col_global, src%values, &
    !     0, 0, global_row_rank, ierr) ! TODO-WARN: here set halo depth & num import rings to 1, 1 in parallel
    call famgx_matrix_upload_all(matrix, nrows, nnz, 1, 1, row_offset, col_global, src%values, ierr)
    INSIST(ierr == 0)

    call famgx_unpin_memory(src%values, ierr); INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix


  function config_string(this) result(config)

    type(pcsr_precon_amgx), intent(in) :: this
    character(:), allocatable :: config

    character(128) :: line

    config = "config_version = 2,"

    write (line,'(a,i2)') "max_iters = ", this%max_iter
    config = config // trim(line) // ","

    write (line,'(a,i2)') "print_solve_stats = ", min(this%print_level, 1)
    config = config // trim(line) // ","

    write (line,'(a,i2)') "print_grid_stats = ", max(this%print_level - 1, 0)
    config = config // trim(line) // ","

    write (line,'(a,i2)') "obtain_timings = ", this%logging_level
    config = config // trim(line) // ","

    config = config // &
        "monitor_residual = 0," // & ! run until max_iters is reached (default)
        "determinism_flag = 1," // & ! use deterministic algorithms
        !"smoother = PCG," // &
        "solver = AMG"               ! top-level solver (default)

  end function config_string

end module pcsr_precon_amgx_type
