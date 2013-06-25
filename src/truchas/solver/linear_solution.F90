#include "f90_assert.fpp"

MODULE LINEAR_SOLUTION
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures and quantities for the linear system solver.
  !
  !   Public Interface:
  !
  !     call LINEAR_SOLVER (Solution, RHS, Ubik, MATVEC, PRECONDITIONER)
  !
  !         call the linear system solver; return the result in Solution.
  !
  ! Contains: LINEAR_SOLVER
  !
  ! Author(s): The Telluridians (telluride-info@lanl.gov)
  !
  !=======================================================================
  use constants_module, only: ipreset
  use kind_module,      only: int_kind, real_kind, log_kind
  use parameter_module, only: string_len
  use UbikSolve_module

  implicit none

  ! Default private module
  private
 
  ! Public procedures
  public :: LINEAR_SOLVER, POST_SOLVE
 
  ! Public variables and data types.
  public :: Ubik_solver, Ubik_control_type, Ubik_type, Ubik_user
 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Number of Ubik_user elements that are used for default parameters.
  ! NOTE: remember to increment this if another default solver is added!
  integer(int_kind), parameter, public :: DEFAULT_UBIK_CONTROLS = 9

  ! Position of default solver parameters in the Ubik_user array.
  integer(int_kind), parameter, public :: UBIK_PRESSURE_DEFAULT     = 1
  integer(int_kind), parameter, public :: UBIK_ENERGY_DEFAULT       = 2
  integer(int_kind), parameter, public :: UBIK_ENERGY_NK_DEFAULT    = 3
  integer(int_kind), parameter, public :: UBIK_NK_DEFAULT           = 4
  integer(int_kind), parameter, public :: UBIK_DISPLACEMENT_DEFAULT = 5
  integer(int_kind), parameter, public :: UBIK_VISCOUS_DEFAULT      = 6 
  integer(int_kind), parameter, public :: UBIK_VIEWFACTOR_DEFAULT   = 7 
  integer(int_kind), parameter, public :: UBIK_SENSITIVITY_DEFAULT   = 8

  ! Solver parameters
  integer(int_kind), parameter, public :: SOLVER_NONE    = -1
  integer(int_kind), parameter, public :: SOLVER_CG      = 0
  integer(int_kind), parameter, public :: SOLVER_GMRES   = 1
  integer(int_kind), parameter, public :: SOLVER_TFQMR   = 2
  integer(int_kind), parameter, public :: SOLVER_BCGSTAB = 3
  integer(int_kind), parameter, public :: SOLVER_DIRECT  = 4
  integer(int_kind), parameter, public :: SOLVER_FGMRES  = 5

  ! Preconditioning parameters
  integer(int_kind), parameter, public :: PRECOND_NONE     = 0
  integer(int_kind), parameter, public :: PRECOND_JACOBI   = 1
  integer(int_kind), parameter, public :: PRECOND_SSOR     = 2
  integer(int_kind), parameter, public :: PRECOND_ILU0     = 3
  integer(int_kind), parameter, public :: PRECOND_LU       = 4
  integer(int_kind), parameter, public :: PRECOND_2LEVEL   = 5
  integer(int_kind), parameter, public :: PRECOND_TM_SSOR  = 6
  integer(int_kind), parameter, public :: PRECOND_TM_DIAG  = 7
  integer(int_kind), parameter, public :: PRECOND_DIAGONAL = 8

  ! Scope parameters.
  integer(int_kind), parameter, public :: PRECOND_SCOPE_GLOBAL = 0
  integer(int_kind), parameter, public :: PRECOND_SCOPE_LOCAL  = 1

  ! Parameter to tell linear solver that it is part of a nonlinear NK solve.
  integer(int_kind), parameter, public :: COMING_FROM_INSIDE_NK = ipreset

  ! Define a UbikSolve derived type which will contain control
  ! parameters and arrays needed by UbikSolve.
  type Ubik_type
 
     ! status flag
     !
     ! with JTpack this was used for JTpack solver status as well as for
     ! signals within Telluride such as whether the linear solve is part
     ! of a nonlinear solve - now it is used only for Telluride signals,
     ! since solver status is now encapsulated with in the control type
     integer(int_kind) :: status
 
     ! solver type
     integer(int_kind) :: solver
 
     ! Name for this particular instance of Ubik_type
     character(string_len) :: name
 
     ! Use preconditioning? If so, then specify solution
     ! method for the preconditioning equation.
     integer(int_kind) :: precond
 
     ! preconditioner used for the preconditioning equation
     integer(int_kind) :: precond_pre
 
     ! number of preconditioning iterations
     integer(int_kind) :: precond_iter
 
     ! scope flag - local (on processor) or global (all processors)
     integer(int_kind) :: precond_scope

     ! factor flag -
     ! TRUE means ILU or LU should factor the preconditioning matrix,
     ! FALSE means just do the back solve
     logical(log_kind) :: factor

     ! UbikSolve control parameters
     type(Ubik_control_type) :: control

  end type Ubik_type
 
  ! Declare an instance of the Ubik type for the solver.
  type(Ubik_type), save :: Ubik_solver
  
  ! Declare an array of the Ubik type to hold
  ! user-input linear solution parameters
  type(Ubik_type), pointer, dimension(:), save :: Ubik_user
 
  ! Number of user-specified linear solution control parameters.
  integer(int_kind), public, save :: linear_solutions

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS
 
  SUBROUTINE LINEAR_SOLVER (Solution, RHS, Ubik, MATVEC, PRECONDITIONER)
    !======================================================================
    ! Purpose:
    !
    !   Linear system solver driver; return the result in Solution
    !
    !======================================================================
    use debug_control_data
    use kind_module,     only: int_kind, real_kind
    use UbikSolve_module
    use truchas_logging_services
#ifdef USE_TBROOK
    use output_data_module, only: enable_tbrook_output
#endif
 
    ! Arguments
    real(real_kind), dimension(:), intent(INOUT), target :: Solution
    real(real_kind), dimension(:), intent(INOUT) :: RHS
    type(Ubik_type), intent(INOUT) :: Ubik

    ! Interface blocks for external routines:
    !   RESIDUAL
    !   MATVEC
    !   PRECONDITIONER
    !   PRECONDITIONER_UPDATE
#include "solver_function_prototypes.fpp"

    ! Local Variables (for residual visualization)
    real(real_kind), dimension(:), allocatable :: Residuals
    type(Ubik_vector_type), target :: Solution_vec

    ! Local Variables
    integer(int_kind) :: nunknowns, status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get number of unknowns for this solution.
    nunknowns = SIZE(Solution,1)

    Ubik_solver = Ubik
    call Ubik_create (Ubik_solver%control)

    ! setup the preconditioning options
    PRECONDITIONER_SETUP: if (Ubik_solver%precond /= PRECOND_NONE) then

       Ubik_solver%precond_iter = 0
       Ubik_solver%Factor = .true.

    end if PRECONDITIONER_SETUP


    if (Ubik_outmode_summary(Ubik_solver%control) .or.  &
        Ubik_outmode_iterates(Ubik_solver%control) .or. &
        Ubik_outmode_full(Ubik_solver%control)) then
       call TLS_info (' Linear solver name: ' // trim(Ubik_Solver%name))
    end if
 
    ! call the appropriate solver
    select case(Ubik_solver%solver)
    case (SOLVER_GMRES)

       ! UbikSolve GMRES
       call Ubik_GMRES (Solution, RHS, Ubik_solver%control, MATVEC, PRECONDITIONER)

    case (SOLVER_FGMRES)

       ! UbikSolve FGMRES
       call Ubik_FGMRES (Solution, RHS, Ubik_solver%control, MATVEC, PRECONDITIONER)

    case (SOLVER_CG)

       ! UbikSolve CG
       call Ubik_CG (Solution, RHS, Ubik_solver%control, MATVEC, PRECONDITIONER)

    case (SOLVER_TFQMR)

       ! UbikSolve TFQMR
       call Ubik_TFQMR (Solution, RHS, Ubik_solver%control, MATVEC, PRECONDITIONER)

    case (SOLVER_BCGSTAB)

       ! UbikSolve Bi-CGSTAB
       call Ubik_BCGSTAB (Solution, RHS, Ubik_solver%control, MATVEC, PRECONDITIONER)

    case DEFAULT

       ! error
       call TLS_panic ('LINEAR_SOLVER: invalid linear solution method chosen in ' // trim(Ubik_solver%name))

    end select

    ! dump residuals in GMV format for visualization if solve failed
    if (.not.Ubik_converged(Ubik_solver%control)) then

       ! temporary array for residuals
       allocate (Residuals(nunknowns), stat=status)
       call TLS_fatal_if_any (status /= 0, 'LINEAR_SOLVER: error allocating array Residuals in ' // trim(Ubik_solver%name)) 

       ! MATVEC takes a Ubik_vector_type as the vector to be multiplied,
       ! so set one up for Solution
       call Ubik_create (Solution_vec, overlap_only=.true.)
       call Ubik_set_values_ptr(Solution_vec, Solution)

       ! compute Residuals = RHS - Coeff*Solution
       call MATVEC (Solution_vec, Residuals, status)
       Residuals = RHS - Residuals

#ifdef USE_TBROOK
       ! write residual and solution iterate to the xml file
       if (enable_tbrook_output) call xml_write_residual (Ubik_solver%name, Residuals, Solution)

#endif
       ! cleanup
       deallocate (Residuals)
    end if

    ! Check for errors in solver; print diagnostics and quit
    call POST_SOLVE (Ubik_solver, 'LINEAR_SOLVER')

    ! copy Ubik_solver out to Ubik before deallocating
    Ubik%status = Ubik_solver%status
    call Ubik_set_iter (Ubik%control, Ubik_iter(Ubik_solver%control))
    Ubik%precond_iter = Ubik_solver%precond_iter
    call Ubik_destroy (Ubik_solver%control)

    return
  END SUBROUTINE LINEAR_SOLVER
 
  SUBROUTINE POST_SOLVE (Ubik, routine)
    !======================================================================
    ! Purpose:
    !
    !   take appropriate action after a linear solve has been completed
    !   (successfully or unsuccessfully)
    !======================================================================
    use UbikSolve_module
    use truchas_logging_services

    implicit none

    ! Arguments
    type (Ubik_type), intent(INOUT) :: Ubik
    character(*), intent(IN) :: routine
    
    character(256) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (Ubik_itmax_exceeded(Ubik%control)) then

       call TLS_info (' UbikSolve: Maximum iterations exceeded')
       call TLS_info ('   Linear solver name: ' // trim(Ubik%name))
       ! Don't punt if we're inside of NK; just warn the user.
       if (Ubik%status == COMING_FROM_INSIDE_NK) then
          call TLS_info ('     since this is a linear solve within a nonlinear Newton iteration,&
              & the nonlinear iteration will continue')
          write(message,'(7x,a,es11.4)') 'residual norm      : ', Ubik_rnorm(Ubik%control)
          call TLS_info (message)
          write(message,'(7x,a,es11.4)') 'error estimate norm: ', Ubik_err(Ubik%control)
          call TLS_info (message)
       else
          call TLS_fatal (trim(routine) // ': UbikSolve: maximum iterations exceeded for linear solver ' &
                // trim(Ubik%name) // '; see the .err file for more info')
       end if

    else if (.not.Ubik_converged(Ubik%control)) then

       ! other kinds of failures
       if (Ubik_invalid_input(Ubik%control)) then
          message = ': UbikSolve: Invalid arguments'
       else if (Ubik_alloc_failure(Ubik%control)) then
          message = ': UbikSolve: Memory allocation failure'
       else if (Ubik_internal_error(Ubik%control)) then
          message = ': UbikSolve: Internal error'
       else if (Ubik_breakdown(Ubik%control)) then
          message = ': UbikSolve: Breakdown'
       else if (Ubik_matmul_error(Ubik%control)) then
          message = ': UbikSolve: Failure in matrix-vector multiplication'
       else if (Ubik_precond_error(Ubik%control)) then
          message = ': UbikSolve: Failure in preconditioner'
       else
          message = ': UbikSolve: Unknown failure'
       end if
       call TLS_fatal (trim(routine) // trim(message) // ' for linear solver ' &
             // trim(Ubik%name) // '; see the .err file for more info')

    end if

    return
  END SUBROUTINE POST_SOLVE

#ifdef USE_TBROOK
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! XML_WRITE_RESIDUAL
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 19 Jul 2005
 !!
 !! Write the passed residual and solution vectors to the XML output file.
 !! This creates a 'LINEAR_RESIDUAL' xml element in the output file; the data
 !! itself is written to a binary lookaside file.
 !!
 !! NB: This is extremely similar to a procedure of the same name from the
 !! NONLINEAR_SOLUTION module; the two ought to be consolidated into a
 !! common code base somehow.
 !!

  subroutine xml_write_residual (name, r, x)

    use brook_module
    use tbrook_module
    use output_module, only: prefix
    use string_utilities, only: i_to_c
    use parameter_module, only: ncells, nnodes

    character(len=*), intent(in) :: name
    real(real_kind),  intent(in) :: r(:), x(:)

    integer :: status, dim
    integer, save :: df_num = 0
    character(len=256) :: df_name
    type(brook), target :: df_brook
    real(real_kind), allocatable :: tmp(:,:)

    ASSERT( size(r) == size(x) )

    status = 0 ! for some insane reason, this is intent(in) for all the tbrook stuff.

    df_num = df_num + 1
    df_name = trim(prefix) // '.linear_res.' // i_to_c(df_num) // '.bin'

    !! Create the binary look-aside file.
    call tbrook_set (df_brook, file=trim(df_name), form='binary', istatus=status)
    if (status /= 0) return

    !! Open the LINEAR_RESIDUAL tag.
    call tbrook_openxmltag (BaseBrook, XMLTag='LINEAR_RESIDUAL', &
        XMLAttributes='SEQ="' // i_to_c(df_num) // '" SOLVER="' // trim(name) // '"', &
        istatus=status)
        if (status /= 0) return

    !! Write the FILE tag which specifies the look-aside data file.
    call tbrook_writexmltag (BaseBrook, XMLTag='FILE', &
        XMLAttributes='FORMAT="binary"', XMLStringData=trim(df_name), istatus=status)
        if (status /= 0) return

    !! Write the residual data.
    if (mesh_based_scalar_data(size(r),ncells)) then ! scalar, cell-based data
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', r, status, map='cell')
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', x, status, map='cell')
    else if (mesh_based_scalar_data(size(r),nnodes)) then ! scalar, node-based data
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', r, status, map='node')
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', x, status, map='node')
    else if (mesh_based_vector_data(size(r),ncells,dim)) then ! vector, cell-based data
      allocate(tmp(dim,ncells))
      call copy_to_rank_2 (r, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', tmp, status, map='cell')
      call copy_to_rank_2 (x, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', tmp, status, map='cell')
      deallocate(tmp)
    else if (mesh_based_vector_data(size(r),nnodes,dim)) then ! vector, node-based data
      allocate(tmp(dim,nnodes))
      call copy_to_rank_2 (r, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', tmp, status, map='node')
      call copy_to_rank_2 (x, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', tmp, status, map='node')
      deallocate(tmp)
    else ! we have no clue what it is ...
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', r, status, map='none')
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', x, status, map='none')
    end if
    if (status /= 0) return

    !! Close the LINEAR_RESIDUAL tag.
    call tbrook_closexmltag (BaseBrook, XMLTag='LINEAR_RESIDUAL', istatus=status)
        if (status /= 0) return

    !! Close the binary look-aside file.
    call tbrook_close   (df_brook, istatus=status)
    call tbrook_destroy (df_brook, istatus=status)

  contains

    !!
    !! Auxillary functions that infer whether the data, which is stored in a
    !! rank-1 array, is mesh-based scalar or vector data.  NDATA is the number
    !! of data elements, and NMESH is the number of mesh objects -- either
    !! NCELLS or NNODES, typically.  For vector data we need to be careful as
    !! it is possible for the number of nodes/cells on a processor to be zero.
    !! These are parallel procedures, returning global results.  The vector
    !! procedure returns the dimension of the vector data in DIM, which is
    !! only meaningful when the function returns the value true.
    !!

    logical function mesh_based_scalar_data (ndata, nmesh)
      use pgslib_module, only: pgslib_global_all
      integer, intent(in) :: ndata, nmesh
      mesh_based_scalar_data = pgslib_global_all(ndata == nmesh)
    end function mesh_based_scalar_data

    logical function mesh_based_vector_data (ndata, nmesh, dim)
      use pgslib_module, only: pgslib_global_all, pgslib_global_maxval
      integer, intent(in) :: ndata, nmesh
      integer, intent(out) :: dim
      integer :: d
      if (nmesh == 0) then
        mesh_based_vector_data = (ndata == 0)
        d = 0
      else
        mesh_based_vector_data = (modulo(ndata,nmesh) == 0)
        d = ndata / nmesh
      end if
      mesh_based_vector_data = pgslib_global_all(mesh_based_vector_data)
      if (.not.mesh_based_vector_data) return
      dim = pgslib_global_maxval(d)
      mesh_based_vector_data = (dim > 1) .and. pgslib_global_all(d == dim .or. d == 0)
    end function mesh_based_vector_data

    !!
    !! Auxillary subroutine to copy the contents of a rank-1 array into a
    !! same-sized rank-2 array.  Works for zero-sized arrays too.
    !!

    subroutine copy_to_rank_2 (in, out)
      real(real_kind), intent(in)  :: in(:)
      real(real_kind), intent(out) :: out(:,:)
      integer :: i, j, n
      ASSERT( size(in) == size(out) )
      n = 0
      do j = 1, size(out,2)
        do i = 1, size(out,1)
          n = n + 1
          out(i,j) = in(n)
        end do
      end do
    end subroutine copy_to_rank_2

  end subroutine xml_write_residual

#endif
END MODULE LINEAR_SOLUTION
