!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use kinds, only: r8
  use input_utilities,  only: NULL_I
  use parameter_module, only: string_len
  use UbikSolve_module
  use truchas_logging_services
  implicit none
  private
 
  ! Public procedures
  public :: LINEAR_SOLVER, POST_SOLVE
 
  ! Public variables and data types.
  public :: Ubik_solver, Ubik_control_type, Ubik_type, Ubik_user
 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Number of Ubik_user elements that are used for default parameters.
  ! NOTE: remember to increment this if another default solver is added!
  integer, parameter, public :: DEFAULT_UBIK_CONTROLS = 4

  ! Position of default solver parameters in the Ubik_user array.
  integer, parameter, public :: UBIK_PRESSURE_DEFAULT     = 1
  integer, parameter, public :: UBIK_NK_DEFAULT           = 2
  integer, parameter, public :: UBIK_DISPLACEMENT_DEFAULT = 3
  integer, parameter, public :: UBIK_VISCOUS_DEFAULT      = 4

  ! Solver parameters
  integer, parameter, public :: SOLVER_NONE    = -1
  integer, parameter, public :: SOLVER_CG      = 0
  integer, parameter, public :: SOLVER_GMRES   = 1
  integer, parameter, public :: SOLVER_TFQMR   = 2
  integer, parameter, public :: SOLVER_BCGSTAB = 3
  integer, parameter, public :: SOLVER_DIRECT  = 4
  integer, parameter, public :: SOLVER_FGMRES  = 5

  ! Preconditioning parameters
  integer, parameter, public :: PRECOND_NONE     = 0
  integer, parameter, public :: PRECOND_JACOBI   = 1
  integer, parameter, public :: PRECOND_SSOR     = 2
  integer, parameter, public :: PRECOND_ILU0     = 3
  integer, parameter, public :: PRECOND_LU       = 4
  integer, parameter, public :: PRECOND_TM_SSOR  = 6
  integer, parameter, public :: PRECOND_TM_DIAG  = 7
  integer, parameter, public :: PRECOND_DIAGONAL = 8

  ! Scope parameters.
  integer, parameter, public :: PRECOND_SCOPE_GLOBAL = 0
  integer, parameter, public :: PRECOND_SCOPE_LOCAL  = 1

  ! Parameter to tell linear solver that it is part of a nonlinear NK solve.
  integer, parameter, public :: COMING_FROM_INSIDE_NK = NULL_I

  ! Define a UbikSolve derived type which will contain control
  ! parameters and arrays needed by UbikSolve.
  type Ubik_type
 
     ! status flag
     !
     ! with JTpack this was used for JTpack solver status as well as for
     ! signals within Telluride such as whether the linear solve is part
     ! of a nonlinear solve - now it is used only for Telluride signals,
     ! since solver status is now encapsulated with in the control type
     integer :: status
 
     ! solver type
     integer :: solver
 
     ! Name for this particular instance of Ubik_type
     character(string_len) :: name
 
     ! Use preconditioning? If so, then specify solution
     ! method for the preconditioning equation.
     integer :: precond
 
     ! number of preconditioning iterations
     integer :: precond_iter
 
     ! scope flag - local (on processor) or global (all processors)
     integer :: precond_scope

     ! factor flag -
     ! TRUE means ILU or LU should factor the preconditioning matrix,
     ! FALSE means just do the back solve
     logical :: factor

     ! UbikSolve control parameters
     type(Ubik_control_type) :: control

  end type Ubik_type
 
  ! Declare an instance of the Ubik type for the solver.
  type(Ubik_type), save :: Ubik_solver
  
  ! Declare an array of the Ubik type to hold
  ! user-input linear solution parameters
  type(Ubik_type), pointer, dimension(:), save :: Ubik_user
 
  ! Number of user-specified linear solution control parameters.
  integer, public, save :: linear_solutions

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
 
    ! Arguments
    real(r8), dimension(:), intent(INOUT), target :: Solution
    real(r8), dimension(:), intent(INOUT) :: RHS
    type(Ubik_type), intent(INOUT) :: Ubik

    ! Interface blocks for external routines:
    !   RESIDUAL
    !   MATVEC
    !   PRECONDITIONER
    !   PRECONDITIONER_UPDATE
#include "solver_function_prototypes.fpp"

    ! Local Variables (for residual visualization)
    real(r8), dimension(:), allocatable :: Residuals
    type(Ubik_vector_type), target :: Solution_vec

    ! Local Variables
    integer :: nunknowns, status

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

  END SUBROUTINE LINEAR_SOLVER
 
  SUBROUTINE POST_SOLVE (Ubik, routine)
    !======================================================================
    ! Purpose:
    !
    !   take appropriate action after a linear solve has been completed
    !   (successfully or unsuccessfully)
    !======================================================================

    ! Arguments
    type(Ubik_type), intent(INOUT) :: Ubik
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

  END SUBROUTINE POST_SOLVE

END MODULE LINEAR_SOLUTION
