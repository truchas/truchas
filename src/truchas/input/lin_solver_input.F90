MODULE LIN_SOLVER_INPUT
  !=======================================================================
  ! Purpose(s):
  !
  !   Define various procedures for the input of global linear and 
  !   nonlinear solution algorithm flags and parameters.
  !
  !   Public Interface:
  !
  !     * call LINEAR_SOLVER_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the LINEAR_SOLVER namelist.
  !
  ! Contains: LINEAR_SOLVER_INPUT
  !           LINEAR_SOLVER_CHECK
  !           LINEAR_SOLVER_DEFAULT
  !           LINEAR_SOLVER_INPUT_PARALLEL
  !           SET_UBIK
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use input_utilities,  only: NULL_I, NULL_R, NULL_C
  use parameter_module, only: string_len
  use truchas_logging_services
  implicit none
  private

  public :: LINEAR_SOLVER_INPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! LINEAR_SOLVER Namelist input variables.

  ! Name for this solver.
  character(string_len) :: name

  ! Solution method.
  character(string_len) :: method

  ! Solver preconditioner.
  character(string_len) :: preconditioning_method

  ! Preconditioner preconditioner.
  character(string_len) :: preconditioning_preconditioner

  ! Preconditioner steps.
  integer :: preconditioning_steps

  ! Solver/preconditioner scope.
  character(string_len) :: preconditioning_scope

  ! Maximum iterations.
  integer :: maximum_iterations

  ! Stopping test.
  character(string_len) :: stopping_criterion 

  ! Solver output unit number.
  character(string_len) :: output_mode

  ! Frequency of solver status reports to the tty.
  integer :: status_frequency

  ! Convergence criterion.
  real(r8) :: convergence_criterion

  ! Krylov subspace vector size.
  integer :: krylov_vectors

  ! Preconditioner relaxation parameter (for Jacobi or SSOR).
  real(r8) :: relaxation_parameter

  ! Derived quantities.

  ! linear solution stopping criterion flag
  integer :: linear_solve_stop

  ! linear solution output mode
  integer :: linear_solve_output

CONTAINS

  SUBROUTINE LINEAR_SOLVER_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read LINEAR_SOLVER namelist.
    !
    !=======================================================================
    use input_utilities,        only: seek_to_namelist
    use linear_solution,        only: Ubik_user, DEFAULT_UBIK_CONTROLS, linear_solutions
    use parallel_info_module,   only: p_info
    use parameter_module,       only: string_len, string_dim
    use pgslib_module,          only: PGSLIB_BCAST

    integer, intent(in) :: lun

    ! Local Variables
    character(string_len) :: line
    logical :: fatal, found
    integer :: ioerror, i, j
    character(128) :: message

    ! Define NUMERICS namelist.
    namelist /LINEAR_SOLVER/ convergence_criterion, maximum_iterations,      &
                             output_mode, method, preconditioning_method,    &
                             relaxation_parameter, name, stopping_criterion, &
                             status_frequency, krylov_vectors,               &
                             preconditioning_steps, preconditioning_scope,   &
                             preconditioning_preconditioner

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant variables.
    fatal = .false.
    linear_solutions   = 0

    ! Inform user we're about to read.
    call TLS_info ('')
    call TLS_info (' Reading LINEAR_SOLVER Namelist(s) ...')

    ! First count the namelists.
    if (p_info%IOP) then
      rewind lun
       LS_NAMELIST_COUNT: do 
          call seek_to_namelist (lun, 'LINEAR_SOLVER', found)
          if (found) then
             linear_solutions = linear_solutions + 1
             read (lun, *) line
          else
             exit LS_NAMELIST_COUNT
          end if
       end do LS_NAMELIST_COUNT
       rewind lun
    end if

    ! Broadcast the linear_solutions counter.
    if (.not. p_info%UseGlobalServices) call PGSLIB_BCAST (linear_solutions)

    ! Set defaults and/or allocate control structure.
    if (linear_solutions <= 0) then
       ! None found; write conditional warning.
       write (message, 15) 
15     format ('LINEAR_SOLVER namelist not found! Using defaults.')
       call TLS_warn (message)
       ALLOCATE (Ubik_user(DEFAULT_UBIK_CONTROLS))
    else
       ALLOCATE (Ubik_user(DEFAULT_UBIK_CONTROLS + linear_solutions))
    end if

    ! Allocate the default linear solution controls.
    do i = 1,DEFAULT_UBIK_CONTROLS
       call SET_UBIK (Ubik_user(i), i)
    end do
   
    ! Loop until all linear solver namelists are found.
    LS_NAMELIST_LOOP: do i = 1,linear_solutions

       ! Running total for control structures.
       j = i + DEFAULT_UBIK_CONTROLS

       ! Default variables in this namelist.
       call LINEAR_SOLVER_DEFAULT ()

       ! Read the namelist only on the I/O Root PE.
       IO_PE_ONLY: if (p_info%IOP) then

          ! Find namelist.
          call seek_to_namelist (lun, 'LINEAR_SOLVER', found)
          if (found) then

             ! Read namelist; abort if errors occur.
             read (lun, NML = LINEAR_SOLVER, IOSTAT = ioerror)
             fatal = (ioerror /= 0)
             ! Inform user that the LINEAR_SOLVER namelist is being read.
             write (message, 10) i
10           format (9x,'Reading LINEAR_SOLVER Namelist #',i2,' ...')
             call TLS_info ('')
             call TLS_info (message)

          else   

             ! Done; exit loop.
             exit LS_NAMELIST_LOOP

          end if

       end if IO_PE_ONLY

       ! Punt if this namelist read was unsuccessful. 
       call TLS_fatal_if_any (fatal, 'error reading LINEAR_SOLVER namelist')

       ! Broadcast all variables in the LINEAR_SOLVER namelist.
       call LINEAR_SOLVER_INPUT_PARALLEL ()

       ! Check for fatal input errors; abort if found.
       call LINEAR_SOLVER_CHECK (fatal)
       call TLS_fatal_if_any (fatal, 'terminating execution due to previous input errors')

       ! Assign this namelist's parameters to a Ubik control structure.
       call SET_UBIK (Ubik_user(j), j)

    end do LS_NAMELIST_LOOP

  END SUBROUTINE LINEAR_SOLVER_INPUT

  SUBROUTINE LINEAR_SOLVER_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check LINEAR_SOLVER namelist input variables for obvious errors
    !   and inconsistencies.
    !
    !=======================================================================
    use parameter_module, only: ncells_tot, nnodes_tot, ndim
    use utilities_module, only: STRING_COMPARE
    use solid_mechanics_data, only: solid_mechanics

    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    logical :: strings_match, this_string_matches
    character(string_len) :: string, string_default
    integer :: i, krylov_vectors_max
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Convergence criterion.
    if (convergence_criterion /= NULL_R) then
       if (convergence_criterion <= 0 .or. &
           convergence_criterion > 0.1_r8) then
          write (message,'(a,es13.5)') 'Invalid convergence criterion: Must be <= 0.10 and > ', 0
          call TLS_error (message)
          fatal = .true.
       end if   
    end if

    ! Krylov vectors.
    if (krylov_vectors /= NULL_I) then
       if (krylov_vectors <= 0) then
          call TLS_error ('Invalid Krylov vector size: must be positive!')
          fatal = .true.
       end if
       krylov_vectors_max = ncells_tot
       if (solid_mechanics) krylov_vectors_max = nnodes_tot * ndim
       if (krylov_vectors > krylov_vectors_max) then
          ! Warning Message; overwrite krylov_vectors
          write (message,'(a,i0)') 'Krylov vector size too large! Setting size = ', krylov_vectors_max
          call TLS_warn (message)
          krylov_vectors = krylov_vectors_max
       end if
   end if
   
    ! Maximum iterations.
    if (maximum_iterations /= NULL_I) then
       if (maximum_iterations <= 0) then
          call TLS_error ('Invalid allowed maximum iterations: must be positive!')
          fatal = .true.
       end if
    end if
   
    ! Solution method.
    if (method /= NULL_C) then

       ! Check for valid strings.
       method = ADJUSTL(method)
       strings_match = .false.

       call STRING_COMPARE (TRIM(method), 'none', this_string_matches)
       if (this_string_matches) method = 'none'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'default', this_string_matches)
       if (this_string_matches) method = 'default'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'cg', this_string_matches)
       if (this_string_matches) method = 'cg'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'gmres', this_string_matches)
       if (this_string_matches) method = 'gmres'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'fgmres', this_string_matches)
       if (this_string_matches) method = 'fgmres'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'tfqmr', this_string_matches)
       if (this_string_matches) method = 'tfqmr'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'bcgstab', this_string_matches)
       if (this_string_matches) method = 'bcgstab'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'sparse_kit', this_string_matches)
       if (this_string_matches) method = 'sparse_kit'
       strings_match = strings_match .or. this_string_matches

       if (.not. strings_match) then
          call TLS_error ('Solution method "' // trim(method) //'" not valid!')
          fatal = .true.
       end if

    end if

    ! Solution name (can be any non-null character string).
    if (name == NULL_C) then
       call TLS_error ('Solution must be given a name (arbitrary string).')
       fatal = .true.
    else  
       ! Make sure the string is non-null.
       name = ADJUSTL(name)
       if (LEN_TRIM(name) == 0) then
          call TLS_error ('Solution name cannot be a null string!')
          fatal = .true.
       end if   
    end if

    ! Output mode.
    if (output_mode /= NULL_C) then

       ! Check for valid strings.
       output_mode = ADJUSTL(output_mode)
       strings_match = .false.

       call STRING_COMPARE (TRIM(output_mode), 'none', this_string_matches)
       if (this_string_matches) then
          output_mode = 'none'
          linear_solve_output = 0
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'errors', this_string_matches)
       if (this_string_matches) then
          output_mode = 'errors'
          linear_solve_output = 1
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'errors+warnings', this_string_matches)
       if (this_string_matches) then
          output_mode = 'errors+warnings'
          linear_solve_output = 2
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'warnings+errors', this_string_matches)
       if (this_string_matches) then
          output_mode = 'warnings+errors'
          linear_solve_output = 2
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'summary', this_string_matches)
       if (this_string_matches) then
          ! warnings, errors, plus a one-line summary consisting of the iteration
          ! number, the norms of the calculated residual and error estimate, and
          ! the norms of the true residual and error estimate (if calculated)
          ! each time convergence is checked
          output_mode = 'summary'
          linear_solve_output = 3
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'iterates', this_string_matches)
       if (this_string_matches) then
          ! same as summary, plus the coefficient, preconditioner (if there is
          ! one), source, initial guess and converged solution, plus the current
          ! iterate and true residual (if computed) at each iteration
          output_mode = 'iterates'
          linear_solve_output = 4
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'full', this_string_matches)
       if (this_string_matches) then
          ! same as iterates, plus vectors at each stage of the iteration, details
          ! about the preconditioner (if available), etc.
          ! WARNING: this generates reams of output
          output_mode = 'full'
          linear_solve_output = 5
       end if   
       strings_match = strings_match .or. this_string_matches

       if (.not. strings_match) then
          call TLS_error ('Output mode "' // trim(output_mode) // '" not valid!')
          fatal = .true.
       end if

    end if

    ! Frequency of status updates
    if (status_frequency /= NULL_I) then
       if (status_frequency < 0) then
          call TLS_error ('Invalid value for frequency of status updates: must be non-negative!')
          fatal = .true.
       end if
    end if

    ! Preconditioner and preconditioning preconditioner method.
    PRECONDITIONER_CHECK: do i = 1,2

       select case (i)
       case (1)
          string         = preconditioning_method
          string_default = NULL_C
       case (2)
          string         = preconditioning_preconditioner
          string_default = NULL_C
       end select

       if (string /= string_default) then

          ! Check for valid strings.
          string = ADJUSTL(string)
          strings_match = .false.

          call STRING_COMPARE (TRIM(string), 'none', this_string_matches)
          if (this_string_matches) string = 'none'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'diagonal', this_string_matches)
          if (this_string_matches) string = 'diagonal'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'jacobi', this_string_matches)
          if (this_string_matches) string = 'jacobi'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'ssor', this_string_matches)
          if (this_string_matches) string = 'ssor'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'ilu0', this_string_matches)
          if (this_string_matches) string = 'ilu0'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'lu', this_string_matches)
          if (this_string_matches) string = 'lu'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), '2level', this_string_matches)
          if (this_string_matches) string = '2level'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'tm_ssor', this_string_matches)
          if (this_string_matches) string = 'tm_ssor'
          strings_match = strings_match .or. this_string_matches

          call STRING_COMPARE (TRIM(string), 'tm_diag', this_string_matches)
          if (this_string_matches) string = 'tm_diag'
          strings_match = strings_match .or. this_string_matches

          if (.not. strings_match) then
             call TLS_error ('Preconditioning method "' // trim(string) // '" not valid!')
             fatal = .true.
          end if

       end if
   
    end do PRECONDITIONER_CHECK

    ! Preconditioning scope.
    if (preconditioning_scope /= NULL_C) then

       ! Check for valid strings.
       preconditioning_scope = ADJUSTL(preconditioning_scope)
       strings_match = .false.

       call STRING_COMPARE (TRIM(preconditioning_scope), 'global', this_string_matches)
       if (this_string_matches) preconditioning_scope = 'global'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(preconditioning_scope), 'local', this_string_matches)
       if (this_string_matches) preconditioning_scope = 'local'
       strings_match = strings_match .or. this_string_matches

       if (.not. strings_match) then
          call TLS_error ('Preconditioning scope "' // trim(preconditioning_scope) // '" not valid!')
          fatal = .true.
       end if

    end if

    ! Preconditioning steps.
    if (preconditioning_steps /= NULL_I) then
       if (preconditioning_steps <= 0) then
          call TLS_error ('Invalid number of preconditioning steps: must be positive!')
          fatal = .true.
       end if
    end if

    ! Relaxation parameter.
    if (relaxation_parameter /= NULL_R) then
       if (relaxation_parameter <= 0 .or. relaxation_parameter >= 2) then
          call TLS_error ('Relaxation parameter must be > 0.0 and < 2.0!')
          fatal = .true.
       end if
    end if

    ! Stopping criterion.
    if (stopping_criterion /= NULL_C) then

       ! Check for valid strings.
       stopping_criterion = ADJUSTL(stopping_criterion)
       strings_match = .false.

       call STRING_COMPARE (TRIM(stopping_criterion), '||x-xold||/||x||', this_string_matches)
       if (this_string_matches) then
          stopping_criterion = '||x-xold||/||x||'
          linear_solve_stop = 0
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(stopping_criterion), '||r||/(||A||*||x||+||b||)', this_string_matches)
       if (this_string_matches) then
          stopping_criterion = '||r||/(||A||*||x||+||b||)'
          linear_solve_stop = 1
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(stopping_criterion), '||r||/||b||', this_string_matches)
       if (this_string_matches) then
          stopping_criterion = '||r||/||b||'
          linear_solve_stop = 2
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(stopping_criterion), '||r||/||x||', this_string_matches)
       if (this_string_matches) then
          stopping_criterion = '||r||/||x||'
          linear_solve_stop = 3
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(stopping_criterion), '||r||/||r0||', this_string_matches)
       if (this_string_matches) then
          stopping_criterion = '||r||/||r0||'
          linear_solve_stop = 4
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(stopping_criterion), '||r||', this_string_matches)
       if (this_string_matches) then
          stopping_criterion = '||r||'
          linear_solve_stop = 5
       end if   
       strings_match = strings_match .or. this_string_matches

       if (.not. strings_match) then
          call TLS_error ('Stopping criterion "' // trim(stopping_criterion) // '" not valid!')
          fatal = .true.
       end if

    end if
   
  END SUBROUTINE LINEAR_SOLVER_CHECK

  SUBROUTINE LINEAR_SOLVER_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default LINEAR_SOLVER namelist variables to nonphysical values
    !   so that it is easy to tell if they have been input into the
    !   namelist.
    !
    !=======================================================================

    convergence_criterion          = NULL_R
    krylov_vectors                 = NULL_I
    maximum_iterations             = NULL_I
    method                         = NULL_C
    name                           = NULL_C
    output_mode                    = NULL_C
    preconditioning_method         = NULL_C
    preconditioning_preconditioner = NULL_C
    preconditioning_scope          = NULL_C
    preconditioning_steps          = NULL_I
    relaxation_parameter           = NULL_R
    stopping_criterion             = NULL_C
    status_frequency               = NULL_I

  END SUBROUTINE LINEAR_SOLVER_DEFAULT

  SUBROUTINE LINEAR_SOLVER_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast all elements of LINEAR_SOLVER namelist.
    !
    !======================================================================
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_BCAST

    if (.NOT. p_info%UseGlobalServices) then
       call PGSLib_BCAST (convergence_criterion)
       call PGSLib_BCAST (krylov_vectors)
       call PGSLib_BCAST (maximum_iterations)
       call PGSLib_BCAST (method)
       call PGSLib_BCAST (name)
       call PGSLib_BCAST (output_mode)
       call PGSLib_BCAST (preconditioning_method)
       call PGSLib_BCAST (preconditioning_preconditioner)
       call PGSLib_BCAST (preconditioning_scope)
       call PGSLib_BCAST (preconditioning_steps)
       call PGSLib_BCAST (relaxation_parameter)
       call PGSLib_BCAST (stopping_criterion)
       call PGSLib_BCAST (status_frequency)
    end if

  END SUBROUTINE LINEAR_SOLVER_INPUT_PARALLEL

  SUBROUTINE SET_UBIK (Ubik, solution)
    !=======================================================================
    ! Purpose:
    !
    !   Initialize parameters and arrays necessary for UbikSolve.
    !
    !   We select a default solver here, but always defer to the users wishes
    !=======================================================================
    use linear_solution,      only: Ubik_type, &
                                    SOLVER_NONE, &
                                    SOLVER_CG, SOLVER_GMRES, SOLVER_FGMRES,     &
                                    SOLVER_TFQMR, SOLVER_BCGSTAB,               &
                                    PRECOND_NONE, PRECOND_DIAGONAL,             &
                                    PRECOND_JACOBI, PRECOND_SSOR,               &
                                    PRECOND_ILU0, PRECOND_LU,                   &
                                    PRECOND_2LEVEL,                             &
                                    PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL,  &
                                    UBIK_PRESSURE_DEFAULT, UBIK_ENERGY_DEFAULT, &
                                    UBIK_ENERGY_NK_DEFAULT, UBIK_NK_DEFAULT,    &
                                    UBIK_DISPLACEMENT_DEFAULT, PRECOND_TM_DIAG, &
                                    UBIK_VIEWFACTOR_DEFAULT,                    &
                                    UBIK_SENSITIVITY_DEFAULT,                   &
                                    PRECOND_TM_SSOR, ubik_viscous_default
    use mesh_input_module,    only: mesh_file
    use parameter_module,     only: ncells_tot
    use two_level_partition,  only: precond_2level_active
    use UbikSolve_module
    use input_utilities, only: NULL_C

    ! Arguments
    type(Ubik_type), intent(INOUT) :: Ubik
    integer, intent(IN) :: solution
 
    ! Input variable defaults.
    real(r8), parameter :: CONVERGENCE_CRITERION_DEFAULT = 1.0d-8
    character(string_len), parameter :: METHOD_DEFAULT                  = 'fgmres'
    character(string_len), parameter :: NAME_DEFAULT                    = 'default'
    character(string_len), parameter :: PRECONDITIONING_METHOD_DEFAULT  = 'none'
    character(string_len), parameter :: PRECONDITIONING_PRECOND_DEFAULT = 'none'
    character(string_len), parameter :: PRECONDITIONING_SCOPE_DEFAULT   = 'global'
    integer,  parameter :: PRECONDITIONING_STEPS_DEFAULT   = 1
    real(r8), parameter :: RELAXATION_PARAMETER_DEFAULT    = 0.90
    integer,  parameter :: LINEAR_SOLVE_STOP_DEFAULT       = 5  ! ||r||
    integer,  parameter :: LINEAR_SOLVE_R_B_STOP_DEFAULT   = 2  ! ||r||/||b||
    integer,  parameter :: LINEAR_SOLVE_OUTPUT_DEFAULT     = 2
    integer,  parameter :: STATUS_FREQUENCY_DEFAULT        = 0
    integer             :: KRYLOV_VECTORS_DEFAULT
    integer             :: MAXIMUM_ITERATIONS_DEFAULT
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Set defaults that can't be parameters.
    MAXIMUM_ITERATIONS_DEFAULT = MAX(20.0,MIN(1000.0,2.0*SQRT(REAL(ncells_tot))))
    KRYLOV_VECTORS_DEFAULT     = MAX(10.0,MIN(100.0,SQRT(REAL(ncells_tot))))

    ! First default the UbikSolve control parameters.
    call Ubik_set_defaults (Ubik%Control)

    ! Set input variable defaults corresponding
    ! to default linear solution controls.
    select case (solution)
    case (UBIK_PRESSURE_DEFAULT)
       name                           = 'pressure'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_STOP_DEFAULT
       convergence_criterion          = 1.0d-10
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = 'jacobi'
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = 1.4
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (UBIK_ENERGY_DEFAULT)
       name                           = 'energy'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_STOP_DEFAULT
       convergence_criterion          = CONVERGENCE_CRITERION_DEFAULT
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = PRECONDITIONING_METHOD_DEFAULT
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (UBIK_ENERGY_NK_DEFAULT)
       name                           = 'energy_nk'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_STOP_DEFAULT
       convergence_criterion          = 1.0d-3
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = PRECONDITIONING_METHOD_DEFAULT
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (UBIK_NK_DEFAULT)
       name                           = 'nk'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_STOP_DEFAULT
       convergence_criterion          = 1.0d-3
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = PRECONDITIONING_METHOD_DEFAULT
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (UBIK_DISPLACEMENT_DEFAULT)
       name                           = 'displacement'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_R_B_STOP_DEFAULT
       convergence_criterion          = CONVERGENCE_CRITERION_DEFAULT
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = PRECONDITIONING_METHOD_DEFAULT
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (ubik_viscous_default)
       name                           = 'viscous'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_R_B_STOP_DEFAULT
       convergence_criterion          = CONVERGENCE_CRITERION_DEFAULT
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       ! Default preconditioning method should be diagonal for viscous terms
       preconditioning_method         = 'diagonal'
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (ubik_viewfactor_default)
       name                           = 'radiation_flux'
       method                         = 'gmres'
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_R_B_STOP_DEFAULT
       convergence_criterion          = CONVERGENCE_CRITERION_DEFAULT
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = PRECONDITIONING_METHOD_DEFAULT
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case (UBIK_SENSITIVITY_DEFAULT)
       name                           = 'energy_sensitivity'
       method                         = METHOD_DEFAULT
       linear_solve_output            = LINEAR_SOLVE_OUTPUT_DEFAULT
       linear_solve_stop              = LINEAR_SOLVE_STOP_DEFAULT
       convergence_criterion          = CONVERGENCE_CRITERION_DEFAULT
       maximum_iterations             = MAXIMUM_ITERATIONS_DEFAULT
       preconditioning_method         = PRECONDITIONING_METHOD_DEFAULT
       preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       preconditioning_scope          = PRECONDITIONING_SCOPE_DEFAULT
       preconditioning_steps          = PRECONDITIONING_STEPS_DEFAULT
       relaxation_parameter           = RELAXATION_PARAMETER_DEFAULT
       krylov_vectors                 = KRYLOV_VECTORS_DEFAULT
       status_frequency               = STATUS_FREQUENCY_DEFAULT
    case DEFAULT
       ! If variables haven't been provided, set defaults.
       if (method == NULL_C) then
          method = METHOD_DEFAULT
       end if
       if (output_mode == NULL_C) then
          linear_solve_output = LINEAR_SOLVE_OUTPUT_DEFAULT
       end if
       if (stopping_criterion == NULL_C) then
          linear_solve_stop = LINEAR_SOLVE_STOP_DEFAULT
       end if
       if (convergence_criterion == NULL_R) then
          convergence_criterion = CONVERGENCE_CRITERION_DEFAULT
       end if
       if (maximum_iterations == NULL_I) then
          maximum_iterations = MAXIMUM_ITERATIONS_DEFAULT
       end if
       if (preconditioning_method == NULL_C) then
          preconditioning_method = PRECONDITIONING_METHOD_DEFAULT
       end if
       if (preconditioning_preconditioner == NULL_C) then
          preconditioning_preconditioner = PRECONDITIONING_PRECOND_DEFAULT
       end if
       if (preconditioning_scope == NULL_C) then
          preconditioning_scope = PRECONDITIONING_SCOPE_DEFAULT
       end if
       if (preconditioning_steps == NULL_I) then
          preconditioning_steps = PRECONDITIONING_STEPS_DEFAULT
       end if
       if (relaxation_parameter == NULL_R) then
          relaxation_parameter = RELAXATION_PARAMETER_DEFAULT
       end if
       if (krylov_vectors == NULL_I) then
          krylov_vectors = KRYLOV_VECTORS_DEFAULT
       end if
       if (status_frequency == NULL_I) then
          status_frequency = STATUS_FREQUENCY_DEFAULT
       end if
    end select
   
    ! Name.
    Ubik%name = name

    ! Logical units for output.
    call Ubik_set_luout (Ubik%Control, TLS_logging_unit())
    call Ubik_set_luerr (Ubik%Control, TLS_logging_unit()) 
    call Ubik_set_lutty (Ubik%Control, TLS_logging_unit())
 
    ! Output mode.
    select case (linear_solve_output)
    case (0)
       call Ubik_set_outmode_none (Ubik%Control)
    case (1)
       call Ubik_set_outmode_errors (Ubik%Control)
    case (2)
       call Ubik_set_outmode_warnings (Ubik%Control)
    case (3)
       call Ubik_set_outmode_summary (Ubik%Control)
    case (4)
       call Ubik_set_outmode_iterates (Ubik%Control)
    case (5)
       call Ubik_set_outmode_full (Ubik%Control)
    end select

    ! Frequency of solver status reports to the tty.
    call Ubik_set_output_frequency (Ubik%Control, status_frequency)
 
    ! Maximum iterations.
    if (maximum_iterations > 0) then
       call Ubik_set_itmax (Ubik%Control, maximum_iterations)
    end if
 
    ! Stopping Test: default is ||r|| / ||b||
    ! - prevent use of stopping test requiring ||A||
    select case (linear_solve_stop)
    case (0)
       call Ubik_set_stopping_relchg (Ubik%Control)
    case (1)
       call TLS_fatal ('SET_UBIK: invalid stopping test')
!       call Ubik_set_stopping_axb (Ubik%Control)
    case (2)
       call Ubik_set_stopping_b (Ubik%Control)
    case (3)
       call Ubik_set_stopping_x (Ubik%Control)
    case (4)
       call Ubik_set_stopping_r0 (Ubik%Control)
    case (5)
       call Ubik_set_stopping_r (Ubik%Control)
    case default
       call TLS_fatal ('SET_UBIK: invalid stopping test')
    end select
 
    ! Max size of Krylov subspace vector stored and used in [F]GMRES.
    call Ubik_set_subspace (Ubik%Control, krylov_vectors)

    ! Update computed residual with b-Ax during iteration.
    call Ubik_set_residual_update (Ubik%Control)

    ! Preconditioning method.
    select case (TRIM(preconditioning_method))
    case ('none')
       Ubik%precond = precond_none
    case ('diagonal')
       Ubik%precond = precond_DIAGONAL
    case ('jacobi')
       Ubik%precond = precond_Jacobi
    case ('ssor')
       Ubik%precond = precond_SSOR
    case ('ilu0')
       Ubik%precond = precond_ILU0
    case ('lu')
       Ubik%precond = precond_LU
    case ('2level')
       Ubik%precond = precond_2Level
       ! This just lets us know that we can output so add'l info
       precond_2level_active = .TRUE.
    case ('tm_ssor')
       Ubik%precond = precond_TM_SSOR
    case ('tm_diag')
       Ubik%precond = precond_TM_DIAG
    case default
       call TLS_fatal ('SET_UBIK: invalid preconditioner choice')
    end select

    ! Number of preconditioning iterations.
    Ubik%precond_iter = 0

    ! prepreconditioners.
    select case (TRIM(preconditioning_preconditioner))
    case ('none')
       Ubik%precond_pre = precond_none
    case ('diagonal')
       Ubik%precond_pre = precond_DIAGONAL
    case ('jacobi')
       Ubik%precond_pre = precond_Jacobi
    case ('ssor')
       Ubik%precond_pre = precond_SSOR
    case ('ilu0')
       Ubik%precond_pre = precond_ILU0
    case ('lu')
       Ubik%precond_pre = precond_LU
    case ('tm_ssor')
       Ubik%precond_pre = precond_TM_SSOR
    case ('tm_diag')
       Ubik%precond_pre = precond_TM_DIAG
    case default
       call TLS_fatal ('SET_UBIK: invalid prepreconditioner choice')
    end select

    ! Select a solution method.
    if (TRIM(method) == 'default') then
       method = 'fgmres'
    end if

    select case (TRIM(method))
    case ('none')
       Ubik%solver = solver_NONE
    case ('cg')
       if (TRIM(mesh_file) /= NULL_C) then
          write (message, 10) 
10        format (9x,'Suggest using FGMRES for linear solutions when a mesh file is specified. User specified CG, using CG')
          call TLS_info (message)
       end if
       if (TRIM(preconditioning_method) == '2level') then
          write (message, 20) 
20        format (9x,'Suggest using FGMRES for linear solutions (2-level preconditioner specified). User specified CG, using CG')
          call TLS_info (message)
       end if
       Ubik%solver = solver_CG
#ifndef LINUX
       call Ubik_set_method (Ubik%control, Ubik_method_CG)
#endif
    case ('gmres')
       Ubik%solver = solver_GMRES
#ifndef LINUX
       call Ubik_set_method (Ubik%control, Ubik_method_GMRES)
#endif
    case ('fgmres')
       Ubik%solver = solver_FGMRES
#ifndef LINUX
       call Ubik_set_method (Ubik%control, Ubik_method_FGMRES)
#endif
    case ('tfqmr')
       Ubik%solver = solver_TFQMR
#ifndef LINUX
       call Ubik_set_method (Ubik%control, Ubik_method_TFQMR)
#endif
    case ('bcgstab')
       Ubik%solver = solver_BCGSTAB
#ifndef LINUX
       call Ubik_set_method (Ubik%control, Ubik_method_BCGSTAB)
#endif
    case default
       call TLS_fatal ('SET_UBIK: Invalid solution method choice')
    end select

    ! Number of preconditioner steps.
    if (preconditioning_steps > 1) then
       call Ubik_set_steps (Ubik%Control, preconditioning_steps)
    end if
 
    ! Convergence criteria.
    if (convergence_criterion > 0 .and. convergence_criterion < 0.1_r8) then
       call Ubik_set_eps (Ubik%Control, convergence_criterion)
    end if
 
    ! Preconditioner relaxation.
    if (relaxation_parameter > 0 .and. relaxation_parameter < 2) then
       call Ubik_set_omega (Ubik%Control, relaxation_parameter)
    end if
 
    ! Preconditioning scope. Default is global.
    Ubik%precond_scope = precond_scope_global
    call Ubik_set_scope_global (Ubik%Control)
    if (TRIM(preconditioning_scope) == 'local') then
       Ubik%precond_scope = precond_scope_local
       call Ubik_set_scope_local (Ubik%Control)
    end if

  END SUBROUTINE SET_UBIK

END MODULE LIN_SOLVER_INPUT
