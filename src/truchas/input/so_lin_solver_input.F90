MODULE SO_LIN_SOLVER_INPUT
  !=======================================================================
  ! Purpose(s):
  !
  !   Define various procedures for the input of global linear and 
  !   nonlinear solution algorithm flags and parameters.
  !
  !   Public Interface:
  !
  !     * call SO_LINEAR_SOLVER_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the SO_LINEAR_SOLVER namelist.
  !
  ! Contains: SO_LINEAR_SOLVER_INPUT
  !           SO_LINEAR_SOLVER_CHECK
  !           SO_LINEAR_SOLVER_DEFAULT
  !           SO_LINEAR_SOLVER_INPUT_PARALLEL
  !           SO_SET_UBIK
  !
  !
  !=======================================================================
  use constants_module, only: preset, ipreset
  use kind_module,      only: int_kind, real_kind
  use parameter_module, only: string_len
  use truchas_logging_services

  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: SO_SOLVER_INPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! SO_LINEAR_SOLVER Namelist input variables.

  ! Name for this solver.
  character(string_len) :: name

  ! Solution method.
  character(string_len) :: method

  ! Maximum iterations.
  integer(int_kind)     :: maximum_iterations

  ! Stopping test.
  character(string_len) :: stopping_criterion 

  ! Solver output unit number.
  character(string_len) :: output_mode

  ! Frequency of solver status reports to the tty.
  integer(int_kind)     :: status_frequency

  ! Convergence criterion.
  real(real_kind)       :: convergence_criterion

  ! Parameters signifying no input.
  real(real_kind),       parameter :: NO_REAL_INPUT    = -preset
  integer(int_kind),     parameter :: NO_INTEGER_INPUT = -ipreset
  character(string_len), parameter :: NO_CHAR_INPUT    = ''
  
  ! Derived quantities.

  ! linear solution stopping criterion flag
  integer(int_kind) :: linear_solve_stop

  ! linear solution output mode
  integer(int_kind) :: linear_solve_output

CONTAINS

  SUBROUTINE SO_SOLVER_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read LINEAR_SOLVER namelist.
    !
    !=======================================================================
    use input_utilities,        only: seek_to_namelist
    use kind_module,            only: int_kind, log_kind
    use parallel_info_module,   only: p_info
    use parameter_module,       only: string_len
    use pgslib_module,          only: PGSLIB_BCAST

    integer, intent(in) :: lun

    ! Local Variables
    character(string_len)                        :: line
    logical(log_kind)                            :: fatal, found
    integer(int_kind)                            :: ioerror, so_linear_solutions
    character(128) :: message

    ! Define NUMERICS namelist.
    namelist /SO_SOLVER/ convergence_criterion, maximum_iterations,      &
                         output_mode, name, stopping_criterion,          &  
                         status_frequency

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    ! Initialize relevant variables.
    fatal = .false.

    ! Inform user we're about to read.
    call TLS_info ('')
    call TLS_info (' Reading SO_SOLVER Namelist ...')

    ! First count the namelists -- there should be either one or none
    ! (in the case that there is no namelist we will use reasonable defaults)
    so_linear_solutions = 0
    if (p_info%IOP) then
       rewind lun
       SO_LS_NAMELIST_COUNT: do 
          call seek_to_namelist (lun, 'SO_SOLVER', found)
          if (found) then
             so_linear_solutions = so_linear_solutions + 1
             read (lun, *) line
          else
             exit SO_LS_NAMELIST_COUNT
          end if
       end do SO_LS_NAMELIST_COUNT
       rewind lun
    end if

    ! Broadcast the linear_solutions counter.
    if (.not. p_info%UseGlobalServices) call PGSLIB_BCAST (so_linear_solutions)

    ! Set defaults and/or allocate control structure.
    if (so_linear_solutions <= 0) then
       ! None found; write conditional warning.
       write (message, 15) 
15     format (9x,'SO_SOLVER namelist not found! Using defaults.')
       call TLS_info (message)
    else if (so_linear_solutions > 1) then
       ! bail out, there should be at most one name list
       WRITE(*,*) 'We should bail here...'
    end if

    ! Default variables in this namelist.
    call SO_SOLVER_DEFAULT ()
    
    call SO_SET_UBIK ()

    IO_PE_ONLY: if (p_info%IOP) then
       
       if (so_linear_solutions == 1) then
          ! Find namelist.
          call seek_to_namelist (lun, 'SO_SOLVER', found)
          if (found) then
             
             ! Read namelist; abort if errors occur.
             read (lun, NML = SO_SOLVER, IOSTAT = ioerror)
             fatal = (ioerror /= 0)
          end if
          
       end if
    end if IO_PE_ONLY
  
    if ( .not. fatal ) then
       ! Punt if this namelist read was unsuccessful. 
       call TLS_fatal_if_any (fatal, 'SO_SOLVER_INPUT: error reading SO_SOLVER namelist')
       
       ! Broadcast all variables in the SO_SOLVER namelist.
       call SO_SOLVER_INPUT_PARALLEL ()
       
       ! Check for fatal input errors; abort if found.
       call SO_SOLVER_CHECK (fatal)
       
       call TLS_fatal_if_any (fatal, 'SO_SOLVER_INPUT: SO_SOLVER namelist input error')
       
       ! Assign this namelist's parameters to a Ubik control structure.
       call SO_SET_UBIK ()
       
    end if
    
    return
  END SUBROUTINE SO_SOLVER_INPUT

  SUBROUTINE SO_SOLVER_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check LINEAR_SOLVER namelist input variables for obvious errors
    !   and inconsistencies.
    !
    !=======================================================================
    use constants_module, only: one_tenth
    use kind_module,      only: log_kind, real_kind
    use utilities_module, only: STRING_COMPARE

    implicit none

    ! Argument List
    logical(log_kind), intent(INOUT) :: fatal

    ! Local Variables
    logical(log_kind)     :: strings_match, this_string_matches
    real(real_kind)       :: x
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Convergence criterion.
    if (convergence_criterion /= NO_REAL_INPUT) then
       if (convergence_criterion <= EPSILON(x) .or. &
           convergence_criterion > one_tenth) then
          write (message, 1) EPSILON(x)
1         format ('Invalid convergence criterion: Must be <= 0.10 and >= ',1pe13.5)
          call TLS_error (message)
          fatal = .true.
       end if   
    end if

    ! Maximum iterations.
    if (maximum_iterations /= NO_INTEGER_INPUT) then
       if (maximum_iterations <= 0) then
          call TLS_error ('Invalid allowed maximum iterations: must be positive!')
          fatal = .true.
       end if
    end if
   
    ! Solution name (can be any non-null character string).
    if (name == NO_CHAR_INPUT) then
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
    if (output_mode /= NO_CHAR_INPUT) then

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
    if (status_frequency /= NO_INTEGER_INPUT) then
       if (status_frequency < 0) then
          call TLS_error ('Invalid value for frequency of status updates: must be non-negative!')
          fatal = .true.
       end if
    end if

    ! Stopping criterion.
    if (stopping_criterion /= NO_CHAR_INPUT) then

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
   
    return
  END SUBROUTINE SO_SOLVER_CHECK

  SUBROUTINE SO_SOLVER_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default LINEAR_SOLVER namelist variables to nonphysical values
    !   so that it is easy to tell if they have been input into the
    !   namelist.
    !
    !=======================================================================
    implicit none

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    convergence_criterion          = NO_REAL_INPUT
    maximum_iterations             = NO_INTEGER_INPUT
    name                           = NO_CHAR_INPUT
    output_mode                    = NO_CHAR_INPUT
    stopping_criterion             = NO_CHAR_INPUT
    status_frequency               = NO_INTEGER_INPUT

    return
  END SUBROUTINE SO_SOLVER_DEFAULT

  SUBROUTINE SO_SOLVER_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast all elements of LINEAR_SOLVER namelist.
    !
    !======================================================================
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_BCAST

    implicit none

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.NOT. p_info%UseGlobalServices) then
       call PGSLib_BCAST (convergence_criterion)
       call PGSLib_BCAST (maximum_iterations)
       call PGSLib_BCAST (name)
       call PGSLib_BCAST (output_mode)
       call PGSLib_BCAST (stopping_criterion)
       call PGSLib_BCAST (status_frequency)
    end if

    return
  END SUBROUTINE SO_SOLVER_INPUT_PARALLEL

  SUBROUTINE SO_SET_UBIK ()
    !=======================================================================
    ! Purpose:
    !
    !   Initialize parameters and arrays necessary for UbikSolve.
    !
    !   We select a default solver here, but always defer to the users wishes
    !=======================================================================
    use constants_module,     only: one_tenth, zero, ten_tominus8
    use kind_module,          only: int_kind, real_kind
    use parameter_module,     only: ncells_tot
    use UbikSolve_module
    use support_operators,      only: SO_Ubik

    implicit none

    ! Input variable defaults.
    real(KIND = real_kind),      parameter :: CONVERGENCE_CRITERION_DEFAULT   = ten_tominus8
    character(LEN = string_len), parameter :: METHOD_DEFAULT                  = 'cg'
    character(LEN = string_len), parameter :: NAME_DEFAULT                    = 'default'
    integer(KIND = int_kind),    parameter :: LINEAR_SOLVE_STOP_DEFAULT       = 5  ! ||r||
    integer(KIND = int_kind),    parameter :: LINEAR_SOLVE_OUTPUT_DEFAULT     = 0
    integer(KIND = int_kind),    parameter :: STATUS_FREQUENCY_DEFAULT        = 0
    integer(KIND = int_kind)               :: MAXIMUM_ITERATIONS_DEFAULT

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Set defaults that can't be parameters.
    MAXIMUM_ITERATIONS_DEFAULT = MAX(20.0,MIN(1000.0,2.0*SQRT(REAL(ncells_tot))))

    ! First default the UbikSolve control parameters.
    call Ubik_set_defaults (SO_Ubik)
    call Ubik_create(SO_Ubik)

    ! Set input variable defaults corresponding
    ! to default linear solution controls.

    name                           = 'so'

    if (method == NO_CHAR_INPUT) then
       method = METHOD_DEFAULT
    end if
    
    if (output_mode == NO_CHAR_INPUT) then
       linear_solve_output = LINEAR_SOLVE_OUTPUT_DEFAULT
    end if
    
    if (stopping_criterion == NO_CHAR_INPUT) then
       linear_solve_stop = LINEAR_SOLVE_STOP_DEFAULT
    end if
    
    if (convergence_criterion == NO_REAL_INPUT) then
       convergence_criterion = CONVERGENCE_CRITERION_DEFAULT
    end if
    
    if (maximum_iterations == NO_INTEGER_INPUT) then
       maximum_iterations = MAXIMUM_ITERATIONS_DEFAULT
    end if
    
    if (status_frequency == NO_INTEGER_INPUT) then
       status_frequency = STATUS_FREQUENCY_DEFAULT
    end if


    ! Logical units for output.
    call Ubik_set_luout (SO_Ubik, TLS_logging_unit())
    call Ubik_set_luerr (SO_Ubik, TLS_logging_unit())
    call Ubik_set_lutty (SO_Ubik, TLS_logging_unit())
 
    ! Output mode.
    select case (linear_solve_output)
    case (0)
       call Ubik_set_outmode_none (SO_Ubik)
    case (1)
       call Ubik_set_outmode_errors (SO_Ubik)
    case (2)
       call Ubik_set_outmode_warnings (SO_Ubik)
    case (3)
       call Ubik_set_outmode_summary (SO_Ubik)
    case (4)
       call Ubik_set_outmode_iterates (SO_Ubik)
    case (5)
       call Ubik_set_outmode_full (SO_Ubik)
    end select

    ! Frequency of solver status reports to the tty.
    call Ubik_set_output_frequency (SO_Ubik, status_frequency)
 
    ! Maximum iterations.
    if (maximum_iterations > 0) then
       call Ubik_set_itmax (SO_Ubik, maximum_iterations)
    end if
 
    ! Stopping Test: default is ||r|| / ||b||
    ! - prevent use of stopping test requiring ||A||
    select case (linear_solve_stop)
    case (0)
       call Ubik_set_stopping_relchg (SO_Ubik)
    case (1)
       call TLS_fatal ('SET_UBIK: invalid stopping test')
!       call Ubik_set_stopping_axb (SO_Ubik)
    case (2)
       call Ubik_set_stopping_b (SO_Ubik)
    case (3)
       call Ubik_set_stopping_x (SO_Ubik)
    case (4)
       call Ubik_set_stopping_r0 (SO_Ubik)
    case (5)
       call Ubik_set_stopping_r (SO_Ubik)
    case default
       call TLS_fatal ('SET_UBIK: invalid stopping test')
    end select
 

    ! Update computed residual with b-Ax during iteration.
    !call Ubik_set_residual_update (Ubik)

    call Ubik_set_method (SO_Ubik, Ubik_method_CG)

    ! Convergence criteria.
    if (convergence_criterion > zero .and. convergence_criterion < one_tenth) then
       call Ubik_set_eps (SO_Ubik, convergence_criterion)
    end if
 

    return
  END SUBROUTINE SO_SET_UBIK

END MODULE SO_LIN_SOLVER_INPUT
