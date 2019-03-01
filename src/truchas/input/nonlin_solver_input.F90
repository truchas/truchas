!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NONLIN_SOLVER_INPUT
  !=======================================================================
  ! Purpose(s):
  !
  !   Define various procedures for the input of global linear and 
  !   nonlinear solution algorithm flags and parameters.
  !
  !   Public Interface:
  !
  !     * call NONLINEAR_SOLVER_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the NONLINEAR_SOLVER namelist.
  !
  ! Contains: NONLINEAR_SOLVER_INPUT
  !           NONLINEAR_SOLVER_CHECK
  !           NONLINEAR_SOLVER_DEFAULT
  !           NONLINEAR_SOLVER_INPUT_PARALLEL
  !           SET_NK_CONTROLS
  !
  ! Author(s): Brian Vanderheyden (wbv@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use input_utilities,    only: NULL_I, NULL_R, NULL_C
  use nonlinear_solution, only: ndampers
  use parameter_module,   only: string_len
  use truchas_logging_services
  implicit none
  private

  public :: NONLINEAR_SOLVER_INPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! NONLINEAR_SOLVER Namelist input variables.

  ! Name for this solver.
  character(string_len) :: name

  ! Solution method
  character(string_len) :: method

  ! Linear solver name.
  character(string_len) :: linear_solver_name

  ! Output mode
  character(string_len) :: output_mode

  ! Flag for using NK damper.
  logical :: use_damper

  ! Maximum iterations.
  integer :: maximum_iterations

  ! NLK method --- maximum number of acceleration vectors
  integer :: NLK_Max_Vectors

  ! NLK method --- tolerance for rejection of vectors
  real(r8) :: NLK_Vector_Tolerance

  ! Convergence criterion.
  real(r8) :: convergence_criterion

  ! NK damper parameters.
  real(r8), dimension(ndampers) :: Damper_Parameters

  ! NK perturbation parameter.
  real(r8) :: perturbation_parameter

CONTAINS

  SUBROUTINE NONLINEAR_SOLVER_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read NONLINEAR_SOLVER namelist.
    !
    !=======================================================================
    use input_utilities,        only: seek_to_namelist
    use nonlinear_solution,     only: NKuser, DEFAULT_NK_CONTROLS, nonlinear_solutions
    use parallel_info_module,   only: p_info
    use parameter_module,       only: string_len
    use pgslib_module,          only: PGSLIB_BCAST
    use string_utilities,       only: i_to_c

    integer, intent(in) :: lun

    ! Local Variables
    character(string_len) :: line
    logical :: fatal, found
    integer :: ioerror, i, j

    ! Define NUMERICS namelist.
    namelist /NONLINEAR_SOLVER/ name, method, linear_solver_name, use_damper, &
                                maximum_iterations, convergence_criterion,    &
                                Damper_Parameters, perturbation_parameter,    &
                                output_mode, NLK_Max_Vectors, &
                                NLK_Vector_Tolerance

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant variables.
    fatal              = .false.
    nonlinear_solutions   = 0

    call TLS_info ('')
    call TLS_info ('Reading NONLINEAR_SOLVER Namelists ...')

    ! First count the namelists.
    if (p_info%IOP) then
       rewind lun
       NLS_NAMELIST_COUNT: do 
          call seek_to_namelist (lun, 'NONLINEAR_SOLVER', found)
          if (found) then
             nonlinear_solutions = nonlinear_solutions + 1
             read (lun, *) line
          else
             exit NLS_NAMELIST_COUNT
          end if
       end do NLS_NAMELIST_COUNT
       rewind lun
    end if

    ! Broadcast the linear_solutions counter.
    if (.not. p_info%UseGlobalServices) call PGSLIB_BCAST (nonlinear_solutions)

    ! Set defaults and/or allocate control structure.
    if (nonlinear_solutions <= 0) then
       ! None found; write conditional warning.
       call TLS_info ('  NONLINEAR_SOLVER namelist not found; using defaults.')
       ALLOCATE (NKuser(DEFAULT_NK_CONTROLS))
    else
       ALLOCATE (NKuser(DEFAULT_NK_CONTROLS + nonlinear_solutions))
    end if

    ! Allocate the default linear solution controls.
    do i = 1,DEFAULT_NK_CONTROLS
       call NONLINEAR_SOLVER_DEFAULT ()
       call SET_NK_CONTROLS (NKuser(i), i)
    end do
   
    ! Loop until all linear solver namelists are found.
    NLS_NAMELIST_LOOP: do i = 1,nonlinear_solutions

       ! Running total for control structures.
       j = i + DEFAULT_NK_CONTROLS

       ! Default variables in this namelist.
       call NONLINEAR_SOLVER_DEFAULT ()

       ! Read the namelist only on the I/O Root PE.
       IO_PE_ONLY: if (p_info%IOP) then

          ! Find namelist.
          call seek_to_namelist (lun, 'NONLINEAR_SOLVER', found)
          if (found) then
             call TLS_info ('  Reading NONLINEAR_SOLVER Namelist #' // i_to_c(i))

             ! Read namelist; abort if errors occur.
             read (lun, NML = NONLINEAR_SOLVER, IOSTAT = ioerror)
             fatal = (ioerror /= 0)

          else   

             ! Done; exit loop.
             exit NLS_NAMELIST_LOOP

          end if

       end if IO_PE_ONLY

       ! Punt if this namelist read was unsuccessful. 
       call TLS_fatal_if_any (fatal, 'error reading NONLINEAR_SOLVER namelist')

       ! Broadcast all variables in the NONLINEAR_SOLVER namelist.
       call NONLINEAR_SOLVER_INPUT_PARALLEL ()

       ! Check for fatal input errors; abort if found.
       call NONLINEAR_SOLVER_CHECK (fatal)
       call TLS_fatal_if_any (fatal, 'terminating execution due to previous errors')

       ! Assign this namelist's parameters to a Ubik control structure.
       call SET_NK_CONTROLS (NKuser(j), j)

    end do NLS_NAMELIST_LOOP

  END SUBROUTINE NONLINEAR_SOLVER_INPUT

  SUBROUTINE NONLINEAR_SOLVER_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check NONLINEAR_SOLVER namelist input variables for obvious errors
    !   and inconsistencies.
    !
    !=======================================================================
    use nonlinear_solution, only: ndampers
    use utilities_module,   only: STRING_COMPARE

    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    logical :: strings_match, this_string_matches
    real(r8) :: x
    integer :: i
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Solution method.
    if (method /= NULL_C) then

       ! Check for valid strings.
       method = ADJUSTL(method)
       strings_match = .false.

       call STRING_COMPARE (TRIM(method), 'default', this_string_matches)
       if (this_string_matches) method = 'default'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'nk', this_string_matches)
       if (this_string_matches) method = 'nk'
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(method), 'nlk', this_string_matches)
       if (this_string_matches) method = 'nlk'
       strings_match = strings_match .or. this_string_matches

       if (.not. strings_match) then
          call TLS_error ('Solution method "' // trim(method) // '" not valid!')
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

    ! Linear solver name (can be any non-null character string).
    if (linear_solver_name == NULL_C) then
       call TLS_error ('Linear solver must be given a valid name.')
       fatal = .true.
    else  
       ! Make sure the string is non-null.
       linear_solver_name = ADJUSTL(linear_solver_name)
       if (LEN_TRIM(name) == 0) then
          call TLS_error ('Linear solver name cannot be a null string!')
          fatal = .true.
       end if   
    end if

    ! Convergence criterion.
    if (convergence_criterion /= NULL_R) then
       if (convergence_criterion <= EPSILON(x) .or. &
           convergence_criterion > 10.0) then
          write (message, 1) EPSILON(x)
1         format ('Invalid convergence criterion: Must be <= 10 and >= ',1pe13.5)
          call TLS_error (message)
          fatal = .true.
       end if   
    end if

    ! Maximum iterations.
    if (maximum_iterations /= NULL_I) then
       if (maximum_iterations <= 0) then
          call TLS_error ('Invalid allowed maximum iterations: must be positive!')
          fatal = .true.
       end if
    end if
   
    ! NLK_Vector_Tolerance
    if (NLK_Vector_Tolerance /= NULL_R) then
       if (NLK_Vector_Tolerance <= 0.0 .or. &
           NLK_Vector_Tolerance >= 1.0) then
          call TLS_error ('NLK_Vector_Tolerance: Must be > 0 and < 1')
          fatal = .true.
       end if   
    end if

    ! NLK_Max_Vectors
    if (NLK_Max_Vectors /= NULL_I) then
       if (NLK_Max_Vectors < 1) then
          call TLS_error ('Invalid allowed NLK_Max_Vectors: must be positive!')
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
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'errors', this_string_matches)
       if (this_string_matches) then
          output_mode = 'errors'
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'errors+warnings', this_string_matches)
       if (this_string_matches) then
          output_mode = 'errors+warnings'
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'warnings+errors', this_string_matches)
       if (this_string_matches) then
          output_mode = 'warnings+errors'
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'verbose0', this_string_matches)
       if (this_string_matches) then
          output_mode = 'verbose0'
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'verbose1', this_string_matches)
       if (this_string_matches) then
          output_mode = 'verbose1'
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'verbose2', this_string_matches)
       if (this_string_matches) then
          output_mode = 'verbose2'
       end if   
       strings_match = strings_match .or. this_string_matches

       call STRING_COMPARE (TRIM(output_mode), 'verbose3', this_string_matches)
       if (this_string_matches) then
          output_mode = 'verbose3'
       end if   
       strings_match = strings_match .or. this_string_matches

       if (.not. strings_match) then
          call TLS_error ('Output mode "' // trim(output_mode) // '" not valid!')
          fatal = .true.
       end if

    end if

    ! Perturbation parameter.
    if (perturbation_parameter /= NULL_R) then
       if (perturbation_parameter <= 0 .or. perturbation_parameter >= 1) then
          call TLS_error ('Perturbation parameter must be > 0.0 and < 1.0!')
          fatal = .true.
       end if
    end if

    ! Damper parameters.
    do i = 1,ndampers
       if (Damper_Parameters(i) /= NULL_R) then
          if (Damper_Parameters(i) <= 0) then
             write (message, 12) i
12           format ('Damper parameter ',i1,' must be > 0.0!')
             call TLS_error (message)
             fatal = .true.
          end if
       end if   
    end do

  END SUBROUTINE NONLINEAR_SOLVER_CHECK

  SUBROUTINE NONLINEAR_SOLVER_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default NONLINEAR_SOLVER namelist variables to nonphysical values
    !   so that it is easy to tell if they have been input into the
    !   namelist.
    !
    !=======================================================================

    convergence_criterion  = NULL_R
    maximum_iterations     = NULL_I
    method                 = NULL_C
    name                   = NULL_C
    linear_solver_name     = NULL_C
    output_mode            = NULL_C
    perturbation_parameter = NULL_R
    Damper_Parameters      = NULL_R
    NLK_Max_Vectors        = NULL_I
    NLK_Vector_Tolerance   = NULL_R
    use_damper             = .true.

  END SUBROUTINE NONLINEAR_SOLVER_DEFAULT

  SUBROUTINE NONLINEAR_SOLVER_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast all elements of LINEAR_SOLVER namelist.
    !
    !======================================================================
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_BCAST

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.NOT. p_info%UseGlobalServices) then
       call PGSLib_BCAST (convergence_criterion)
       call PGSLib_BCAST (maximum_iterations)
       call PGSLib_BCAST (method)
       call PGSLib_BCAST (name)
       call PGSLib_BCAST (linear_solver_name)
       call PGSLib_BCAST (use_damper)
       call PGSLib_BCAST (output_mode)
       call PGSLib_BCAST (perturbation_parameter)
       call PGSLib_BCAST (Damper_Parameters)
       call PGSLib_BCAST (NLK_Max_Vectors)
       call PGSLib_BCAST (NLK_Vector_Tolerance)
    end if

  END SUBROUTINE NONLINEAR_SOLVER_INPUT_PARALLEL

  SUBROUTINE SET_NK_CONTROLS (NK, solution)
    !=======================================================================
    ! Purpose:
    !
    !   Initialize parameters and arrays necessary for UbikSolve.
    !
    !   We select a default solver here, but always defer to the users wishes
    !=======================================================================
    use linear_solution,    only: UBIK_NK_DEFAULT
    use nonlinear_solution, only: NK_CONTROL, NK_DEFAULT, ndampers

    ! Arguments
    type(NK_CONTROL), intent(INOUT) :: NK
    integer, intent(IN) :: solution
 
    ! Input variable defaults.
    character(string_len), parameter :: METHOD_DEFAULT             = 'nk'
    character(string_len), parameter :: NAME_DEFAULT               = 'default'
    character(string_len), parameter :: LINEAR_SOLVER_NAME_DEFAULT = 'default'
    character(string_len), parameter :: OUTPUT_MODE_DEFAULT        = 'none'
    integer,  parameter :: LINEAR_SOLVER_INDEX_DEFAULT    = UBIK_NK_DEFAULT
    integer,  parameter :: MAXIMUM_ITERATIONS_DEFAULT     = 30
    real(r8), parameter :: CONVERGENCE_CRITERION_DEFAULT  = 1.0d-5
    real(r8), parameter :: PERTURBATION_PARAMETER_DEFAULT = 1.0d-6
    real(r8), dimension(ndampers), parameter :: DAMPER_PARAMETERS_DEFAULT  = 1
    real(r8), parameter :: NLK_Vector_Tolerance_DEFAULT = 0.01d0
    integer,  parameter :: NLK_Max_Vectors_DEFAULT      = 20
    integer :: i, linear_solver_index
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Set input variable defaults corresponding
    ! to default linear solution controls.
    select case (solution)

       case (NK_DEFAULT)

          name                   = NAME_DEFAULT
          linear_solver_name     = LINEAR_SOLVER_NAME_DEFAULT
          method                 = METHOD_DEFAULT
          output_mode            = OUTPUT_MODE_DEFAULT
          convergence_criterion  = CONVERGENCE_CRITERION_DEFAULT
          maximum_iterations     = MAXIMUM_ITERATIONS_DEFAULT
          perturbation_parameter = PERTURBATION_PARAMETER_DEFAULT
          Damper_Parameters      = DAMPER_PARAMETERS_DEFAULT
          NLK_Vector_Tolerance   = NLK_Vector_Tolerance_DEFAULT
          NLK_Max_Vectors        = NLK_Max_Vectors_DEFAULT
          linear_solver_index    = LINEAR_SOLVER_INDEX_DEFAULT

       case DEFAULT

          ! If variables haven't been provided, set defaults.
          if (name == NULL_C)                   name                   = NAME_DEFAULT
          if (linear_solver_name == NULL_C)     linear_solver_name     = LINEAR_SOLVER_NAME_DEFAULT
          if (method == NULL_C)                 method                 = METHOD_DEFAULT
          if (output_mode == NULL_C)            output_mode            = OUTPUT_MODE_DEFAULT
          if (convergence_criterion == NULL_R)  convergence_criterion  = CONVERGENCE_CRITERION_DEFAULT
          if (maximum_iterations == NULL_I)     maximum_iterations     = MAXIMUM_ITERATIONS_DEFAULT
          if (perturbation_parameter == NULL_R) perturbation_parameter = PERTURBATION_PARAMETER_DEFAULT
          if (ALL(Damper_Parameters == NULL_R)) Damper_Parameters      = DAMPER_PARAMETERS_DEFAULT
          if (NLK_Vector_Tolerance == NULL_R)   NLK_Vector_Tolerance   = NLK_Vector_Tolerance_DEFAULT
          if (NLK_Max_Vectors == NULL_I)        NLK_Max_Vectors        = NLK_Max_Vectors_DEFAULT
          linear_solver_index = LINEAR_SOLVER_INDEX_DEFAULT

    end select

    ! Assign appropriate input variables for this nonlinear
    ! solver into components of the NK derived type.

    ! Nonlinear solution name.
    NK%name = name

    ! Nonlinear solution method.
    NK%method = method

    ! Linear solution name. 
    NK%linear_solver_name = linear_solver_name

    ! Linear solution index.
    NK%linear_solver_index = linear_solver_index

    ! Convergence criterion.
    NK%tolnewt = convergence_criterion

    ! Perturbation parameter.
    NK%eps_nk = perturbation_parameter

    ! Maximum iterations.
    NK%newton_itmax = maximum_iterations

    ! Damper flag.
    NK%use_damper = use_damper

    ! Damper parameters.
    NK%limit_low  = Damper_Parameters(1)
    NK%limit_high = Damper_Parameters(2)

    ! NLK parameters
    NK%NLK_Vector_Tolerance = NLK_Vector_Tolerance
    NK%NLK_Max_Vectors      = NLK_Max_Vectors

    ! Output mode
    NK%output_mode = output_mode

    ! Non-input initializations; allocate arrays and set various quantities to zero.
    ALLOCATE (NK%L2(0:maximum_iterations))
    ALLOCATE (NK%LI(0:maximum_iterations))
    ALLOCATE (NK%LI_Location(0:maximum_iterations))

    NK%linear_tot  = 0
    NK%newton_tot = 0
    do i = 0,maximum_iterations
       NK%L2(i)          = 0
       NK%LI(i)          = 0
       NK%LI_Location(i) = 0
    end do    

  END SUBROUTINE SET_NK_CONTROLS

END MODULE NONLIN_SOLVER_INPUT
