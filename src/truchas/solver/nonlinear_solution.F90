#include "f90_assert.fpp"

MODULE NONLINEAR_SOLUTION
  !=======================================================================
  !
  ! Nonlinear Solver Module
  !
  !=======================================================================
  use kind_module,      only: int_kind, real_kind, log_kind
  use parameter_module, only: string_len
  use truchas_logging_services

  implicit none

  private

  ! Public procedures.
  public :: NONLINEAR_SOLVE

  public :: NK_GET_SOLUTION_FIELD, NK_SET_SOLUTION_FIELD, &
            NK_INITIALIZE, NK_FINALIZE

  ! Public data and types.
  public :: NK_CONTROL, NK_FIELD, NK_STATE, NK_SOLUTION_FIELD, NKuser     

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
  ! Number of NKuser elements needed for default control parameters.
  integer (int_kind), parameter, public :: DEFAULT_NK_CONTROLS = 1
  integer (int_kind), parameter, public :: NK_DEFAULT = 1

  ! Max number of NK damper parameters
  integer (int_kind), parameter, public :: ndampers = 4

  integer (int_kind), parameter :: MAX_FIELDS = 10
  integer (int_kind), parameter :: MAX_TIMES  = 3

  type NK_CONTROL

     character (string_len) :: name
     character (string_len) :: method
     character (string_len) :: linear_solver_name
     character (string_len) :: output_mode

     real (real_kind) :: tolnewt
     real (real_kind) :: eps_NK
     real (real_kind) :: limit_high
     real (real_kind) :: limit_low

     real (real_kind)   :: NLK_Vector_Tolerance
     integer (int_kind) :: NLK_Max_Vectors

     real (real_kind),   pointer, dimension(:) :: L2
     real (real_kind),   pointer, dimension(:) :: LI
     integer (int_kind), pointer, dimension(:) :: LI_Location

     integer (int_kind) :: linear_tot
     integer (int_kind) :: newton_tot
     integer (int_kind) :: newton_itmax
     integer (int_kind) :: linear_solver_index
   
     logical (log_kind) :: use_damper
   
  end type NK_CONTROL

  real (real_kind), public, save, pointer, dimension(:) :: P_Residual
  real (real_kind), public, save, pointer, dimension(:) :: P_Past
  real (real_kind), public, save, pointer, dimension(:) :: P_Present
  real (real_kind), public, save, pointer, dimension(:) :: P_Future

  type (NK_CONTROL), public, save, pointer :: P_Control

  type NK_FIELD
     real (real_kind), pointer, dimension(:) :: Field
  end type

  type NK_STATE
     type (NK_FIELD), dimension(MAX_FIELDS) :: CC_Field
  end type

  type NK_SOLUTION_FIELD
     integer (int_kind)                      :: vectorsize
     integer (int_kind)                      :: unknowns_per_element
     integer (int_kind)                      :: time_levels
     type (NK_STATE),  dimension(MAX_TIMES)  :: t
     type (NK_STATE)                         :: residual
     real (real_kind), dimension(MAX_FIELDS) :: cc_norms
     real (real_kind)                        :: tot_norm
  end type

  ! Declare an array of the NK_CONTROL type to hold
  ! user-input NK nonlinear solution parameters
  type (NK_CONTROL), pointer, dimension(:) :: NKuser
 
  ! Number of user-specified NK nonlinear solution control parameters.
  integer (int_kind), public, save :: nonlinear_solutions

CONTAINS

  
  SUBROUTINE NONLINEAR_SOLVE (NK_Data, RESIDUAL, MATVEC, PRECONDITIONER, &
                    PRECONDITIONER_UPDATE, NLS, LS, status)
    !======================================================================
    ! Purpose:
    !
    !    Solution of Nonlinear System of Equations
    !
    !======================================================================
    use linear_solution,       only: Ubik_type
    use UbikSolve_module

    implicit none

    ! Arguments
    type (NK_SOLUTION_FIELD), intent(INOUT) :: NK_Data

    ! Interface blocks for external routines:
    !   RESIDUAL
    !   MATVEC
    !   PRECONDITIONER
    !   PRECONDITIONER_UPDATE
#include "solver_function_prototypes.fpp"

    optional :: PRECONDITIONER_UPDATE

    type (NK_CONTROL),  target, intent(INOUT) :: NLS
    type (Ubik_type),           intent(INOUT) :: LS
    integer (int_kind),         intent(INOUT) :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    select case (NLS%method)

    case ("nk","default")
       call NK_BB (NK_Data, RESIDUAL, MATVEC, PRECONDITIONER, &
            PRECONDITIONER_UPDATE, NLS, LS, status)
    case ("nlk")
       call nlk (NK_Data, RESIDUAL, MATVEC, PRECONDITIONER, &
            PRECONDITIONER_UPDATE, NLS, LS, status)
    end select

    return

  END SUBROUTINE NONLINEAR_SOLVE

  SUBROUTINE NK_BB (NK_Data, RESIDUAL, MATVEC, PRECONDITIONER, &
                    PRECONDITIONER_UPDATE, NK, LS, status)
    !======================================================================
    ! Purpose:
    !
    !    Newton-Krylov "Black-Box" nonlinear PDE solver driver.
    !
    !======================================================================
    use constants_module,      only: zero
    use kind_module,           only: real_kind, int_kind
    use linear_solution,       only: LINEAR_SOLVER, Ubik_type, COMING_FROM_INSIDE_NK
    use lnorm_module,          only: L2NORM, LINORM
    use pgslib_module,         only: PGSLIB_GLOBAL_MAXLOC, PGSLIB_GLOBAL_MAXVAL,    &
                                     PGSLib_Global_All
    use timing_tree
    use UbikSolve_module
#ifdef USE_TBROOK
    use output_data_module, only: enable_tbrook_output
#endif

    implicit none

    ! Arguments
    type (NK_SOLUTION_FIELD), intent(INOUT) :: NK_Data

    ! Interface blocks for external routines:
    !   RESIDUAL
    !   MATVEC
    !   PRECONDITIONER
    !   PRECONDITIONER_UPDATE
#include "solver_function_prototypes.fpp"

    optional :: PRECONDITIONER_UPDATE

    type (NK_CONTROL),  target, intent(INOUT) :: NK
    type (Ubik_type),           intent(INOUT) :: LS
    integer (int_kind),         intent(INOUT) :: status

    ! Local Variables
    logical (log_kind) :: TINY_INITIAL_RESIDUAL, fixed_point
    integer (int_kind) :: i, unknowns, OLD, CURRENT, NEW
    real (real_kind)   :: avg, alpha
    real (real_kind)   :: max_delta_norm, max_delta_norm_old, future_state_norm
    real (real_kind)   :: relative_max_delta_norm, convergence_rate, criterion
    integer (int_kind), dimension(1)        :: Location
    real (real_kind), dimension(:), pointer :: Solution_Delta, Solution_Res, RHS 
    character(256) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the NK timer.
    call start_timer ("NK")
    
    fixed_point = .false.

    ! Get the system size and set time level information.
    unknowns = NK_Data%vectorsize*NK_Data%unknowns_per_element
    select case (NK_Data%time_levels)
    case (2) ! First order in time
       OLD     = 1
       CURRENT = 1
       NEW     = 2
    case (3) ! Second order in time
       OLD     = 1
       CURRENT = 2
       NEW     = 3
    case default
       call TLS_panic ('NK_BB: invalid number of time levels specified in ' // trim(NK%name))
    end select

    ! Allocate working vectors.
    ALLOCATE (Solution_Delta(unknowns), STAT = status)
    call TLS_fatal_if_any (status /= 0, 'NK_BB: error allocating array Solution_Delta in ' // trim(NK%name))
    ALLOCATE (Solution_Res(unknowns), STAT = status)
    call TLS_fatal_if_any (status /= 0, 'NK_BB: error allocating array Solution_Res in ' // trim(NK%name))
    ALLOCATE (RHS(unknowns), STAT = status)
    call TLS_fatal_if_any (status /= 0, 'NK_BB: error allocating array RHS in ' // trim(NK%name))

    ! Point to the NK data supplied.
    P_Control  => NK
    P_Past     => NK_Data%t(OLD)%cc_field(1)%field
    P_Present  => NK_Data%t(CURRENT)%cc_field(1)%field
    P_Future   => NK_Data%t(NEW)%cc_field(1)%field
    P_Residual => Solution_Res
    NK_Data%Residual%cc_field(1)%field => Solution_Res

    ! Initialize various arrays and values.
    do i = 0,NK%newton_itmax
       NK%L2(i) = zero
       NK%LI(i) = zero
       NK%LI_Location(i) = 0
    end do    
    NK%newton_tot = 0
    NK%linear_tot = 0

    ! Get the initial residual.
    call start_timer ("NK Residual")
    call RESIDUAL (P_Past, P_Future, P_Residual)
    call stop_timer("NK Residual")

    ! Compute norms of the initial residual.
    NK%L2(0) = L2NORM(P_Residual)
    NK%LI(0) = LINORM(P_Residual)
    Location = PGSLIB_GLOBAL_MAXLOC(ABS(P_Residual))
    NK%LI_Location(0) = Location(1)

    ! check to see if initial residual is tiny - if so, we determine
    ! convergence by change in solution rather than reduction in
    ! nonlinear residual
    !
    ! note that TINY(1.0) is typically of the order 1e-38
    if (NK%L2(0) < TINY(1.0)) then
       TINY_INITIAL_RESIDUAL = .true.
    else
       TINY_INITIAL_RESIDUAL = .false.
    end if

    if (TLS_verbosity >= TLS_VERB_NOISY) then
       call TLS_info ('')
       write(message,'(30x,a)') 'Newton-Krylov Performance Diagnostics'
       call TLS_info (message)
       call TLS_info ('')
       call TLS_info ('     Iter Linear Iters L2(Residual) Linf(Residual) &
          &Location  Alpha    Linf(dX)  Location  Conv. Rate  Conv. Crit.')
       call TLS_info ('     ---- ------------ ------------ -------------- &
          &--------  -----    --------  --------  ----------  -----------')
       write(message,'(6x,i2,5x,i3,8x,1pe9.2,5x,1pe9.2,3x,i6)') 0, 0, NK%L2(0), NK%LI(0), NK%LI_Location(0)
       call TLS_info (message)
    end if

    ! Start the Newton loop counter.
    call start_timer ("Newton Loop")

    max_delta_norm_old = 1.0

    ! Newton loop.
    NEWTON_LOOP: do i = 1, NK%newton_itmax

      ! Call the linear solver; solution is the incremental change in the answer.
      Solution_Delta = -P_Residual            ! Initial guess
      RHS            = -P_Residual            ! RHS
      LS%status = COMING_FROM_INSIDE_NK
      
      call start_timer ("NK LS")
      call LINEAR_SOLVER (Solution_Delta, RHS, LS, MATVEC, PRECONDITIONER)
      call stop_timer ("NK LS")

      ! Initialize damping coefficient to unity.
      call start_timer ("NK Damper")
      call NK_DAMPER (P_Future, P_Present, Solution_Delta, unknowns, NK, alpha)
      call stop_timer ("NK Damper")

      !! A quick fix -- this needs to be redone properly (NNC)
      if (all(P_Future == P_Future + Solution_Delta)) fixed_point = .true.

      ! Add damped update.
      P_Future = P_Future + alpha*Solution_Delta
     
      ! Evaluate the new nonlinear residual.
      call start_timer ("NK Residual")
      call RESIDUAL (P_Past, P_Future, P_Residual)
      call stop_timer ("NK Residual")

      ! Test for nonlinear convergence via the L2 and L-infinity norm of the residual
      ! vector and the L-infinity norm of the ratio of the update to the solution.
      NK%L2(i) = L2NORM(P_Residual)
      NK%LI(i) = LINORM(P_Residual)
      Location = PGSLIB_GLOBAL_MAXLOC(ABS(P_Residual))
      NK%LI_Location(i) = Location(1)

      ! Also monitor relative change in the solution as a barometer on progress.
      future_state_norm = LINORM(P_Future)
      max_delta_norm = LINORM(Solution_Delta)
      if (future_state_norm > 0.0d0) then
        relative_max_delta_norm = max_delta_norm / future_state_norm
      else
        relative_max_delta_norm = 1.0d0
      end if



      ! Use convergence rate to modify user-specified tolerance and compute a value
      ! that will actually be used to determine convergence.
      !
      ! The idea is that convergence rate, computed as:
      !
      !         || x(k+1) - x(k) ||_inf
      !         -----------------------
      !             || x(k+1) ||_inf
      !     --------------------------------
      !         || x(k) - x(k-1) ||_inf
      !         -----------------------
      !              || x(k) ||_inf
      !
      ! is an estimate of the spectral radius
      if ((i.gt.1).and.(max_delta_norm_old .gt. 0.0)) then
         convergence_rate = max_delta_norm / max_delta_norm_old
      else
         ! Convergence rate not really defined on the first iteration
         convergence_rate = 0.d0
      endif

      criterion = (1.0 - convergence_rate) * NK%tolnewt

      if (TLS_verbosity >= TLS_VERB_NOISY) then
         Location = PGSLIB_GLOBAL_MAXLOC(ABS(Solution_Delta))
         write (message,10) i, Ubik_iter(LS%Control), NK%L2(i), NK%LI(i), NK%LI_Location(i), &
              alpha, PGSLIB_GLOBAL_MAXVAL(ABS(Solution_Delta)), Location, convergence_rate, criterion
10       format(6x,i2,5x,i3,8x,1pe9.2,5x,1pe9.2,3x,i6,3x,1pe9.2,1x,1pe9.2,2x,i6,4x,1pe9.2,4x,1pe9.2)
         call TLS_info (message)
      end if

      NK%newton_tot = NK%newton_tot + 1
      NK%linear_tot = NK%linear_tot + Ubik_iter(LS%Control)

      ! test for convergence of the nonlinear iteration
      !
      ! if either:
      !
      ! o relative change in the max-norm of the solution meets criterion
      ! o reduction in 2-norm of the nonlinear residual meets criterion
      !
      ! then we declare victory
    if (relative_max_delta_norm < criterion) then
         exit NEWTON_LOOP
      else if (.not. TINY_INITIAL_RESIDUAL) then
         if (NK%L2(i)/NK%L2(0) < criterion) then
            exit NEWTON_LOOP
         end if
      !! A quick fix -- this needs to be redone properly (NNC)
      else if (PGSLib_Global_All(fixed_point)) then
         exit NEWTON_LOOP
      end if

      ! Update the nonlinear preconditioner if called for.
      if (PRESENT(PRECONDITIONER_UPDATE) .and. LS%precond /= 0) then
         call PRECONDITIONER_UPDATE (P_Future)
      end if
      max_delta_norm_old = max_delta_norm

    end do NEWTON_LOOP

    ! Stop the Newton loop counter.
    call stop_timer ("Newton Loop")

    ! Check for convergence failure; if so, write out residual and normalized
    ! delta solution and punt.
    if (i >= NK%newton_itmax + 1) then
#ifdef USE_TBROOK
       if (enable_tbrook_output) call xml_write_residual ('NK', p_residual, p_future, solution_delta)
#endif
       call TLS_fatal ('NK_BB: Newton-Krylov iteration limit! Examine .err file for details')
    end if
    if (TLS_verbosity >= TLS_VERB_NOISY) then
       avg = zero
       if (NK%Newton_tot > 0) then
          avg = REAL(NK%linear_tot)/REAL(NK%Newton_tot)
       end if
       call TLS_info ('')
       write (message,20) avg
       call TLS_info (message)
       call TLS_info ('')
20     format(15x,'Average Number of Linear Iterations per Nonlinear Iteration: ',1pe11.4)
    end if

    ! This continue statement is sometimes useful for setting breakpoints
    ! when debugging.
100 continue

    ! Deallocate and nullify working arrays.
    NULLIFY (P_Residual)
    NULLIFY (P_Past)
    NULLIFY (P_Present)
    NULLIFY (P_Future)
    NULLIFY (P_Control)
    DEALLOCATE (Solution_Delta)
    DEALLOCATE (Solution_Res)
    DEALLOCATE (RHS)
 
    ! Stop the NK timer.
    call stop_timer("NK")

    return

  END SUBROUTINE NK_BB

  SUBROUTINE NLK (NK_Data, RESIDUAL, MATVEC, PRECONDITIONER, &
                    PRECONDITIONER_UPDATE, NLS, LS, status)
    !======================================================================
    ! Purpose:
    !
    !    Approximate Inexact Newton Method for Solving Nonlinear Equations
    !
    !======================================================================
    use constants_module,      only: zero
    use debug_control_data,    only: verbose, VERBOSE_NOISY
    use fixed_point_accelerator, only: fpa_create, fpa_destroy, fpa_correction, fpa_state 
    use kind_module,           only: real_kind, int_kind
!    use linear_solution,       only: LINEAR_SOLVER, Ubik_type, COMING_FROM_INSIDE_NK
!    use linear_solution
    use linear_solution,       only: Ubik_type, Ubik_solver, PRECOND_NONE, &
                                     SOLVER_NONE
    use lnorm_module,          only: L2NORM, LINORM
    use pgslib_module,         only: PGSLIB_GLOBAL_MAXLOC, PGSLIB_GLOBAL_MAXVAL,    &
                                     PGSLib_Global_All
    use timing_tree
    use UbikSolve_module
#ifdef USE_TBROOK
    use output_data_module, only: enable_tbrook_output
#endif

    implicit none

    ! Arguments
    type (NK_SOLUTION_FIELD), intent(INOUT) :: NK_Data

    ! Interface blocks for external routines:
    !   RESIDUAL
    !   MATVEC
    !   PRECONDITIONER
    !   PRECONDITIONER_UPDATE
#include "solver_function_prototypes.fpp"

    optional :: PRECONDITIONER_UPDATE

    type (NK_CONTROL),  target, intent(INOUT) :: NLS
    type (Ubik_type),           intent(INOUT) :: LS
    type (Ubik_vector_type)                   :: ubik_vec
    integer (int_kind),         intent(INOUT) :: status

    ! Local Variables
    logical (log_kind) :: TINY_INITIAL_RESIDUAL, fixed_point
    integer (int_kind) :: i, unknowns, OLD, CURRENT, NEW, n
    real (real_kind)   :: avg, alpha
    real (real_kind)   :: max_delta_norm, max_delta_norm_old, future_state_norm
    real (real_kind)   :: relative_max_delta_norm, convergence_rate, criterion
    integer (int_kind), dimension(1)        :: Location
    real (real_kind), dimension(:), pointer :: Solution_Delta, Solution_Res, RHS 
    type (fpa_state)   :: accel_state
    integer(int_kind) :: mvec
    character(256) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the NK timer.
    call start_timer ("NK")

    fixed_point = .false.

    ! Get the system size and set time level information.
    unknowns = NK_Data%vectorsize*NK_Data%unknowns_per_element
    select case (NK_Data%time_levels)
    case (2) ! First order in time
       OLD     = 1
       CURRENT = 1
       NEW     = 2
    case (3) ! Second order in time
       OLD     = 1
       CURRENT = 2
       NEW     = 3
    case default
       call TLS_panic ('NLK: invalid number of time levels specified in ' // trim(NLS%name))
    end select

    ! Allocate working vectors.
    ALLOCATE (Solution_Delta(unknowns), STAT = status)
    call TLS_fatal_if_any (status /= 0, 'NLK: error allocating array Solution_Delta in ' // trim(NLS%name))
    ALLOCATE (Solution_Res(unknowns), STAT = status)
    call TLS_fatal_if_any (status /= 0, 'NLK: error allocating array Solution_Res in ' // trim(NLS%name))
    ALLOCATE (RHS(unknowns), STAT = status)
    call TLS_fatal_if_any (status /= 0, 'NLK: error allocating array RHS in ' // trim(NLS%name))

    ! Point to the NK data supplied.
    P_Control  => NLS
    P_Past     => NK_Data%t(OLD)%cc_field(1)%field
    P_Present  => NK_Data%t(CURRENT)%cc_field(1)%field
    P_Future   => NK_Data%t(NEW)%cc_field(1)%field
    P_Residual => Solution_Res
    NK_Data%Residual%cc_field(1)%field => Solution_Res

    ! Initialize various arrays and values.
    do i = 0,NLS%newton_itmax
       NLS%L2(i) = zero
       NLS%LI(i) = zero
       NLS%LI_Location(i) = 0
    end do    
    NLS%newton_tot = 0
    NLS%linear_tot = 0

    ! Get the initial residual.
    call start_timer ("NK Residual")
    call RESIDUAL (P_Past, P_Future, P_Residual)
    call stop_timer ("NK Residual")

    ! Compute norms of the initial residual.
    NLS%L2(0) = L2NORM(P_Residual)
    NLS%LI(0) = LINORM(P_Residual)
    Location = PGSLIB_GLOBAL_MAXLOC(ABS(P_Residual))
    NLS%LI_Location(0) = Location(1)

    ! check to see if initial residual is tiny - if so, we determine
    ! convergence by change in solution rather than reduction in
    ! nonlinear residual
    !
    ! note that TINY(1.0) is typically of the order 1e-38
    if (NLS%L2(0) < TINY(1.0)) then
       TINY_INITIAL_RESIDUAL = .true.
    else
       TINY_INITIAL_RESIDUAL = .false.
    end if

    if (TLS_verbosity >= TLS_VERB_NOISY) then
       call TLS_info ('')
       write(message,'(30x,a)') 'Accelerated Inexact Newton Performance Diagnostics'
       call TLS_info (message)
       call TLS_info ('')
       call TLS_info ('     Iter Linear Iters L2(Residual) Linf(Residual) &
          &Location  Alpha    Linf(dX)  Location  Conv. Rate  Conv. Crit.')
       call TLS_info ('     ---- ------------ ------------ -------------- &
          &--------  -----    --------  --------  ----------  -----------')
       write(message,'(6x,i2,5x,i3,8x,1pe9.2,5x,1pe9.2,3x,i6)') 0, 0, NLS%L2(0), NLS%LI(0), NLS%LI_Location(0)
       call TLS_info (message)
    end if

    ! Start the Newton loop counter.
    call start_timer ("Newton Loop")

    mvec=min(NLS%NLK_Max_Vectors,NLS%newton_itmax)
    ubik_vec%values => Solution_Delta

    ! Initialize fixed point accelerator
    call fpa_create(accel_state,size(Solution_Delta),mvec,NLS%NLK_Vector_Tolerance)

    if (LS%solver.ne.SOLVER_NONE) then
       LS%solver=SOLVER_NONE
       call TLS_warn ('LINEAR_SOLVER METHOD changed to NONE when NONLINEAR_SOLVER METHOD is NLK')
    endif

    ! Copy and allocate the linear solver control structure (Ubik).
    Ubik_solver = LS
    call Ubik_create (Ubik_solver%control)

    ! setup the preconditioning options
    PRECONDITIONER_SETUP: if (Ubik_solver%precond /= PRECOND_NONE) then
 
       Ubik_solver%precond_iter = 0
       Ubik_solver%Factor = .true.

    end if PRECONDITIONER_SETUP

    ! Newton loop.
    NEWTON_LOOP: do i = 1, NLS%newton_itmax

      ! Solution is the incremental change in the answer.
 
       Solution_Delta = 0.d0  ! (using -P_Residual as initial guess does not work)
       RHS            = -P_Residual  

! This timer effectively records the amount of time used in doing a SOLVE
! with the preconditioner.
      call start_timer ("NK LS")

! Do not make the following following call to TIMER_START, because
! it is already made inside of the PRECONDITIONER call, and it 
! would actually cause double-counting.
!!$      call TIMER_START (TIMER_NK_LS_PRECONDITIONER)

      call PRECONDITIONER(RHS,ubik_vec,status)

!!$      call TIMER_STOP (TIMER_NK_LS_PRECONDITIONER)
      call stop_timer ("NK LS")
      

      ! The Newton accelerator of Carlson/Miller SIAM Sci. Comp. Vol 19 pp 728-765 1998
      call fpa_correction(accel_state,Solution_Delta,dp=nlk_dot_product)
     
! FPA is not currently set up to handle damping.
      alpha=1.d0

      !! A quick fix -- this needs to be redone properly (NNC)
      if (all(P_Future == P_Future + Solution_Delta)) fixed_point = .true.

      P_Future = P_Future + Solution_Delta
     
      ! Evaluate the new nonlinear residual.
      call start_timer ("NK Residual")
      call RESIDUAL (P_Past, P_Future, P_Residual)
      call stop_timer ("NK Residual")
      
      ! Monitor the L2 and L-infinity norm of the residual vector.
      NLS%L2(i) = L2NORM(P_Residual)
      NLS%LI(i) = LINORM(P_Residual)
      Location = PGSLIB_GLOBAL_MAXLOC(ABS(P_Residual))
      NLS%LI_Location(i) = Location(1)

      ! Also monitor relative change in the solution.
      future_state_norm = LINORM(P_Future)
      max_delta_norm = LINORM(Solution_Delta)
      !! A quick fix -- this needs to be redone properly (NNC)
      if (future_state_norm > 0.0d0) then
        relative_max_delta_norm = max_delta_norm / future_state_norm
      else
        relative_max_delta_norm = 1.0d0
      end if

      ! Use convergence rate to modify user-specified tolerance and compute a value
      ! that will actually be used to determine convergence.
      !
      ! The idea is that convergence rate, computed as:
      !
      !         || x(k+1) - x(k) ||_inf
      !     --------------------------------
      !         || x(k) - x(k-1) ||_inf
      !
      ! is an estimate of the spectral radius

      if (i.gt.1) then
         convergence_rate = max_delta_norm / max_delta_norm_old
      else
         ! Convergence rate not really defined on the first iteration
         convergence_rate = 0.d0
      endif

!     If we assume linear convergence (which is pessimistic since Newton's method
!     is quadratically convergent and the secant method is superlinearly convergent),
!     then the error in the last iterate would be roughly
!         convergence_rate/(1-convergence_rate) * change in last iterate
!     Thus we should set
!         criterion = (1-convergence_rate)/convergence_rate * NLS%tolnewt
!     If convergence_rate < 1/2, then this would give a criterion larger
!     than NLS%tolnewt.
!     We take the conservative approach and set
!         criterion = (1-convergence_rate)/convergence_rate * NLS%tolnewt  , if convergence_rate >= 1/2
!         criterion = NLS%tolnewt                                          , if convergence_rate < 1/2

      if (convergence_rate.ge.0.5d0) then
        criterion = (1.d0-convergence_rate)/convergence_rate * NLS%tolnewt
      else
        criterion = NLS%tolnewt
      endif

      if (TLS_verbosity >= TLS_VERB_NOISY) then
         Location = PGSLIB_GLOBAL_MAXLOC(ABS(Solution_Delta))
         write (message,10) i, Ubik_iter(LS%Control), NLS%L2(i), NLS%LI(i), NLS%LI_Location(i), &
              alpha, PGSLIB_GLOBAL_MAXVAL(ABS(Solution_Delta)), Location, convergence_rate, criterion
10       format(6x,i2,5x,i3,8x,1pe9.2,5x,1pe9.2,3x,i6,3x,1pe9.2,1x,1pe9.2,2x,i6,4x,1pe9.2,4x,1pe9.2)
         call TLS_info (message)
      end if

      NLS%newton_tot = NLS%newton_tot + 1
      NLS%linear_tot = NLS%linear_tot + Ubik_iter(LS%Control)

!  At least do two iterations...
      if (i.gt.1) then
! If relative change in the max-norm of the solution meets criterion
! then we declare victory
         if (relative_max_delta_norm < criterion) then
            exit NEWTON_LOOP
!  We take a conservative approach and do not kick out just
!  because the residual has been greatly reduced.
         else if (.not. TINY_INITIAL_RESIDUAL) then
            if (NLS%L2(i)/NLS%L2(0) < criterion) then
               exit NEWTON_LOOP
            end if
         end if
      !! A quick fix -- this needs to be redone properly (NNC)
      else if (PGSLib_Global_All(fixed_point)) then
         exit NEWTON_LOOP
      endif

! We currently do not allow updating of the preconditioner, since the
! usual presentation of the FPA has the preconditioner being 
! an unchanging linear operator.
!!$      ! Update the nonlinear preconditioner if called for.
!!$      if (PRESENT(PRECONDITIONER_UPDATE) .and. LS%precond /= 0) then
!!$         call PRECONDITIONER_UPDATE (P_Future)
!!$      end if

      max_delta_norm_old = max_delta_norm

    end do NEWTON_LOOP

    ! Check for errors in solver; print diagnostics and quit
!    call POST_SOLVE (Ubik_solver, 'LINEAR_SOLVER')

    ! Deallocate acceleration vectors
    call fpa_destroy(accel_state)

    ! copy Ubik_solver out to Ubik before deallocating
    LS%status = Ubik_solver%status
    call Ubik_set_iter (LS%control, Ubik_iter(Ubik_solver%control))
    LS%precond_iter = Ubik_solver%precond_iter
    call Ubik_destroy (Ubik_solver%control)

    ! Stop the Newton loop counter. 
    call stop_timer ("Newton Loop")

    ! Check for convergence failure; if so, write out residual and punt.
    if (i >= NLS%newton_itmax + 1) then
#ifdef USE_TBROOK
       if (enable_tbrook_output) call xml_write_residual ('nlk', p_residual, p_future, solution_delta)
#endif

! Print residual norms on failure
       call TLS_info ('')
       call TLS_info ('                    Iter  L2(Residual)  Linf(Residual)  Linf Location')
       call TLS_info ('                    ----  ------------  --------------  -------------')

       do n = 0,NLS%Newton_tot

          write (message, 18) n, NLS%L2(n), NLS%LI(n), &
                                    NLS%LI_Location(n)
18        format (21x,i4,3x,1pe11.4,4x,1pe11.4,6x,i6)
          call TLS_info (message)

       end do   
! End print residual norms

       call TLS_fatal ('NLK: Newton iteration limit! Examine .err file for details')
    end if
    if (TLS_verbosity >= TLS_VERB_NOISY) then
       avg = zero
       if (NLS%Newton_tot > 0) then
          avg = REAL(NLS%linear_tot)/REAL(NLS%Newton_tot)
       end if
       write (message,20) avg
20     format(15x,'Average Number of Linear Iterations per Nonlinear Iteration: ',1pe11.4)
       call TLS_info ('')
       call TLS_info (message)
       call TLS_info ('')
    end if

    ! This continue statement is sometimes useful for setting breakpoints
    ! when debugging.
100 continue

    ! Deallocate and nullify working arrays.
    NULLIFY (P_Residual)
    NULLIFY (P_Past)
    NULLIFY (P_Present)
    NULLIFY (P_Future)
    NULLIFY (P_Control)
    DEALLOCATE (Solution_Delta)
    DEALLOCATE (Solution_Res)
    DEALLOCATE (RHS)
 
    ! Stop the NK timer.
    call stop_timer ("NK")

    return

  END SUBROUTINE NLK

  !!
  !! We need to pass the dot product routine to FPA_CORRECTION, but it
  !! is illegal to pass a generic name like PGSLIB_GLOBAL_DOT_PRODUCT
  !! (or an internal procedure name for that matter).  Hence this routine.
  !!
  
  function nlk_dot_product (a, b) result (dp)
    use pgslib_module, only: pgslib_global_dot_product
    real(real_kind), intent(in) :: a(:), b(:)
    real(real_kind) :: dp
    dp = pgslib_global_dot_product(a, b)
  end function nlk_dot_product

  SUBROUTINE NK_DAMPER (Future_state, Present_state, Solution_Delta, unknowns, NK, alpha_min)
    !===========================================================================
    ! Purpose:
    !    NK Solution update damper. WARNING: This damper *only* works if
    !    the solution is positive definite everywhere.
    !
    !===========================================================================
    use constants_module,      only: one, one_tenth
    use kind_module,           only: real_kind, int_kind
    use pgslib_module,         only: PGSLIB_GLOBAL_MINVAL

    implicit none

    ! Arguments
    real (real_kind), dimension(:), intent(IN)    :: Future_state
    real (real_kind), dimension(:), intent(IN)    :: Present_state
    real (real_kind), dimension(:), intent(IN)    :: Solution_Delta
    integer (int_kind),             intent(IN)    :: unknowns
    type (NK_CONTROL),              intent(IN)    :: NK
    real (real_kind),               intent(INOUT) :: alpha_min

    ! Local Variables
    integer (int_kind) :: i
    real (real_kind)   :: alpha, factor, alpha_floor = one_tenth, &
                          alpha_ceiling = one

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the damping scaler (alpha) to unity.
    alpha = one

    if (NK%use_damper) then
       do i = 1,unknowns
          if (Future_state(i)+Solution_Delta(i) < NK%limit_low*Future_state(i)) then

             factor = -(one - NK%limit_low ) * Future_state(i) / Solution_Delta(i)
             if (factor < alpha) then
                alpha = factor
             end if

          else if (Future_state(i)+Solution_Delta(i) > NK%limit_high*Future_state(i)) then

             factor = (NK%limit_high - one ) * Future_state(i) / Solution_Delta(i)
             if (factor < alpha ) then
                alpha = factor
             end if

          end if
       end do
    end if

    ! Make sure alpha_floor <= alpha <= alpha_ceiling.
    alpha = MIN(alpha,alpha_ceiling)
    alpha = MAX(alpha,alpha_floor)

    ! Take the global minimum over the entire mesh.
    alpha_min = PGSLIB_GLOBAL_MINVAL(alpha)

    return

  END SUBROUTINE NK_DAMPER

  SUBROUTINE NK_GET_SOLUTION_FIELD (Solution_Field, Vector, time_level)
    !===========================================================================
    ! Purpose:
    !
    !   Retrieve a real rank-1 vector from the NK_Solution_Field data type. 
    !
    !===========================================================================
    use kind_module, only: int_kind, real_kind

    implicit none

    ! Arguments
    type (NK_SOLUTION_FIELD),       intent(IN)    :: Solution_Field
    real (real_kind), dimension(:), intent(INOUT) :: Vector
    integer (int_kind),             intent(IN)    :: time_level
  
    ! Local Variables
    integer (int_kind) :: j, lb, ub

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Put the cc quantities into the vector.
    do j = 1,Solution_Field%unknowns_per_element
      lb = (j-1)*Solution_Field%vectorsize + 1
      ub = j*Solution_Field%vectorsize
      Vector(lb:ub) = Solution_Field%t(time_level)%cc_field(j)%field(:)   
    end do
    
    return

  END SUBROUTINE NK_GET_SOLUTION_FIELD

  SUBROUTINE NK_SET_SOLUTION_FIELD (Vector, Solution_Field, time_level)
    !============================================================================
    ! Purpose:
    !
    !   Set the NK_Solution_Field data type by pointing
    !   it to a real rank-1 Vector.
    !
    !===========================================================================
    use kind_module, only: int_kind, real_kind

    implicit none
 
    ! Arguments
    real (real_kind), dimension(:), pointer :: Vector
    type (NK_SOLUTION_FIELD), intent(INOUT) :: Solution_Field
    integer (int_kind),          intent(IN) :: time_level
    
    ! Local Variables.
    integer (int_kind) :: j, lb, ub
    real (real_kind), dimension(:), pointer :: Tmp => NULL()

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Put the cc quantities into the vector
    do j = 1,Solution_Field%unknowns_per_element
      lb = (j-1)*Solution_Field%vectorsize + 1
      ub = j*Solution_Field%vectorsize
      Tmp => Vector(lb:ub)
      Solution_Field%t(time_level)%cc_field(j)%field => Tmp
      NULLIFY (Tmp)
    end do
    
    return

  END SUBROUTINE NK_SET_SOLUTION_FIELD
  
  SUBROUTINE NK_INITIALIZE (Solution_Field, vectorsize, unknowns_per_element, &
                            time_levels, Solution_Old, Solution_Current)
    !============================================================================
    ! Purpose:
    !
    !   Initialize the NK_Solution_Field data type.
    !
    !===========================================================================
    use kind_module,  only: int_kind, real_kind

    implicit none
 
    ! Arguments
    type (NK_SOLUTION_FIELD),       intent(INOUT) :: Solution_Field
    integer (int_kind),             intent(IN)    :: vectorsize
    integer (int_kind),             intent(IN)    :: unknowns_per_element
    integer (int_kind),             intent(IN)    :: time_levels
    real (real_kind), dimension(:), intent(IN)    :: Solution_Old
    real (real_kind), dimension(:), intent(IN)    :: Solution_Current
    
    ! Local Variables.
    integer (int_kind) :: vlength, j, lb, ub, OLD, CURRENT, NEW, status, n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Store the vectorsize, number of unknowns per element, and time levels.
    Solution_Field%vectorsize           = vectorsize
    Solution_Field%unknowns_per_element = unknowns_per_element
    Solution_Field%time_levels          = time_levels

    ! Get the vector length and make sure it is equal to the size of Solution
    vlength = Solution_Field%vectorsize * Solution_Field%unknowns_per_element
    if (vlength /= SIZE(Solution_Old)) then
       call TLS_panic ('NK_INITIALIZE: inconsistent vector length size')
    end if

    ! Check and set all time levels.
    select case (time_levels)
       case (2) ! First order in time
          OLD     = 1
          CURRENT = 1
          NEW     = 2
       case (3) ! Second order in time
          OLD     = 1
          CURRENT = 2
          NEW     = 3
       case default
          call TLS_panic ('NK_INITIALIZE: Invalid number of time levels specified')
    end select

    ! Allocate and initialize the solution field vectors.
    do n = 1, Solution_Field%time_levels
       do j = 1, Solution_Field%unknowns_per_element
          lb = (j-1)*Solution_Field%vectorsize + 1
          ub = j*Solution_Field%vectorsize
          ALLOCATE(Solution_Field%t(n)%cc_field(j)%field(vectorsize), STAT = status)
          call TLS_fatal_if_any (status /= 0, 'NK_INITIALIZE: Error allocating Solution_Field')
          if (n == OLD .or. n == CURRENT) then
             Solution_Field%t(n)%cc_field(j)%field = Solution_Old(lb:ub)
          else if (n == NEW) then
             Solution_Field%t(n)%cc_field(j)%field = Solution_Current(lb:ub)
          end if
       end do
    end do
    
    return

  END SUBROUTINE NK_INITIALIZE
  
  SUBROUTINE NK_FINALIZE (Solution_Field)
    !============================================================================
    ! Purpose:
    !
    !   Destroy and deallocate all data associated with Solution_Field
    !
    !===========================================================================
    use kind_module, only: int_kind

    implicit none
 
    ! Arguments
    type (NK_SOLUTION_FIELD), intent(INOUT) :: Solution_Field
    
    ! Local Variables.
    integer (int_kind) :: j, lb, ub, n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Allocate and initialize the solution field vectors.
    do n = 1, Solution_Field%time_levels
       do j = 1, Solution_Field%unknowns_per_element
          lb = (j-1)*Solution_Field%vectorsize + 1
          ub = j*Solution_Field%vectorsize
          DEALLOCATE(Solution_Field%t(n)%cc_field(j)%field)
       end do
    end do
    
    return

  END SUBROUTINE NK_FINALIZE

#ifdef USE_TBROOK
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! XML_WRITE_RESIDUAL
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 19 Jul 2005
 !!
 !! Write the passed residual vector to the XML output file.  This creates
 !! a 'NONLIN_RESIDUAL' xml element in the output file; the data itself is
 !! written to a binary lookaside file.
 !!
 !! NB: This is extremely similar to a procedure of the same name from the
 !! LINEAR_SOLUTION module; the two ought to be consolidated into a common
 !! code base somehow.
 !!

  subroutine xml_write_residual (name, r, x, d)

    use brook_module
    use tbrook_module
    use output_module, only: prefix
    use string_utilities, only: i_to_c
    use parameter_module, only: ncells, nnodes

    character(len=*), intent(in) :: name
    real(real_kind),  intent(in) :: r(:), x(:), d(:)

    integer :: status, dim
    integer, save :: df_num = 0
    character(len=256) :: df_name
    type(brook), target :: df_brook
    real(real_kind), allocatable :: tmp(:,:)

    status = 0 ! for some insane reason, this is intent(in) for all the tbrook stuff.

    df_num = df_num + 1
    df_name = trim(prefix) // '.nonlin_res.' // i_to_c(df_num) // '.bin'

    !! Create the binary look-aside file.
    call tbrook_set (df_brook, file=trim(df_name), form='binary', istatus=status)
    if (status /= 0) return

    !! Open the NONLIN_RESIDUAL tag.
    call tbrook_openxmltag (BaseBrook, XMLTag='NONLIN_RESIDUAL', &
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
      call tbu_make_file_entry (BaseBrook, df_brook, 'delta',    d, status, map='cell')
    else if (mesh_based_scalar_data(size(r),nnodes)) then ! scalar, node-based data
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', r, status, map='node')
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', x, status, map='node')
      call tbu_make_file_entry (BaseBrook, df_brook, 'delta',    d, status, map='node')
    else if (mesh_based_vector_data(size(r),ncells,dim)) then ! vector, cell-based data
      allocate(tmp(dim,ncells))
      call copy_to_rank_2 (r, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', tmp, status, map='cell')
      call copy_to_rank_2 (x, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', tmp, status, map='cell')
      call copy_to_rank_2 (d, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'delta',    tmp, status, map='cell')
      deallocate(tmp)
    else if (mesh_based_vector_data(size(r),nnodes,dim)) then ! vector, node-based data
      allocate(tmp(dim,nnodes))
      call copy_to_rank_2 (r, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', tmp, status, map='node')
      call copy_to_rank_2 (x, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', tmp, status, map='node')
      call copy_to_rank_2 (d, tmp)
      call tbu_make_file_entry (BaseBrook, df_brook, 'delta',    tmp, status, map='node')
      deallocate(tmp)
    else ! we have no clue what it is ...
      call tbu_make_file_entry (BaseBrook, df_brook, 'residual', r, status, map='none')
      call tbu_make_file_entry (BaseBrook, df_brook, 'solution', x, status, map='none')
      call tbu_make_file_entry (BaseBrook, df_brook, 'delta',    d, status, map='none')
    end if
    if (status /= 0) return

    !! Close the NONLIN_RESIDUAL tag.
    call tbrook_closexmltag (BaseBrook, XMLTag='NONLIN_RESIDUAL', istatus=status)
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
END MODULE NONLINEAR_SOLUTION
