!!
!! BDF2_DAE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 14 Aug 2006
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL CREATE_STATE (THIS, N, MITR, NTOL, MVEC, VTOL)
!!
!!  CALL SET_INITIAL_STATE (THIS, T, U, UDOT)
!!
!!  CALL BDF2_STEP (THIS, PCFUN, UPDPC, ENORM, H, HMIN, MTRY, U, HNEXT, ERRC)
!!
!!  CALL BDF2_STEP_SIMPLE (THIS, PCFUN, UPDPC, ENORM, H, U, HNEXT, ERRC, ECTRL)
!!
!!  CALL COMMIT_SOLUTION (THIS, H, U)
!!
!!  CALL SET_VERBOSE_STEPPING (THIS, UNIT)
!!  CALL QUIET_STEPPING (THIS)
!!
!!  CALL BDF2_STEP_DRIVER (THIS, PCFUN, UPDPC, ENORM, HNEXT, STATUS &
!!                     [, NSTEP] [, TOUT] [, HMIN] [, HMAX] [, MTRY])
!!
!!  LAST_TIME(THIS) returns the lastest time.
!!  LAST_SOLUTION(THIS) returns a pointer to latest solution vector.
!!  LAST_STEP_SIZE(THIS) returns the last (successful) step size used.
!!  INTERPOLATED_SOLUTION(THIS, T) returns the solution interpolated to time T
!!    from the recent history of solution vectors stored in the state THIS.
!!    This routine should really only be called when the last two step times
!!    straddle the time value T.
!!

#include "f90_assert.fpp"

module bdf2_dae

  use kinds
  use solution_history
  use nka_type
  implicit none
  private

  public :: create_state, set_initial_state, destroy_state, reset_state
  public :: bdf2_step, commit_solution, bdf2_step_driver
  public :: set_verbose_stepping, set_quiet_stepping
  public :: last_time, last_solution, last_step_size, interpolated_solution
  public :: write_bdf2_stepping_statistics, get_bdf2_stepping_statistics
  
  public :: bdf1_step, bdf1_step_new

  type, public :: state
    private
    integer  :: n                   ! number of unknowns
    integer  :: seq = -1            ! number of steps taken
    real(r8) :: hlast               ! last step size
    real(r8) :: hpc                 ! step size built into the current preconditioner
    logical  :: usable_pc = .false. ! whether the current preconditioner is usable
    integer  :: freeze_count        ! don't increase step size for this number of steps
    integer  :: mitr = 5            ! maximum number of nonlinear iterations
    real(r8) :: ntol = 0.1_r8       ! nonlinear solver error tolerance (relative to 1)
    type(nka) :: fpa                ! nonlinear solver (NKA) accelerator structure
    type(history)   :: uhist        ! solution history structure

    !! Perfomance counters
    integer :: pcfun_calls = 0      ! number of calls to PCFUN
    integer :: updpc_calls = 0      ! number of calls to UPDPC
    integer :: updpc_failed = 0     ! number of UPDPC calls returning an error
    integer :: retried_bce = 0      ! number of retried BCE steps
    integer :: failed_bce = 0       ! number of completely failed BCE steps
    integer :: rejected_steps = 0   ! number of steps rejected on error tolerance
    real(r8) :: hmin = huge(1.0_r8) ! minimum step size used on a successful step
    real(r8) :: hmax = tiny(1.0_r8) ! maximum step size used on a successful step

    !! Diagnostics
    integer :: unit = 0
    logical :: verbose = .false.
  end type state

  real(r8), parameter, private :: RMIN = 0.25_r8
  real(r8), parameter, private :: RMAX = 4.0_r8
  real(r8), parameter, private :: MARGIN = 3.0_r8

  !! Successful STATUS return codes:
  integer, parameter, public :: SOLVED_TO_TOUT = 1
  integer, parameter, public :: SOLVED_TO_NSTEP = 2

  !! Unsuccessful STATUS return codes:
  integer, parameter, public :: BAD_INPUT = -1
  integer, parameter, public :: STEP_FAILED = -2
  integer, parameter, public :: STEP_SIZE_TOO_SMALL = -3

contains

  subroutine create_state (this, n, mitr, ntol, mvec, vtol)

    type(state), intent(out) :: this
    integer,  intent(in) :: n
    integer,  intent(in) :: mitr, mvec
    real(r8), intent(in) :: ntol, vtol
    optional :: mitr, ntol, mvec, vtol

    integer :: maxv

    INSIST( n > 0 )
    this%n = n

    if (present(mitr)) then
      INSIST( mitr > 1 )
      this%mitr = mitr
    end if

    if (present(ntol)) then
      INSIST( ntol > 0.0_r8 .and. ntol <= 1.0_r8 )
      this%ntol = ntol
    end if

    maxv = this%mitr - 1
    if (present(mvec)) then
      INSIST( mvec > 0 )
      maxv = min(maxv, mvec)
    end if

    !! Initialize the NKA structure.
    call nka_init (this%fpa, this%n, maxv)
    call nka_set_vec_tol (this%fpa, vtol)

    !! We need to maintain 3 solution vectors for quadratic extrapolation.
    call create_history (this%uhist, 3, this%n)

  end subroutine create_state

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_INITIAL_STATE
 !!

  subroutine set_initial_state (this, t, u, udot)

    type(state), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), udot(:)

    ASSERT( size(u) == this%n )
    ASSERT( size(udot) == this%n )

    call flush_history (this%uhist, t, u, udot)
    this%seq  = 0

  end subroutine set_initial_state
  
  subroutine reset_state (this, index, u, udot)
  
    type(state), intent(inout) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: u
    real(r8), intent(in), optional :: udot
    
    ASSERT( index >= 1 .and. index <= this%n )
    
    call revise_history (this%uhist, index, u, udot)
    
    this%usable_pc = .false.
    
  end subroutine reset_state

  subroutine destroy_state (this)
    type(state), intent(inout) :: this
    call destroy (this%uhist)
    call nka_delete (this%fpa)
  end subroutine destroy_state

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  BDF2_STEP_DRIVER
 !!

  subroutine bdf2_step_driver (this, pcfun, updpc, enorm, hnext, status, &
                               nstep, tout, hmin, hmax, mtry)

    type(state), intent(inout) :: this
    real(r8), intent(inout) :: hnext
    integer,  intent(out)   :: status
    integer,  intent(in) :: nstep, mtry
    real(r8), intent(in) :: tout, hmin, hmax
    optional :: nstep, tout, hmin, hmax, mtry

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    integer  :: max_step, max_try, step
    real(r8) :: tlast, h, t_out, h_min, h_max, u(this%n)

    ASSERT( this%seq >= 0 )

   !!!
   !!! PROCESS THE INPUT AND CHECK IT FOR CORRECTNESS

    status = 0

    !! Set the maximum number of time steps; default is unlimited.
    max_step = huge(1)
    if (present(nstep)) max_step = nstep
    if (max_step < 1) status = BAD_INPUT

    !! Set the target integration time; default is +infinity.
    t_out = huge(1.0_r8)
    if (present(tout)) t_out = tout
    if (t_out <= most_recent_time(this%uhist)) status = BAD_INPUT

    !! Verify that at least one of NSTEP and TOUT were specified.
    if (.not.present(nstep) .and. .not.present(tout)) status = BAD_INPUT

    !! Set the bounds on the step size; default is none.
    h_min = tiny(1.0_r8)
    h_max = huge(1.0_r8)
    if (present(hmin)) h_min = hmin
    if (present(hmax)) h_max = hmax
    if (h_min < 0.0_r8 .or. hnext < h_min .or. hnext > h_max) status = BAD_INPUT

    !! Set the maximum number of attempts allowed for a step.
    max_try = 10
    if (present(mtry)) max_try = mtry
    if (max_try < 1) status = BAD_INPUT

    if (status == BAD_INPUT) return

   !!!
   !!! BEGIN TIME STEPPING

    step = 0
    STEP_LOOP: do

      h = hnext
      tlast = most_recent_time(this%uhist)

      !! Check for a normal return before proceeding.
      if (t_out <= tlast) then
        status = SOLVED_TO_TOUT
        exit STEP_LOOP
      end if
      if (step >= max_step) then
        status = SOLVED_TO_NSTEP
        exit STEP_LOOP
      end if

      step = step + 1
      
      call bdf2_step (this, pcfun, updpc, enorm, h, h_min, max_try, u, hnext, status)
      if (status /= 0) exit STEP_LOOP   ! step failed

      !! BDF2 step was successful; commit the solution.
      call commit_solution (this, h, u)

      !! Set the next step size.
      hnext = min(h_max, hnext)
      if (this%verbose) write(this%unit,fmt=1) hnext/h

    end do STEP_LOOP

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine bdf2_step_driver

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COMMIT_SOLUTION
 !!

  subroutine commit_solution (this, h, u)

    type(state), intent(inout) :: this
    real(r8), intent(in) :: h, u(:)

    real(r8) :: t

    t = h + most_recent_time(this%uhist)
    call record_solution (this%uhist, t, u)

    this%hlast = h
    this%seq = this%seq + 1
    this%freeze_count = max(0, this%freeze_count - 1)

    this%hmin = min(h, this%hmin)
    this%hmax = max(h, this%hmax)

  end subroutine commit_solution

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF2_STEP
 !!

  subroutine bdf2_step (this, pcfun, updpc, enorm, h, hmin, mtry, u, hnext, errc)

    type(state), intent(inout) :: this
    real(r8), intent(inout) :: h
    real(r8), intent(in)    :: hmin
    integer,  intent(in)    :: mtry
    real(r8), intent(out)   :: u(:), hnext
    integer,  intent(out)   :: errc

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    integer :: try

    errc = 0
    if (hmin < 0.0_r8 .or. h < hmin) errc = BAD_INPUT
    if (mtry < 1) errc = BAD_INPUT
    if (errc == BAD_INPUT) return

    try = 0
    do
      try = try + 1

      !! Check for too many attempts at a single step.
      if (try > mtry) then
        errc  = STEP_FAILED
        exit
      end if

      !! Check for a too-small step size.
      if (h < hmin) then
        errc  = STEP_SIZE_TOO_SMALL
        exit
      end if

      !! Attempt a BDF2 step.
      call bdf2_step_simple (this, pcfun, updpc, enorm, h, u, hnext, errc)
      if (errc == 0) exit

      !! Step failed; try again with the suggested step size.
      if (this%verbose) write(this%unit,fmt=1) hnext/h
      h = hnext

    end do

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine bdf2_step

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF2_STEP_SIMPLE
 !!

  subroutine bdf2_step_simple (this, pcfun, updpc, enorm, h, u, hnext, errc, ectrl)

    type(state), intent(inout) :: this
    real(r8), intent(in)  :: h
    real(r8), intent(out) :: u(:)
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc
    logical,  intent(in), optional :: ectrl

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    ASSERT( this%seq >= 0 )
    ASSERT( size(u) == this%n )
    INSIST( h > 0 )

    select case (this%seq)
    case (0)
      call trap_step_one (this, h, u, hnext, errc, pcfun, updpc, enorm)
    case (1:2)
      call bdf2_step_gen (this, h, u, hnext, errc, pcfun, updpc, enorm, .false.)
    case (3:)
      call bdf2_step_gen (this, h, u, hnext, errc, pcfun, updpc, enorm, ectrl)
    end select

  end subroutine bdf2_step_simple

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF2_STEP_GEN
 !!
 !! This auxillary routine advances the DAE system one step using the second
 !! order backwards difference method with a step size of H.  This routine is
 !! used to compute a general step once the BDF2 integrator is bootstrapped
 !! with sufficient history data.  A local truncation error is estimated and
 !! the step fails if this error is too large.  When the step is successful
 !! a new step size is computed and returned in HNEXT.
 !!
 !! The DAE state is not modified except for data associated with the
 !! preconditioner.  It is up to the caller to commit the computed solution
 !! and advance the state.
 !!

  subroutine bdf2_step_gen (this, h, u, hnext, errc, pcfun, updpc, enorm, ectrl)

    type(state), intent(inout) :: this
    real(r8), intent(in)  :: h
    real(r8), intent(out) :: u(:)
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc
    logical, intent(in), optional :: ectrl

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    real(r8) :: eta, etah, t, t0, tlast, perr, dt(3)
    real(r8) :: u0(size(u)), up(size(u))
    logical  :: fresh_pc, predictor_error

    tlast = most_recent_time(this%uhist)
    t = tlast + h
    eta = (this%hlast + h) / (this%hlast + 2.0_r8 * h)
    etah = eta * h
    t0 = tlast + (1.0_r8 - eta)*h

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h, etah

    fresh_pc = .false.

    !! If the PC step size is too different than the current step size we tag
    !! it as unusable in order to preempt a possible nonlinear solve failure.
    if (this%usable_pc) then
      if (this%hpc/etah > 1.0_r8 + MARGIN) this%usable_pc = .false.
      if (etah/this%hpc > 1.0_r8 + MARGIN) this%usable_pc = .false.
    end if

    !! Predicted solution and base point for BCE step.
    call interpolate_solution (this%uhist, t,  up, order=2)
    call interpolate_solution (this%uhist, t0, u0, order=1)

    BCE: do

      !! Update the preconditioner if necessary.
      if (.not.this%usable_pc) then
        this%updpc_calls = this%updpc_calls + 1
        call updpc (t, up, etah, errc)
        if (errc /= 0) then ! update failed; cut h and return error condition.
          this%updpc_failed = this%updpc_failed + 1
          if (this%verbose) write(this%unit,fmt=2) t, etah
          hnext = 0.5_r8 * h
          errc = -1
          return
        end if
        if (this%verbose) write(this%unit,fmt=3) t, etah
        this%hpc = etah
        this%usable_pc = .true.
        fresh_pc = .true.
      end if

      !! Solve the nonlinear BCE system.
      u = up ! Initial solution guess is the predictor.
      call solve_bce_ain (this, t, etah, u0, u, errc, pcfun, enorm)
      if (errc == 0) exit BCE ! the BCE step was successful.

      if (fresh_pc) then ! preconditioner was fresh; cut h and return error condition.
        this%failed_bce = this%failed_bce + 1
        hnext = 0.5_r8 * h
        this%freeze_count = 1
        errc = 1
        return
      else ! update the preconditioner and retry the nonlinear solve.
        this%retried_bce = this%retried_bce + 1
        this%usable_pc = .false.
        cycle BCE
      end if

    end do BCE

    predictor_error = .true.
    if (present(ectrl)) predictor_error = ectrl

    if (predictor_error) then

      !! Predictor error control.
      u0 = u - up
      perr = enorm(u, u0)
      if (perr < 4.0_r8) then ! accept the step.
        if (this%verbose) write(this%unit,fmt=4) perr
        errc = 0
      else ! reject the step; cut h and return error condition.
        this%rejected_steps = this%rejected_steps + 1
        if (this%verbose) write(this%unit,fmt=5) perr
        hnext = 0.5_r8 * h
        this%freeze_count = 1
        errc = 2
        return
      end if

      !! Select the next step size based on the predictor error and past step
      !! size history, but don't change the step size by too great a factor.
      dt(1) = h
      dt(2:) = h + time_deltas(this%uhist)
      call select_step_size (dt, perr, hnext)
      hnext = max(RMIN*h, min(RMAX*h, hnext))
      if (this%freeze_count /= 0) hnext = min(h, hnext)

    else

      if (this%verbose) write(this%unit,fmt=6)
      hnext = h
      errc = 0

    end if

    1 format(/,'BDF2 step ',i6,': T=',es12.5,', H=',es12.5,', ETAH=',es12.5)
    2 format(2x,'Preconditioner update FAILED at T=',es12.5,', ETAH=',es12.5)
    3 format(2x,'Preconditioner updated at T=',es12.5,', ETAH=',es12.5)
    4 format(2x,'Step accepted: perr=',es12.5)
    5 format(2x,'Step REJECTED: perr=',es12.5)
    6 format(2x,'Step accepted: no local error control')

  end subroutine bdf2_step_gen

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! TRAP_STEP_ONE
 !!
 !! This auxillary routine advances the DAE system one step using the second
 !! order trapezoid method, and is intended to be used to compute the initial
 !! step when bootstrapping the BDF2 integrator.  It does not try to estimate
 !! and control the local truncation error; if the nonlinear solve succeeds
 !! the step is considered successful.
 !!
 !! The DAE state is not modified except for data associated with the
 !! preconditioner.  It is up to the caller to commit the computed solution
 !! and advance the state.
 !!

  subroutine trap_step_one (this, h, u, hnext, errc, pcfun, updpc, enorm)

    type(state), intent(inout) :: this
    real(r8), intent(in)  :: h
    real(r8), intent(out) :: u(:)
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    real(r8) :: etah, t, t0, tlast
    real(r8) :: u0(size(u))

    ASSERT( size(u) == this%n )

    tlast = most_recent_time(this%uhist)
    t = tlast + h
    etah = 0.5_r8 * h
    t0 = tlast + etah

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h, etah

    !! Predicted solution and base point for the BCE step.
    call interpolate_solution (this%uhist, t,  u,  order=1)
    call interpolate_solution (this%uhist, t0, u0, order=1)

    !! Update the preconditioner.
    this%updpc_calls = this%updpc_calls + 1
    call updpc (t, u, etah, errc)
    if (errc /= 0) then
      this%updpc_failed = this%updpc_failed + 1
      if (this%verbose) write(this%unit,fmt=2) t, etah
      hnext = 0.1_r8 * h  ! want to quickly find an acceptably small step size.
      errc = -1
      return
    end if
    if (this%verbose) write(this%unit,fmt=3) t, etah
    this%hpc = etah
    this%usable_pc = .true.

    !! Solve the nonlinear BCE system.
    call solve_bce_ain (this, t, etah, u0, u, errc, pcfun, enorm)
    if (errc /= 0) then
      this%failed_bce = this%failed_bce + 1
      hnext = 0.1_r8 * h  ! want to quickly find an acceptably small step size.
      this%freeze_count = 1
      errc = 1
      return
    end if

    if (this%verbose) write(this%unit,fmt=4)
    this%freeze_count = 2
    hnext = h
    errc = 0

    1 format(/,'TRAP step ',i6,': T=',es12.5,', H=',es12.5,', ETAH=',es12.5)
    2 format(2x,'Preconditioner update FAILED at T=',es12.5,', ETAH=',es12.5)
    3 format(2x,'Preconditioner updated at T=',es12.5,', ETAH=',es12.5)
    4 format(2x,'Step accepted: no local error control')

  end subroutine trap_step_one

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  SELECT_STEP_SIZE -- Choose a new time step.
 !!

  subroutine select_step_size (dt, perr, h)

    real(r8), intent(in)  :: dt(:), perr
    real(r8), intent(out) :: h

    real(r8), parameter :: tol = 0.001_r8
    real(r8) :: a, dh, phi, dphi

    ASSERT( size(dt) == 3 )

    a = 0.5_r8*dt(1)*dt(2)*dt(3)/max(perr,0.001_r8)
    h = dt(1)

    do ! until converged -- DANGEROUS!

      phi  = h*(h + dt(1))*(h + dt(2)) - a
      dphi = (2.0_r8*h + dt(1))*(h + dt(2)) + h*(h + dt(1))

      dh = phi / dphi
      h = h - dh
      if (abs(dh) / h < tol) exit

    end do

  end subroutine select_step_size
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF1_STEP
 !!

  subroutine bdf1_step (this, h, u, errc, pcfun, updpc, enorm)

    type(state), intent(inout) :: this
    real(r8), intent(in)  :: h
    real(r8), intent(out) :: u(:)
    integer,  intent(out) :: errc

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    real(r8) :: t, t0, tlast
    real(r8) :: up(size(u))
    real(r8), pointer :: u0(:)
    logical :: fresh_pc

    tlast = most_recent_time(this%uhist)
    t = tlast + h
    t0 = tlast

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h

    fresh_pc = .false.

    !! If the PC step size is too different than the current step size we tag
    !! it as unusable in order to preempt a possible nonlinear solve failure.
    if (this%usable_pc) then
      if (this%hpc/h > 1.0_r8 + MARGIN) this%usable_pc = .false.
      if (h/this%hpc > 1.0_r8 + MARGIN) this%usable_pc = .false.
    end if
    
    !! Force an update of the preconditioner.  This is in lieu of an effective
    !! strategy for balancing the gain of reusing an old preconditioner for
    !! several steps against the loss of a slowly convergent nonlinear iteration.
    this%usable_pc = .false.

    !! Predicted solution and base point for BCE step.
    call interpolate_solution (this%uhist, t, up, order=1)
    u0 => most_recent_solution(this%uhist)

    BCE: do

      !! Update the preconditioner if necessary.
      if (.not.this%usable_pc) then
        this%updpc_calls = this%updpc_calls + 1
        call updpc (t, up, h, errc)
        if (errc /= 0) then ! update failed; cut h and return error condition.
          this%updpc_failed = this%updpc_failed + 1
          if (this%verbose) write(this%unit,fmt=2) t
          errc = -1
          return
        end if
        if (this%verbose) write(this%unit,fmt=3) t
        this%hpc = h
        this%usable_pc = .true.
        fresh_pc = .true.
      end if

      !! Solve the nonlinear BCE system.
      u = up ! Initial solution guess is the predictor.
      call solve_bce_var (this, t, h, u0, u, errc, pcfun, enorm)
      if (errc == 0) exit BCE ! the BCE step was successful.

      if (fresh_pc) then ! preconditioner was fresh; return error condition.
        this%failed_bce = this%failed_bce + 1
        if (this%verbose) write(this%unit,fmt=4)
        errc = 1
        return
      else ! update the preconditioner and retry the nonlinear solve.
        this%retried_bce = this%retried_bce + 1
        this%usable_pc = .false.
        cycle BCE
      end if

    end do BCE

    if (this%verbose) write(this%unit,fmt=5)
    errc = 0

    1 format(/,'BDF1 step ',i6,': T=',es12.5,', H=',es12.5)
    2 format(2x,'Preconditioner update FAILED at T=',es12.5)
    3 format(2x,'Preconditioner updated at T=',es12.5)
    4 format(2x,'Step FAILED')
    5 format(2x,'Step accepted')

  end subroutine bdf1_step

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOLVE_BCE_AIN -- Solve the Backward Cauchy-Euler system using AIN.
 !!
 !! The backward Cauchy-Euler (BCE) method applied to the implicit DAE
 !!
 !!     f(t,u,u') = 0
 !!
 !! yields a nonlinear system of equations for advancing the solution from a
 !! given state u0 at time t - h to the unknown solution u at time t,
 !!
 !!     f(t,u,(u-u0)/h) = 0.
 !!
 !! This subroutine solves this nonlinear system using an accelerated fixed
 !! point iteration [1] for the preconditioned system
 !! g(u) = pc(f(t,u,(u-u0)/h)) = 0:
 !!
 !!    u given
 !!    Do until converged:
 !!      du <-- g(u)
 !!      du <-- NKA(du)
 !!      u  <-- u - du
 !!    End do
 !!
 !! The procedure NKA uses information about g' gleaned from the unaccelerated
 !! correction du=g(u) and previous g values to compute an improved correction.
 !! The preconditioning function pc() is typically an approximate solution
 !! of the Newton correction equation  J*du = f(t,u,(u-u0)/h) where J is an
 !! approximation to the Jacobian of f(t,u,(u-u0)/h) as a function of u.  Thus
 !! this method can be regarded as an accelerated inexact Newton (AIN) method.
 !!
 !! The dummy procedure PCFUN evaluates the preconditioned function g.
 !!
 !! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
 !!     weighted moving finite element code I: in one dimension", SIAM J.
 !!     Sci. Comput;, 19 (1998), pp. 728-765..
 !!

  subroutine solve_bce_ain (this, t, h, u0, u, errc, pcfun, enorm)

    type(state), intent(inout) :: this
    real(r8), intent(in)    :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out)   :: errc

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    integer  :: itr
    real(r8) :: error, du(size(u))

    call nka_restart (this%fpa)

    itr = 0
    do

      !! Check for too many nonlinear iterations.
      if (itr >= this%mitr) then
        if (this%verbose) write(this%unit,fmt=1) itr, error
        errc = 1
        exit
      end if

      itr = itr + 1

      !! Evaluate the preconditioned nonlinear function.
      this%pcfun_calls = this%pcfun_calls + 1
      call pcfun (t, u, (u-u0)/h, du)

      !! Accelerated correction.
      call nka_accel_update (this%fpa, du, dp=pardp)

      !! Next solution iterate and error estimate.
      u  = u - du
      error = enorm(u, du)
      if (this%verbose) write(this%unit,fmt=3) itr, error

      !! Check for convergence.
      if (((error < this%ntol) .and. (itr > 1)) .or. (error < 0.01_r8 * this%ntol)) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        errc = 0
        exit
      end if

    end do

    1 format(2x,'AIN BCE solve FAILED: ',i3,' iterations (max), error=',es12.5)
    2 format(2x,'AIN BCE solve succeeded: ',i3,' iterations, error=',es12.5)
    3 format(2x,i3,': error=',es12.5)

  end subroutine solve_bce_ain

  function pardp (a, b) result (dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function pardp
  
  !!
  !! SOLVE_BCE_VAR -- a test variation that uses different convergence
  !! criteria more suited to the BDF1 stepping without trunctation
  !! error control.
  !! 

  subroutine solve_bce_var (this, t, h, u0, u, errc, pcfun, enorm)

    type(state), intent(inout) :: this
    real(r8), intent(in)    :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out)   :: errc

#include "bdf2_dae_pcfun_if.fpp"
#include "bdf2_dae_enorm_if.fpp"

    integer  :: itr
    real(r8) :: error, error0, du(size(u)), conv_rate

    call nka_restart (this%fpa)

    itr = 0

    !! Evaluate the initial preconditioned nonlinear function.
    this%pcfun_calls = this%pcfun_calls + 1
    call pcfun (t, u, (u-u0)/h, du)
    
    !! Initial error.
    error0 = enorm(u, du)
    if (this%verbose) write(this%unit,fmt=3) itr, error0
    if (error0 == 0.0_r8) return

    do

      !! Check for too many nonlinear iterations.
      if (itr >= this%mitr) then
        if (this%verbose) write(this%unit,fmt=1) itr, error
        errc = 1
        exit
      end if

      itr = itr + 1

      !! Accelerated correction and next solution iterate.
      call nka_accel_update (this%fpa, du, dp=pardp)
      u  = u - du
      
      !! Evaluate the preconditioned nonlinear function.
      this%pcfun_calls = this%pcfun_calls + 1
      call pcfun (t, u, (u-u0)/h, du)

      !! Residual error and convergence check.  Recall that ENORM returns scaled
      !! values; anything less than 1 is considered acceptably small.  We stop
      !! iterating when the residual error is less than 1 or when the initial
      !! residual error has been reduced by the factor NTOL.
      error = enorm(u, du)
      conv_rate = (error/error0)**(1.0_r8/itr)
      if (this%verbose) write(this%unit,fmt=3) itr, error, conv_rate
      if ((itr >= 2) .and. ((error < 1.0_r8) .or. (error < this%ntol * error0))) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        errc = 0
        exit
      end if
      
      !! Exit with an error if the convergence rate is too slow
      !if (conv_rate > 0.95_r8) then
      !  if (this%verbose) write(this%unit,fmt=4) itr, conv_rate
      !  errc = 1
      !  exit
      !end if

    end do

    1 format(2x,'AIN BCE solve FAILED: ',i3,' iterations (max), error=',es12.5)
    2 format(2x,'AIN BCE solve succeeded: ',i3,' iterations, error=',es12.5)
    3 format(2x,i3,': error=',es12.5,', convergence rate=',es12.5)
    4 format(2x,'AIN BCE solve FAILED: ',i3,' iterations, convergence rate=',es12.5)

  end subroutine solve_bce_var

  
  subroutine bdf1_step_new (this, h, u, errc, fun, pc, updpc, fnorm)

    type(state), intent(inout) :: this
    real(r8), intent(in)  :: h
    real(r8), intent(out) :: u(:)
    integer,  intent(out) :: errc

#include "bdf2_dae_fun_if.fpp"
#include "bdf2_dae_pc_if.fpp"
#include "bdf2_dae_updpc_if.fpp"
#include "bdf2_dae_fnorm_if.fpp"

    real(r8) :: t, tlast
    real(r8), pointer :: u0(:)

    tlast = most_recent_time(this%uhist)
    t = tlast + h

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h

    !! Predicted solution and base point for BCE step.
    call interpolate_solution (this%uhist, t, u, order=1)
    u0 => most_recent_solution(this%uhist)

    !! Update the preconditioner.
    this%updpc_calls = this%updpc_calls + 1
    call updpc (t, u, h, errc)
    if (errc /= 0) then
      this%updpc_failed = this%updpc_failed + 1
      if (this%verbose) write(this%unit,fmt=2)
      errc = -1
      return
    end if

    !! Solve the nonlinear BCE system.
    call solve_bce3 (this, t, h, u0, u, errc, fun, pc, fnorm)
    if (errc == 0) then ! the BCE step was successful.
      if (this%verbose) write(this%unit,fmt=5)
    else if (errc /= 0) then ! preconditioner was fresh; return error condition.
      this%failed_bce = this%failed_bce + 1
      if (this%verbose) write(this%unit,fmt=4)
      errc = -1
      return
    end if

    1 format(/,'BDF1 step ',i6,': T=',es12.5,', H=',es12.5)
    2 format(2x,'Preconditioner update FAILED')
    4 format(2x,'Step FAILED')
    5 format(2x,'Step accepted')

  end subroutine bdf1_step_new

  !!
  !! SOLVE_BCE3 -- yet another BCE step solver that uses an entirely different
  !! convergence criteria.
  !! 

  subroutine solve_bce3 (this, t, h, u0, u, errc, fun, pc, fnorm)

    type(state), intent(inout) :: this
    real(r8), intent(in)    :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out)   :: errc

#include "bdf2_dae_fun_if.fpp"
#include "bdf2_dae_pc_if.fpp"
#include "bdf2_dae_fnorm_if.fpp"

    integer  :: itr
    real(r8) :: f(size(u)), udot(size(u)), error

    call nka_restart (this%fpa)

    itr = 0
    
    !! Compute the initial function value and norm.
    this%pcfun_calls = this%pcfun_calls + 1
    udot = (u - u0)/h
    call fun (t, u, udot, f)
    call fnorm (t, u, udot, f)
    
    do

      itr = itr + 1
    
      !! Compute the next solution iterate.
      call pc (t, u, udot, f)
      call nka_accel_update (this%fpa, f, dp=pardp)
      u = u - f
      udot = (u - u0)/h
      
      !! Compute the function value and norm.
      this%pcfun_calls = this%pcfun_calls + 1
      call fun (t, u, udot, f)
      call fnorm (t, u, udot, f, error)
      if (this%verbose) write(this%unit,fmt=3) itr, error
      
      !! Convergence check and iteration control.
      if (itr >= 2 .and. error < 1.0_r8) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        errc = 0
        exit
      else if (itr >= this%mitr) then
        if (this%verbose) write(this%unit,fmt=1) itr, error
        errc = 1
        exit
      end if

    end do

    1 format(2x,'AIN BCE solve FAILED: ',i3,' iterations (max), error=',es12.5)
    2 format(2x,'AIN BCE solve succeeded: ',i3,' iterations, error=',es12.5)
    3 format(2x,i3,': error=',es12.5)

  end subroutine solve_bce3

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_VERBOSE_STEPPING / SET_QUIET_STEPPING
 !!

  subroutine set_verbose_stepping (this, unit)
    use parallel_communication, only: is_IOP
    type(state), intent(inout) :: this
    integer, intent(in) :: unit
    this%unit = unit
    this%verbose = is_IOP
  end subroutine set_verbose_stepping

  subroutine set_quiet_stepping (this)
    type(state), intent(inout) :: this
    this%verbose = .false.
  end subroutine set_quiet_stepping

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  State inquiry procedures...
 !!

  function last_time (this) result (t)
    type(state), intent(in) :: this
    real(r8) :: t
    t = most_recent_time(this%uhist)
  end function last_time

  function last_solution (this) result (uptr)
    type(state), intent(in) :: this
    real(r8), pointer :: uptr(:)
    uptr => most_recent_solution(this%uhist)
  end function last_solution

  function last_step_size (this) result (h)
    type(state), intent(in) :: this
    real(r8) :: h
    h = this%hlast
  end function last_step_size

  function interpolated_solution (this, t) result (u)
    type(state), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8) :: u(this%n)
    call interpolate_solution (this%uhist, t, u)
  end function interpolated_solution

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_BDF2_STEPPING_STATISTICS
 !!

  subroutine write_bdf2_stepping_statistics (this, unit)
    use parallel_communication, only: is_IOP
    type(state), intent(in) :: this
    integer, intent(in) :: unit
    if (.not.is_IOP) return
    write(unit,fmt='(/,a,i6,a,es11.5,a,es9.3)') &
      'STEP=', this%seq, ', T=', most_recent_time(this%uhist), ', H=', this%hlast
    write(unit,fmt='(a,i7.7,":",i5.5,a,5(i4.4,:,":"))') &
      'NFUN:NPC=', this%pcfun_calls, this%updpc_calls, &
      ', NPCF:NNR:NNF:NSR=', this%updpc_failed, &
      this%retried_bce, this%failed_bce, this%rejected_steps
  end subroutine write_bdf2_stepping_statistics

  subroutine get_bdf2_stepping_statistics (this, counters)
    type(state), intent(in) :: this
    integer, intent(out) :: counters(:)
    ASSERT( size(counters) == 6 )
    counters(1) = this%pcfun_calls
    counters(2) = this%updpc_calls
    counters(3) = this%updpc_failed
    counters(4) = this%retried_bce
    counters(5) = this%failed_bce
    counters(6) = this%rejected_steps
  end subroutine get_bdf2_stepping_statistics

end module bdf2_dae
