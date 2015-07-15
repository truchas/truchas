#include "f90_assert.fpp"

module bdf2_integrator

  use bdf2_kinds
  use bdf2_controller
  use bdf2_profiling
  use solution_history
  use nka_type
  implicit none
  private

  public :: bdf2_create_state, bdf2_init_state
  public :: bdf2_integrate, bdf2_inquire, bdf2_write_profile
  public :: bdf2_interpolate_solution
  
  public :: bdf2_solution, bdf2_solution_time
  
  public :: bdf2_profile, set_solver_messages
  public :: destroy
  
  !! Objects from bdf2_controller to export.
  public :: bdf2_control, bdf2_set_param
  
  ! Return codes.
  integer, parameter, public :: FAIL_ON_STEP = -1
  integer, parameter, public :: SMALL_H_FAIL = -2
  integer, parameter, public :: FAIL_ON_START = -3
  integer, parameter, public :: SOLN_AT_TOUT = 2
  integer, parameter, public :: SOLN_AT_STEP = 3
  integer, parameter, public :: BAD_INPUT = -4

  !! Integrator state identifiers
  integer, parameter, private :: STATE_UNDEF = -1
  integer, parameter, private :: STATE_INIT  = 0
  integer, parameter, private :: STATE_START = 1

  type, public :: bdf2_state
    private
    integer :: n
    type(history) :: uhist ! solution history
    real(kind=rk) :: hnext        ! next step size
    real(kind=rk) :: hlast        ! last step size
    real(kind=rk) :: htwo         ! sum of previous two step sizes

    real(kind=rk), pointer :: dfdu(:,:) => null()
    integer       :: freeze_count   ! Don't increase step size for this number of steps

    type(nka), pointer :: fpa => null()
    
    integer :: state = STATE_UNDEF
    type(bdf2_profile), pointer :: profile => null()
  end type bdf2_state

  real(kind=rk), parameter, private :: RMIN = 0.25_rk
  real(kind=rk), parameter, private :: RMAX = 4.0_rk
  
  interface destroy
    module procedure destroy_bdf2_state
  end interface

contains

  subroutine bdf2_create_state (state, n)
    type(bdf2_state), intent(out) :: state
    integer, intent(in) :: n
    ASSERT( n > 0 )
    state%n = n
    call create_history (state%uhist, 3, n)
  end subroutine bdf2_create_state
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF2_INITIALIZE
 !!
 
  subroutine bdf2_init_state (state, u, t, hstart, profile)
  
    type(bdf2_state), intent(inout) :: state
    real(kind=rk), intent(in) :: u(:), t, hstart
    logical, intent(in), optional :: profile
  
    ASSERT( size(u) == state%n )
    
    call flush_history (state%uhist, t, u)
    state%hnext = hstart
    state%state = STATE_INIT
    
    !allocate(state%dfdu(state%n,state%n))
    
    if (present(profile)) then
      if (profile) then
        if (associated(state%profile)) then
          call reset_profile (state%profile)
        else
          allocate(state%profile) ! default initialized
        end if
      else if (associated(state%profile)) then
        deallocate(state%profile)
      end if
    else if (associated(state%profile)) then
      deallocate(state%profile)
    end if

  end subroutine bdf2_init_state
  
  subroutine destroy_bdf2_state (state)
    type(bdf2_state), intent(inout) :: state
    call destroy (state%uhist)
    if (associated(state%fpa)) deallocate(state%fpa)
    if (associated(state%profile)) deallocate(state%profile)
  end subroutine destroy_bdf2_state
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  BDF2_INTEGRATE
 !!

  subroutine bdf2_integrate (state, control, nstep, tout, stat, &
    rhs, jac, schk, user)

    type(bdf2_state),   intent(inout) :: state
    type(bdf2_control), intent(in)    :: control
    integer,       intent(in), optional :: nstep
    real(kind=rk), intent(in), optional :: tout
    integer, intent(out)    :: stat

#include "rhs_interface.fpp"
#include "jac_interface.fpp"
#include "schk_interface.fpp"
#include "user_interface.fpp"
    optional :: jac, schk, user

    logical  :: renew_jacobian, stale_jacobian
    integer  :: max_step, try, step, mv
    real(kind=rk) :: perr, eta, etah, t, t0, h, tlast
    real(kind=rk), dimension(state%n) :: u0, up, u
    
    !! Set the maximum number of time steps.
    if (present(nstep)) then  ! use this value
      if (nstep < 1) then
        stat = BAD_INPUT
        return
      end if
      max_step = nstep
    else if (present(tout)) then  ! no limit to the number of steps
      max_step = huge(1)
    else  ! neither NSTEP or TOUT was specified, so just take a single step
      max_step = 1
    end if
    
    !! Check the target integration time (if any) for appropriateness.
    if (present(tout)) then
      if (tout <= most_recent_time(state%uhist)) then
        stat = BAD_INPUT
        return
      end if
    end if
    
    !! 
    if (present(jac)) then
      if (.not.associated(state%dfdu)) then
        allocate(state%dfdu(state%n,state%n))
        renew_jacobian = .true.
      end if
    else  ! jacobian-free procedure
      if (associated(state%dfdu)) then
        deallocate(state%dfdu)
      end if
      stale_jacobian = .false.
      renew_jacobian = .false.
    end if
    
    !! Initialize the NKA structure.
    if (control%mvec > 0) then
      if (.not.associated(state%fpa)) then
        allocate(state%fpa)
        mv = min(control%mvec, control%mitr-1, state%n)
        call state%fpa%init (size(u), mv)
        call state%fpa%set_vec_tol (control%vtol)
      end if
    else
      if (associated(state%fpa)) deallocate(state%fpa)
    end if
    
    step = 0
    
    if (state%state == STATE_INIT) then ! starting procedure
      t = most_recent_time(state%uhist)
      u = most_recent_solution(state%uhist)
      if (present(user)) call user (t, u)
      call start_bdf2 (state, control, stat, rhs, jac, schk)
      if (stat /= 0) then           ! Starting procedure failed.
        stat = FAIL_ON_START
        return
      end if
      state%state = STATE_START
      if (present(jac)) renew_jacobian = .false.
      state%freeze_count = 1
      step = 1
    end if
    
    t = most_recent_time(state%uhist)
    u = most_recent_solution(state%uhist)
    
    BDF2_STEP: do
    
      !! Check for a normal return before proceding.
      if (present(tout)) then
        if (t >= tout) then
          stat = SOLN_AT_TOUT
          exit BDF2_STEP
        end if
      end if
      if (step >= max_step) then
        stat = SOLN_AT_STEP
        exit BDF2_STEP
      end if
      
      if (present(user)) call user (t, u)

      tlast = t
      h = state%hnext
      try = 0
      step = step + 1

     !!!
     !!! Atempt a BDF2 Step

      ATTEMPT: do

        try = try + 1

        !! Check for too many attempts at a single BDF2 step.
        if (try > control%mtry) then
          state%hnext = h
          stat  = FAIL_ON_STEP
          exit BDF2_STEP
        end if

        !! Check for a too-small step size.
        if (h < control%hmin) then
          state%hnext = h
          stat  = SMALL_H_FAIL
          exit BDF2_STEP
        end if

        if (associated(state%profile)) call event_bdf2_step (state%profile, try, tlast, h)
        
        t = tlast + h
        eta = (state%hlast + h) / (state%hlast + 2.0_rk * h)
        etah = eta * h
        t0 = tlast + (1.0_rk - eta)*h
        if (present(jac)) stale_jacobian = .true.

        !! Predicted solution and base point for BCE step
        call interpolate_solution (state%uhist, t,  up, order=2)
        call interpolate_solution (state%uhist, t0, u0, order=1)

        !! Check the predicted solution for admissibility.
        if (present(schk)) then
          call schk (up, 0, stat)
          if (stat /= 0) then ! it's bad; cut h and retry.
            h = 0.5_rk * h
            state%freeze_count = 1
            if (associated(state%profile)) call event_bad_predictor (state%profile, h)
            cycle ATTEMPT
          end if
        end if

        BCE: do

          !! Evaluate the Jacobian.
          if (renew_jacobian) then
            if (associated(state%profile)) call event_jacobian_evaluated (state%profile)
            call jac (t, up, state%dfdu, stat)
            if (stat /= 0) then ! evaluation failed; cut h and retry.
              h = 0.5_rk * h
              renew_jacobian = .true.
              if (associated(state%profile)) call event_failed_jacobian (state%profile, h)
              cycle ATTEMPT
            end if
            renew_jacobian = .false.
            stale_jacobian = .false.
          end if

          !! Compute the BCE step: the nonlinear solve.
          call bce_step (state, control, t, etah, u0, up, u, stat, rhs, schk)
          if (stat == 0) exit BCE ! the BCE step was successful.

          if (stale_jacobian) then ! update jacobian and retry BCE step.
            renew_jacobian = .true.
            if (associated(state%profile)) call event_retry_bce_step (state%profile)
            cycle BCE
          else ! jacobian was fresh; cut h and retry.
            h = 0.5_rk * h
            state%freeze_count = 1
            if (associated(state%profile)) call event_failed_bce_step (state%profile, h)
            cycle ATTEMPT
          end if

        end do BCE

        !! Predictor error control.
        u0 = u - up
        perr = error_norm(control, u, u0)
        if (perr < 4.0_rk) then ! accept the step.
          if (associated(state%profile)) call event_step_accepted (state%profile, perr)
          exit ATTEMPT
        else ! reject the step; cut h and retry.
          h = 0.5_rk * h
          state%freeze_count = 1
          if (associated(state%profile)) call event_step_rejected (state%profile, perr, h)
          cycle ATTEMPT
        end if

      end do ATTEMPT

      !! BDF2 step accepted; commit the solution.
      call record_solution (state%uhist, t, u)
      state%hnext = h

      !! Select the next step size.
      call select_step_size (state, perr)
      state%hnext = max(RMIN*state%hlast, min(control%hmax, RMAX*state%hlast, state%hnext))
      if (state%freeze_count /= 0) state%hnext = min(state%hlast, state%hnext)
      state%freeze_count = max(0, state%freeze_count - 1)
      if (associated(state%profile)) call event_new_step_size (state%profile, state%hlast, state%hnext)
      
    end do BDF2_STEP

  end subroutine bdf2_integrate

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  SELECT_STEP_SIZE -- Choose a new time step.
 !!
 
  subroutine select_step_size (state, error)
  
    type(bdf2_state), intent(inout) :: state
    real(kind=rk),    intent(in)    :: error  ! predictor error
    
    real(kind=rk), parameter :: tol = 0.001_rk
    real(kind=rk) :: a, dh, phi, dphi, h

    a = 0.5_rk*state%hnext*(state%hnext+state%hlast)*(state%hnext+state%htwo)/max(error,0.001_rk)
    h = state%hnext

    do ! until converged -- DANGEROUS!

      phi  = h * (h + state%hnext) * (h + state%hlast + state%hnext) - a
      dphi = (2.0_rk * h + state%hnext) * (h + state%hlast + state%hnext) + h * (h + state%hnext)

      dh = phi / dphi
      h = h - dh
      if (abs(dh) / h < tol) exit

    end do
    
    state%htwo  = state%hlast + state%hnext
    state%hlast = state%hnext
    state%hnext = h
    
  end subroutine select_step_size

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BCE_STEP -- Backward Cauchy-Euler Step.                                    
 !!
 !! The backward Cauchy-Euler method applied to the ODE system u' = f(t,u)     
 !! yields a nonlinear system of equations of the form                         
 !!
 !!    R(u) = u - u_0 - h*f(t,u) = 0,
 !!
 !! for advancing the solution from a given state u_0 at time t - h to the
 !! solution u at time t.  This auxillary procedure solves this nonlinear
 !! system using an accelerated inexact Newton (AIN) method [1],
 !!
 !!    u given
 !!    Do until converged
 !!      du <-- R(u)
 !!      du <-- NKA(du)                      
 !!      u  <-- u - du
 !!    End do
 !!
 !! The Newton correction equation R'(u) du = R(u), R'(u) = I - h * df/du,
 !! has been approximated by entirely neglecting the Jacobian matrix of f.
 !! The inferior correction du that results is accelerated by the procedure
 !! NKA which uses information about the true R' gleaned from the sequence
 !! of previous R values.
 !!
 !! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
 !!     weighted moving finite element code I: in one dimension", SIAM J.
 !!     Sci. Comput;, 19 (1998), pp. 728-765..
 !!

  subroutine bce_step (state, control, t, h, u0, up, u, stat, rhs, schk)

    type(bdf2_state),   intent(inout) :: state
    type(bdf2_control), intent(in)    :: control
    real(kind=rk),      intent(in)    :: t, h, u0(:), up(:)
    real(kind=rk),      intent(out)   :: u(:)
    integer,            intent(out)   :: stat

#include "rhs_interface.fpp"
#include "schk_interface.fpp"
    optional :: schk

    integer  :: itr
    real(kind=rk) :: error, du(size(u)), ddu(size(u))
    
    if (associated(state%fpa)) call state%fpa%restart
 
    itr = 0
    u = up
    du = 0.0_rk

    do

      itr = itr + 1

      !! Check for too many nonlinear iterations.
      if (itr > control%mitr) then
        if (associated(state%profile)) call event_newton_iteration (state%profile, itr, error)
        stat = 1
        exit
      end if

      !! Evaluate the nonlinear function.
      if (associated(state%profile)) call event_residual_evaluated (state%profile)
      call rhs (t, u, ddu)
      ddu = u - u0 - h*ddu

      !! Accelerated correction.
      if (associated(state%fpa)) call state%fpa%accel_update (ddu)


      !! Next solution iterate.
      du = du + ddu
      u  = up - du

      !! Check the solution iterate for admissibility.
      if (present(schk)) then
        call schk (u, 1, stat)
        if (stat /= 0) then ! iterate is bad; bail.
          stat = 2
          if (associated(state%profile)) call event_failed_newton (state%profile, itr)
          exit
        end if
      end if

      !! Error estimate.
      error = error_norm(control, u, ddu)

      ! Check for convergence.
      if ((error < 0.01_rk * control%ntol) .or. ((error < control%ntol) .and. (itr > 1))) then
        if (associated(state%profile)) call event_newton_iteration (state%profile, itr, error)
        stat = 0
        exit
      end if

    end do

  end subroutine bce_step
  
  
  subroutine ssor_relaxation (dfdu, h, r, omega, nsweep)
  
    real(kind=rk), intent(in)    :: dfdu(:,:), h, omega
    real(kind=rk), intent(inout) :: r(:)
    integer,       intent(in)    :: nsweep
    
    integer :: i, j, n, l
    real(kind=rk) :: s, b(size(r))
    
    n = size(r)
    
    b = r
    !r = 0.0_rk
    
    do l = 1, nsweep
      !! Forward sweep
      do i = 1, n
        s = 0.0_rk
        do j = 1, n
          s = dfdu(i,j)*r(j)
        end do
        s = b(i) - r(i) + h*s
        r(i) = r(i) + omega * s / (1.0_rk - h*dfdu(i,i))
      end do

      !! Backward sweep
      do i = n, 1, -1
        s = 0.0_rk
        do j = 1, n
          s = dfdu(i,j)*r(j)
        end do
        s = b(i) - r(i) + h*s
        r(i) = r(i) + omega * s / (1.0_rk - h*dfdu(i,i))
      end do
    end do
  
  end subroutine ssor_relaxation
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  START_BDF2 -- Starting procedure for BDF2.
 !!

  subroutine start_bdf2 (state, control, stat, rhs, jac, schk)

    type(bdf2_state),   intent(inout) :: state
    type(bdf2_control), intent(in)    :: control
    integer, intent(out)    :: stat

#include "rhs_interface.fpp"
#include "jac_interface.fpp"
#include "schk_interface.fpp"
    optional :: jac, schk


   ! local variable.
    real(kind=rk) :: etah, t, t0, h
    real(kind=rk), dimension(state%n) :: u, u0, up
    
    !! Recreate the solution history
    u = most_recent_solution(state%uhist)
    t = most_recent_time(state%uhist)
    h = state%hnext
    
    call rhs (t, u, u0)
    call flush_history (state%uhist, t-h, u-h*u0)
    call record_solution (state%uhist, t, u)
    
   !!!
   !!!  Step 1:  Trapezoid method with FCE as the predictor.

    if (associated(state%profile)) call event_bdf2_step (state%profile, 1, t, h)

    etah = 0.5_rk * h
    t0 = t + etah
    t = t + h

    ! Predicted solution and backward difference for the BCE step.
    call interpolate_solution (state%uhist, t,  up, order=1)
    call interpolate_solution (state%uhist, t0, u0, order=1)

    ! Check the predicted solution.
    if (present(schk)) then
      call schk (up, 0, stat)
      if (stat /= 0) return
    end if

    ! Evaluate the Jacobian at the predicted solution.
    if (present(jac)) then
      if (associated(state%profile)) call event_jacobian_evaluated (state%profile)
      call jac (t, up, state%dfdu, stat)
      if (stat /= 0) return
    end if

    call bce_step (state, control, t, etah, u0, up, u, stat, rhs, schk)
    if (stat /= 0) return

    ! Commit the solution.
    call record_solution (state%uhist, t, u)
    
    state%hlast = h
    state%htwo = 2.0_rk * h

    stat = 0

  end subroutine start_bdf2

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  BDF2_INQUIRE -- Read internal solver variables.
 !!


  subroutine bdf2_inquire (state, t, h_next, h_last)

    type(bdf2_state), intent(in) :: state
    real(kind=rk), intent(out), optional :: t, h_next, h_last

    if (present(t)) t = most_recent_time(state%uhist)
    if (present(h_next)) h_next = state%hnext
    if (present(h_last)) h_last = state%hlast

  end subroutine bdf2_inquire
  
  function bdf2_solution (state) result (ptr)
    type(bdf2_state), intent(in) :: state
    real(kind=rk), pointer :: ptr(:)
    ptr => most_recent_solution(state%uhist)
  end function bdf2_solution
  
  function bdf2_solution_time (state) result (t)
    type(bdf2_state), intent(in) :: state
    real(kind=rk) :: t
    t = most_recent_time(state%uhist)
  end function bdf2_solution_time
  
  function bdf2_interpolate_solution (state, t) result (u)
    type(bdf2_state), intent(in) :: state
    real(kind=rk), intent(in) :: t
    real(kind=rk) :: u(state%n)
    call interpolate_solution (state%uhist, t, u)
  end function bdf2_interpolate_solution

  subroutine bdf2_write_profile (state, unit)
    type(bdf2_state), intent(in) :: state
    integer, intent(in) :: unit
    if (.not.associated(state%profile)) return
    write(unit,fmt='(/,a,i5,a,es12.5,a,es10.3)') &
      'STEP=', state%profile%step, ', T=', most_recent_time(state%uhist), ', H=', state%hnext
    write(unit,fmt='(a,i6.6,":",i4.4,a,5(i3.3,:,":"))') &
      'NRES:NJAC=', state%profile%residual_calls, state%profile%jacobian_calls, &
      ', NBP:NJF:NNR:NNF:NSR=', state%profile%bad_predictors, state%profile%failed_jacobians, &
      state%profile%retried_bce_steps, state%profile%failed_bce_steps, state%profile%rejected_steps
  end subroutine bdf2_write_profile

end module bdf2_integrator
