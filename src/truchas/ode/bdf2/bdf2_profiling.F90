!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bdf2_profiling

  use bdf2_kinds
  implicit none
  private
  
  public :: set_solver_messages
  public :: reset_profile
  public :: event_bdf2_step
  public :: event_bad_predictor
  public :: event_jacobian_evaluated
  public :: event_failed_jacobian
  public :: event_retry_bce_step
  public :: event_failed_bce_step
  public :: event_step_accepted
  public :: event_step_rejected
  public :: event_new_step_size
  public :: event_newton_iteration
  public :: event_residual_evaluated
  public :: event_failed_newton
  

  type, public :: bdf2_profile
    integer :: msgs = 0
    integer :: unit = 0
    integer :: step = 0
    integer :: residual_calls = 0      ! Number of calls to EVAL_RESIDUAL
    integer :: jacobian_calls = 0      ! Number of calls to EVAL_JACOBIAN
    integer :: bad_predictors = 0      ! Number of bad predicted solutions (CHECK_SOLN)
    integer :: failed_jacobians = 0    ! Number of failed Jacobian evaluations 
    integer :: retried_bce_steps = 0   ! Number of retried BCE steps
    integer :: failed_bce_steps = 0    ! Number of completely failed BCE steps
    integer :: rejected_steps = 0      ! Number of steps rejected on error tolerance
  end type bdf2_profile
      
  ! Diagnostic message formats.
  character(len=*), parameter, private ::                                  &
        FMT_N = "('N:',i4.4,':',i1,':',e12.6,':',e12.6,':')", &
        FMT_J = "('J:')",                                     &
        FMT_A = "('A:',e12.6,':')",                           &
        FMT_H = "('H:',e12.6,':',e12.6,':')",                 &
       FMT_PF = "('PF:',e12.6,':')",                          &
       FMT_JF = "('JF:',e12.6,':')",                          &
       FMT_NR = "('NR:')",                                    &
       FMT_NF = "('NF:',e12.6,':')",                          &
      FMT_PEF = "('PEF:',e12.6,':',e12.6,':')",               &
       FMT_NI = "('NI:',i1,':',e12.6,':',e12.6,':')",         &
      FMT_NIF = "('NIF:',i1,':',e12.6,':')"
      
contains

  subroutine set_solver_messages (profile, level, unit)
  
    type(bdf2_profile), intent(inout) :: profile
    integer, intent(in) :: level, unit

    if (level >= 0) profile%msgs = level
    if (unit  >= 0) profile%unit = unit

  end subroutine set_solver_messages
  
  subroutine reset_profile (profile)
    type(bdf2_profile), intent(out) :: profile
  end subroutine reset_profile

  subroutine event_bdf2_step (profile, try, t, h)
    type(bdf2_profile), intent(inout) :: profile
    integer, intent(in) :: try
    real(kind=rk), intent(in) :: t, h
    if (try == 1) profile%step = 1 + profile%step
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_N) profile%step, try, t, h
  end subroutine event_bdf2_step
  
  subroutine event_bad_predictor (profile, hnext)
    type(bdf2_profile), intent(inout) :: profile
    real(kind=rk), intent(in) :: hnext
    profile%bad_predictors = profile%bad_predictors + 1
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_PF) hnext
  end subroutine event_bad_predictor
  
  subroutine event_failed_jacobian (profile, hnext)
    type(bdf2_profile), intent(inout) :: profile
    real(kind=rk), intent(in) :: hnext
    profile%failed_jacobians = profile%failed_jacobians + 1
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_JF) hnext
  end subroutine event_failed_jacobian
  
  subroutine event_jacobian_evaluated (profile)
    type(bdf2_profile), intent(inout) :: profile
    profile%jacobian_calls = profile%jacobian_calls + 1
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_J)
  end subroutine event_jacobian_evaluated
  
  subroutine event_retry_bce_step (profile)
    type(bdf2_profile), intent(inout) :: profile
    profile%retried_bce_steps = profile%retried_bce_steps + 1
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_NR)
  end subroutine event_retry_bce_step
  
  subroutine event_failed_bce_step (profile, hnext)
    type(bdf2_profile), intent(inout) :: profile
    real(kind=rk), intent(in) :: hnext
    profile%failed_bce_steps = profile%failed_bce_steps + 1
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_NF) hnext
  end subroutine event_failed_bce_step
  
  subroutine event_step_accepted (profile, error)
    type(bdf2_profile), intent(in) :: profile
    real(kind=rk), intent(in) :: error
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_A) error
  end subroutine event_step_accepted
  
  subroutine event_step_rejected (profile, error, hnext)
    type(bdf2_profile), intent(inout) :: profile
    real(kind=rk), intent(in) :: error, hnext
    profile%rejected_steps = profile%rejected_steps + 1
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_PEF) error, hnext
  end subroutine event_step_rejected
  
  subroutine event_new_step_size (profile, hlast, hnext)
    type(bdf2_profile), intent(in) :: profile
    real(kind=rk), intent(in) :: hlast, hnext
    if (profile%msgs > 0) write (profile%unit,fmt=FMT_H) hlast, hnext
  end subroutine event_new_step_size
  
  subroutine event_newton_iteration (profile, itr, error)
    type(bdf2_profile), intent(in) :: profile
    integer, intent(in) :: itr
    real(kind=rk), intent(in) :: error
    if (profile%msgs > 1) write (profile%unit,fmt=FMT_NI) itr, error
  end subroutine event_newton_iteration
  
  subroutine event_residual_evaluated (profile)
    type(bdf2_profile), intent(inout) :: profile
    profile%residual_calls = profile%residual_calls + 1
  end subroutine event_residual_evaluated
  
  subroutine event_failed_newton (profile, itr)
    type(bdf2_profile), intent(in) :: profile
    integer, intent(in) :: itr
    if (profile%msgs > 1) write (profile%unit,fmt=FMT_NIF) itr
  end subroutine event_failed_newton
  
end module bdf2_profiling
