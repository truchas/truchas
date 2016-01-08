!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vanderpol

  use bdf2_kinds
  implicit none
  private
  
  public :: rhs, user, jac

  real(kind=rk), parameter :: A = 1000.0_rk

contains

  subroutine rhs (t, u, f)
  
    real(kind=rk), intent(in)  :: t, u(:)
    real(kind=rk), intent(out) :: f(:)
    
    f(1) = u(2)
    f(2) = A * (1.0_rk - u(1)**2)*u(2) - u(1)
    
  end subroutine rhs
  
  subroutine jac (t, u, dfdu, errc)
    real(kind=rk), intent(in) :: t, u(:)
    real(kind=rk), intent(out) :: dfdu(:,:)
    integer, intent(out) :: errc
    dfdu(1,1) = 0.0_rk
    dfdu(1,2) = 1.0_rk
    dfdu(2,1) = - 2.0_rk * A * u(1) * u(2) - 1.0_rk
    dfdu(2,2) = A * (1.0_rk - u(1)**2)
    errc = 0
  end subroutine jac
  
  subroutine user (t, u)
    real(kind=rk), intent(in) :: t, u(:)
    print *, t, u(:)
  end subroutine user
  
end module vanderpol

program vanderpol_test

  use bdf2_kinds
  use vanderpol
  use bdf2_integrator
  implicit none
  
  integer :: nstep, stat
  real(kind=rk) :: tout, t, u(2), atol(2)
  type(bdf2_control) :: control
!  type(bdf2_profile) :: profile
  type(bdf2_state)   :: state
  
  call bdf2_create_state (state, 2)
  atol = (/ 1.0d-6, 1.0d-6 /)
  call bdf2_set_param (control, atol=atol, rtol=1.0d-3, mvec=2, ntol=0.01d0)
  
  t = 0.0d0
  u = (/ 2.0d0, 0.0d0 /)
  
  !profile%msgs = 1
  
  call bdf2_init_state (state, u, t, hstart=1.0d-6, profile=.true.)
  
  nstep = 50000
  tout = 3000.0d0
  
  call bdf2_integrate (state, control, tout=tout, stat=stat, rhs=rhs, user=user)
  select case (stat)
  case (SOLN_AT_TOUT)
    call user (tout, bdf2_interpolate_solution (state, tout))
    t = bdf2_solution_time(state)
    u = bdf2_solution(state)
    call user (t, u)
  case (SOLN_AT_STEP)
    print *, 'Failed to reach final time within', nstep, ' time steps.'
  case default
    print *, 'BDF2_INTEGRATE returned an error condition: stat=', stat
    print *, t
  end select
  
  !print *, state%profile
  
end program vanderpol_test
