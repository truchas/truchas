!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module nodal_FEM_1D

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private
  
  public :: init_simple
  public :: eval_f, update_pc, user, enorm, udot
  
  public :: x
  
  integer :: n, prob = 0
  real(r8) :: u_left, u_right
  real(r8), allocatable :: x(:), dx(:), jac(:,:)
  
  real(r8), public :: atol, rtol
  
  !! Parameters for problem 1.
  real(r8), public :: d
  
contains

  subroutine init_simple (nnode, x0, x1, u0, u1, pnum)
  
    integer,  intent(in) :: nnode   ! Number of mesh nodes.
    real(r8), intent(in) :: x0, x1  ! End points of the domain interval.
    real(r8), intent(in) :: u0, u1  ! Dirichlet boundary values.
    integer,  intent(in) :: pnum    ! Problem number (for PDE parameters).
    
    integer :: j
    real(r8) :: h
    
    ASSERT( nnode > 2 )
    ASSERT( x0 < x1 )
    
    n = nnode
    prob = pnum  ! set the problem number.
    
    u_left = u0
    u_right = u1
    
    allocate(x(n), dx(n-1), jac(3,n))
    
    !! Create the equi-spaced mesh.
    h = (x1 - x0) / real(n-1,kind=r8)
    dx = h
    x(1) = x0
    do j = 2, n-1
      x(j) = x(j-1) + h
    end do
    x(n) = x1
    
  end subroutine init_simple

  subroutine eval_diff_coef (u, a)
  
    real(r8), intent(in)  :: u(:)
    real(r8), intent(out) :: a(:)
    
    integer :: j
    
    ASSERT( size(u) == n )
    ASSERT( size(a) == n-1 )
  
    select case (prob)
    case (1)  ! constant coefficient
      a = d
    case (2)
      do j = 1, n-1
        a(j) = d + max(0.0_r8, 0.5_r8*(u(j) + u(j+1)))
      end do
    case default
      INSIST( .false. )
    end select

  end subroutine eval_diff_coef

  subroutine eval_mass_matrix (dx, m)
  
    real(r8), intent(in)  :: dx(:)
    real(r8), intent(out) :: m(:,:)
    
    integer :: j
    
    ASSERT( size(m,1) == 3 )
    ASSERT( size(m,2) == size(dx)+1 )
    
    m(2,1) = 0.0_r8
    do j = 1, size(dx)
      m(2,j) = m(2,j) + dx(j)/3.0_r8
      m(3,j) = dx(j)/6.0_r8
      m(1,j+1) = dx(j)/6.0_r8
      m(2,j+1) = dx(j)/3.0_r8
    end do
    
  end subroutine eval_mass_matrix
  
  subroutine tdfactor (a)
  
    real(r8), intent(inout) :: a(:,:)
    
    integer :: j
    
    ASSERT( size(a,1) == 3 )
    
    do j = 2, size(a,2)
      a(3,j-1) = a(3,j-1)/a(2,j-1)
      a(2,j) = a(2,j) - a(1,j)*a(3,j-1)
    end do
    
  end subroutine tdfactor

  subroutine tdsolve (a, x)
  
    real(r8), intent(in) :: a(:,:)
    real(r8), intent(inout) :: x(:)
    
    integer :: j
    
    ASSERT( size(a,1) == 3 )
    ASSERT( size(a,2) == size(x) )
    
    !! Forward substitution.
    x(1) = x(1)/a(2,1)
    do j = 2, size(x)
      x(j) = (x(j) - a(1,j)*x(j-1))/a(2,j)
    end do
    
    !! Backward substitution.
    do j = size(x)-1, 1, -1
      x(j) = x(j) - a(3,j)*x(j+1)
    end do
    
  end subroutine tdsolve

  subroutine eval_f (t, u, udot, f)
  
    real(r8), intent(in)  :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    
    integer :: j
    real(r8) :: a(n-1)
    
    ASSERT( size(u) == n )
    ASSERT( size(udot) == n )
    ASSERT( size(f) == n )
    
    call eval_diff_coef (u, a)
    f(2) = (dx(1)/3.0_r8)*udot(2) + (a(1)/dx(1))*(u(2) - u_left)
    do j = 2, n-2
      f(j)   = f(j) + (dx(j)/6.0_r8)*(2.0_r8*udot(j) + udot(j+1)) - (a(j)/dx(j))*(u(j+1) - u(j))
      f(j+1) =        (dx(j)/6.0_r8)*(udot(j) + 2.0_r8*udot(j+1)) + (a(j)/dx(j))*(u(j+1) - u(j))
    end do
    f(n-1) = f(n-1) + (dx(n-1)/3.0_r8)*udot(n-1) - (a(n-1)/dx(n-1))*(u_right - u(n-1))
    
    f(1) = u(1) - u_left
    f(n) = u(n) - u_right
    
    !! Precondition the function.
    call tdsolve (jac, f)
    
  end subroutine eval_f
  
  function udot (t, u)
  
    real(r8), intent(in)  :: t, u(:)
    real(r8) :: udot(n)
    
    integer :: j
    real(r8) :: m(3,n), a(n-1)
    
    ASSERT( size(u) == n )
    
    !! Evaluate the rhs of the linear system into UDOT.
    call eval_diff_coef (u, a)
    udot(2) = - (a(1)/dx(1))*(u(2) - u_left)
    do j = 2, n-2
      udot(j)   = udot(j) + (a(j)/dx(j))*(u(j+1) - u(j))
      udot(j+1) =         - (a(j)/dx(j))*(u(j+1) - u(j))
    end do
    udot(n-1) = udot(n-1) + (a(n-1)/dx(n-1))*(u_right - u(n-1))
    
    !! Evaluate the mass matrix.
    call eval_mass_matrix (dx, m)
    
    !! Solve for UDOT on the interior nodes only.
    call tdfactor (m(:,2:n-1))
    call tdsolve (m(:,2:n-1), udot(2:n-1))
    
    !! Time independent Dirichlet BV.
    udot(1) = 0.0_r8
    udot(n) = 0.0_r8
      
  end function udot
  
  subroutine update_pc (t, u, h, errc)
  
    real(r8), intent(in) :: t, u(:), h
    integer, intent(out) :: errc
    
    integer :: j
    real(r8) :: a(n-1), tmp
    
    !! Jacobian of the linear term in udot.
    call eval_mass_matrix (dx, jac)
    jac = jac / h
    
    !! Jacobian of the linear diffusion term in u.
    call eval_diff_coef (u, a)
    do j = 1, n-1
      tmp = a(j)/dx(j)
      jac(2,j) = jac(2,j) + tmp
      jac(3,j) = jac(3,j) - tmp
      jac(1,j+1) = jac(1,j+1) - tmp
      jac(2,j+1) = jac(2,j+1) + tmp
    end do
    
    !! Dirichlet BC at the left-most node.
    jac(2,1) = 1.0_r8
    jac(3,1) = 0.0_r8
    jac(1,2) = 0.0_r8
    
    !! Dirichlet BC at the right-most node.
    jac(2,n) = 1.0_r8
    jac(1,n) = 0.0_r8
    jac(3,n-1) = 0.0_r8
    
    call tdfactor (jac) 

    errc = 0
    
  end subroutine update_pc
  
  subroutine user (t, u)
    real(r8), intent(in) :: t, u(:)
    integer :: j
    print *, 'T=', t
    print '(2es13.5)', (x(j), u(j), j=1,n)
  end subroutine user
  
  function enorm (u, du)
    real(r8),      intent(in) :: u(:), du(:)
    real(r8) :: enorm
    enorm = maxval(abs(du)/(atol + rtol*abs(u)))
  end function enorm
  
end module nodal_FEM_1D


program main

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nodal_FEM_1D
  use bdf2_dae
  implicit none
  
  integer :: nstep, status, nnode
  real(r8) :: tout, t, h
  real(r8), allocatable :: u(:)
  type(state) :: s
  
  nnode = 201
  allocate(u(nnode))
  
  call init_simple (nnode, x0=0.0_r8, x1=1.0_r8, u0=0.0_r8, u1=0.0_r8, pnum=2)
  d = 0.0002_r8
  
  
  rtol = 0.0
  atol = 1.0d-5
  
  call create_state (s, size(u), mvec=2, ntol=0.01d0)
  
  u = sin(4.0_r8*atan(1.0_r8)*x)
  t = 0.0d0
  
  call set_initial_state (s, t, u, udot(t,u))
  
  nstep = 10
  tout = 0.2d0
  h = 1.0d-5
  
  call bdf2_step_driver (s, eval_f, update_pc, enorm, h, status, tout=tout)
  select case (status)
  case (SOLVED_TO_TOUT)
    call user (tout, interpolated_solution (s, tout))
  case default
    print *, 'BDF2_STEP_DRIVER returned an unknown status: ', status
    print *, t
  end select
  
  call write_bdf2_stepping_statistics (s, 6)
  
end program main
