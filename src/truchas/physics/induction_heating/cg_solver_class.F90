!!
!! CG_SOLVER_CLASS
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

module cg_solver_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: cg_solver
    integer :: minitr, maxitr ! min and max number of CG iterations
    integer :: output_level
    real(r8) :: tol  ! CG stopping tolerance on the residual error
    real(r8) :: red  ! CG stopping tolerance on the residual error reduction
    !! Pointers to the system matvec and PC procedures belong here...
  contains
    procedure :: solve
    procedure(matvec), deferred :: ax
    procedure(matvec), deferred :: pc
  end type
  
  abstract interface
    subroutine matvec(this, x, y)
      import cg_solver, r8
      class(cg_solver), intent(in) :: this
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: y(:)
    end subroutine
  end interface

contains

  subroutine solve(this, b, x, dim, stat)
  
    use parallel_communication, only: global_dot_product
    use truchas_logging_services
  
    class(cg_solver), intent(in) :: this
    real(r8), intent(in) :: b(:) ! Equation RHS
    real(r8), intent(inout) :: x(:) ! Solution (initial guess/final)
    integer, intent(in)  :: dim   ! True dimension of the system
    integer, intent(out) :: stat
    
    real(r8) :: r(size(x)), p(size(x)), q(size(x))
    real(r8) :: rho, rho_init, rho_last, s
    integer :: itr, n, j
    character(len=80) :: message
    
    ASSERT(size(x) == size(b))
    ASSERT(dim <= size(x))
    ASSERT(size(x) == 0 .or. dim >= 1)
    
    n = size(b)
    
    !r = b - ax(x)
    call this%ax(x, r)
    do j = 1, n
      r(j) = b(j) - r(j)
    end do
    
    !p = pc(r)
    call this%pc(r, p)
    
    rho = global_dot_product(r(:dim), p(:dim))
    rho_init = rho
    INSIST( rho_init >= 0.0_r8 )
    stat = -1
    
    if (rho == 0.0_r8) then
      stat = 0
      return
    end if
    
    do itr = 1, this%maxitr
    
      !q = ax(p)
      call this%ax(p, q)
      s = rho / global_dot_product(p(:dim), q(:dim))
      INSIST( s >= 0.0_r8 )
      do j = 1, n
        x(j) = x(j) + s * p(j)
        r(j) = r(j) - s * q(j)
      end do
      
      !q = pc(r)
      call this%pc(r, q)
      rho_last = rho
      rho = global_dot_product(r(:dim), q(:dim))
      INSIST( rho >= 0.0_r8 )
      s = rho / rho_last
      do j = 1, n
        p(j) = q(j) + s * p(j)
      end do
      
      if (this%output_level >= 4) then
        write(message,fmt="(t8,a,i4,2(a,es10.3))") 'this iterate', itr, &
            ': |r|=', sqrt(rho), ', |r|/|r_prev|=', sqrt(s)
        call TLS_info (message)
      end if
      
      if (itr >= this%minitr .and. (rho < this%tol**2 .or. rho < this%red**2 * rho_init)) then
        stat = 0
        exit
      end if
      
    end do
    
    if (this%output_level >= 2) then
      write(message,fmt="(t6,a,i4,2(a,es10.3))") 'step this summary:', itr, &
          ' iterations, |r|=', sqrt(rho), ', |r0|=', sqrt(rho_init)
      call TLS_info (message)
    end if
    
  end subroutine solve

end module cg_solver_class
