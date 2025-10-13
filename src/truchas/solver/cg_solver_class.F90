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
    integer :: num_iter, min_iter, max_iter, print_level
    real(r8) :: abs_tol  ! CG stopping tolerance on the residual error
    real(r8) :: rel_tol  ! CG stopping tolerance on the residual error reduction
    real(r8) :: initial_rnorm, rel_rnorm
  contains
    procedure :: solve
    procedure(matvec), deferred :: ax
    procedure(matvec), deferred :: pc
  end type

  abstract interface
    subroutine matvec(this, x, y)
      import cg_solver, r8
      class(cg_solver), intent(inout) :: this
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: y(:)
    end subroutine
  end interface

contains

  subroutine solve(this, b, x, dim, stat, errmsg)

    use parallel_communication, only: global_dot_product
    use truchas_logging_services

    class(cg_solver), intent(inout) :: this
    real(r8), intent(in) :: b(:) ! Equation RHS
    real(r8), intent(inout) :: x(:) ! Solution (initial guess/final)
    integer, intent(in)  :: dim   ! True dimension of the system
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    real(r8) :: r(size(x)), p(size(x)), q(size(x))
    real(r8) :: rho, rho_init, rho_last, s
    integer :: iter, n, j
    character(80) :: message

    ASSERT(size(x) == size(b))
    ASSERT(dim <= size(x))
    ASSERT(size(x) == 0 .or. dim >= 1)

    n = size(b)

    call this%ax(x, r)
    do j = 1, n
      r(j) = b(j) - r(j)
    end do

    call this%pc(r, p)

    rho = global_dot_product(r(:dim), p(:dim))
    rho_init = rho
    this%initial_rnorm = sqrt(rho_init)
    if (rho_init < 0.0_r8) then
      stat = -2
      errmsg = 'non-symmetric preconditioner'
      return
    end if

    if (rho == 0.0_r8) then ! solved exactly!
      stat = 0
      return
    end if

    do iter = 1, this%max_iter

      call this%ax(p, q)
      s = rho / global_dot_product(p(:dim), q(:dim))
      if (s < 0.0_r8) then
        stat = -1
        errmsg = 'non-symmetric matrix'
        exit
      end if
      do j = 1, n
        x(j) = x(j) + s * p(j)
        r(j) = r(j) - s * q(j)
      end do

      call this%pc(r, q)
      rho_last = rho
      rho = global_dot_product(r(:dim), q(:dim))
      if (rho < 0.0_r8) then
        stat = -2
        errmsg = 'non-symmetric preconditioner'
        exit
      end if
      s = rho / rho_last
      do j = 1, n
        p(j) = q(j) + s * p(j)
      end do

      if (this%print_level > 0) then
        write(message,fmt="(t6,a,i4,2(a,es10.3))") 'cg iterate', iter, &
            ': |r|=', sqrt(rho), ', |r|/|r_prev|=', sqrt(s)
        call TLS_info(message)
      end if

      if (iter >= this%min_iter .and. &
          (rho < this%abs_tol**2 .or. rho < this%rel_tol**2 * rho_init)) then
        stat = 0
        exit
      end if

    end do

    if (iter > this%max_iter) then
      stat = 1
      errmsg = 'failed to converge'
    end if

    this%num_iter = min(iter, this%max_iter)
    this%rel_rnorm = sqrt(rho/rho_init)

  end subroutine solve

end module cg_solver_class
