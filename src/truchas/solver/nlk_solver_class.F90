!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module nlk_solver_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nka_vec_type
  use vector_class
  use truchas_logging_services
  implicit none
  private

  type, abstract, public :: nlk_solver
    private
    integer, public :: num_iter = 0
    integer  :: max_iter
    real(r8) :: abs_tol, rel_tol
    type(nka_vec), allocatable :: nka
    class(vector), allocatable :: r

    !! Diagnostics
    integer :: verbose = 1
  contains
    procedure :: init
    procedure :: solve
    procedure(compute_f), deferred :: compute_f
    procedure(apply_precon), deferred :: apply_precon
    procedure(compute_precon), deferred :: compute_precon
  end type

  abstract interface
    subroutine compute_f(this, u, f)
      import nlk_solver, vector
      class(nlk_solver) :: this
      class(vector), intent(in) :: u
      class(vector), intent(inout) :: f
    end subroutine
    subroutine apply_precon(this, u, f)
      import nlk_solver, vector
      class(nlk_solver) :: this
      class(vector), intent(in) :: u
      class(vector), intent(inout) :: f
    end subroutine
    subroutine compute_precon(this, u)
      import nlk_solver, vector
      class(nlk_solver) :: this
      class(vector), intent(in) :: u
    end subroutine
  end interface

contains

  subroutine init(this, vec, params, stat, errmsg)

    use parameter_list_type

    class(nlk_solver), intent(out) :: this
    class(vector), intent(in) :: vec
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    integer :: maxv
    real(r8) :: vtol

    stat = 0

    context = 'processing ' // params%path() // ': '
    call params%get('nlk-max-iter', this%max_iter, stat, errmsg, default=100)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%max_iter < 2) then
      stat = 1
      errmsg = context//'"nlk-max-iter" must be > 1'
      return
    end if

    call params%get('nlk-abs-tol', this%abs_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%abs_tol < 0.0_r8) then
      stat = 1
      errmsg = context//'"nlk-abs-tol" must be >= 0.0'
      return
    end if

    call params%get('nlk-rel-tol', this%rel_tol, stat, errmsg, default=1.0e-8_r8)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%rel_tol < 0.0_r8) then
      stat = 1
      errmsg = context//'"nlk-rel-tol" must be >= 0.0'
      return
    end if

    call params%get('nlk-max-vec', maxv, stat, errmsg, default=20)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (maxv < 0) then
      stat = 1
      errmsg = context//'"nlk-max-vec" must be >= 0'
      return
    end if
    maxv = min(maxv, this%max_iter-1)

    call params%get('nlk-vec-tol', vtol, stat, errmsg, default=0.01_r8)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (vtol <= 0.0_r8) then
      stat = 1
      errmsg = context//'"nlk-vec-tol" must be > 0.0'
      return
    end if

    call params%get('verbosity', this%verbose, default=1, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if

    call vec%clone(this%r)

    !! Initialize the NKA structure.
    if (maxv > 0) then
      allocate(this%nka)
      call this%nka%init(vec, maxv)
      call this%nka%set_vec_tol(vtol)
    end if

  end subroutine init

  !! Solve the nonlinear system f(u) = 0 using an accelerated fixed point
  !! iteration [1] for the preconditioned system g(u) = PC(f(u)) = 0:
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
  !! of the Newton correction equation  J*du = f(u) where J is an approximation
  !! to the Jacobian of f(u).
  !!
  !! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
  !!     weighted moving finite element code I: in one dimension", SIAM J.
  !!     Sci. Comput;, 19 (1998), pp. 728-765..

  subroutine solve(this, u, stat)

    class(nlk_solver), intent(inout) :: this
    class(vector), intent(inout) :: u
    integer, intent(out) :: stat

    integer :: iter
    real(r8) :: r0norm, rnorm
    character(256) :: msg

    stat = 0

    call this%compute_f(u, this%r)
    r0norm = this%r%norm2()

    if (this%verbose >= 2) then
      write(msg,fmt=3) 0, r0norm
      call tls_info(trim(msg))
    end if

    ! convergence check?

    if (allocated(this%nka)) call this%nka%restart
    call this%compute_precon(u)

    do iter = 1, this%max_iter
      call this%apply_precon(u, this%r)
      if (allocated(this%nka)) call this%nka%accel_update(this%r)
      !u = u + r
      call u%update(1.0_r8, this%r) !FIXME? Why is this "+" -- it ought to be "-" unless precon is doing something weird

      call this%compute_f(u, this%r)
      rnorm = this%r%norm2()

      if (this%verbose >= 2) then
        write(msg,fmt=3) iter, rnorm
        call tls_info(trim(msg))
      end if

      !! Check for convergence.
      if (rnorm <= max(this%abs_tol, this%rel_tol*r0norm)) exit
    end do

    if (iter > this%max_iter) then
      stat = 1 ! too many iterations -- failed to converge
      iter = this%max_iter
    end if
    this%num_iter = iter ! expose the number of iterations performed

    if (this%verbose >= 1) then
      if (stat == 0) then
        write(msg,fmt=2) iter, rnorm
      else
        write(msg,fmt=1) iter, rnorm
      end if
      call tls_info(trim(msg))
    end if

1   format(2x,'nlk_solver: convergence failure: ',i0,' iterations (max), 2-norm =',es10.3)
2   format(2x,'nlk_solver: converged: ',i0,' iterations, 2-norm =',es10.3)
3   format(2x,i6,': 2-norm =',1es10.3)

  end subroutine solve

end module nlk_solver_class
