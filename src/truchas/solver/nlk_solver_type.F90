!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module nlk_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nka_type
  use truchas_logging_services
  implicit none
  private

  type, public :: nlk_solver
    private
    class(nlk_solver_model), pointer :: model => null() ! unowned reference
    integer  :: n                   ! number of unknowns
    integer, public :: itr = 0
    integer  :: mitr            ! maximum number of nonlinear iterations
    real(r8) :: ntol            ! nonlinear solver error tolerance (relative to 1)
    type(nka), allocatable :: nka   ! nonlinear Krylov accelerator
    ! type(state_history)   :: uhist        ! solution history structure

    !! Perfomance counters
    integer :: pcfun_calls = 0      ! number of calls to PCFUN

    !! Diagnostics
    integer :: verbose = 1
  contains
    procedure :: init
    procedure :: solve
  end type nlk_solver

  type, abstract, public :: nlk_solver_model
    real(r8), allocatable :: rhs(:) ! made accessible to solvers
  contains
    procedure(model_size), deferred :: size
    procedure(compute_f), deferred :: compute_f
    procedure(apply_precon), deferred :: apply_precon
    procedure(compute_precon), deferred :: compute_precon
    procedure(du_norm), deferred :: du_norm
    procedure(is_converged), deferred :: is_converged
  end type nlk_solver_model

  abstract interface
    integer function model_size(this)
      import nlk_solver_model
      class(nlk_solver_model), intent(in) :: this
    end function
    subroutine compute_f(this, u, f, ax)
      import nlk_solver_model, r8
      class(nlk_solver_model) :: this
      real(r8), intent(in) :: u(:)
      real(r8), intent(out) :: f(:)
      logical, intent(in), optional :: ax
    end subroutine
    subroutine apply_precon(this, u, f)
      import nlk_solver_model, r8
      class(nlk_solver_model) :: this
      real(r8), intent(in) :: u(:)
      real(r8), intent(inout) :: f(:)
    end subroutine
    subroutine compute_precon(this, u)
      import nlk_solver_model, r8
      class(nlk_solver_model) :: this
      real(r8), intent(in) :: u(:)
    end subroutine
    real(r8) function du_norm(this, u, du)
      import :: nlk_solver_model, r8
      class(nlk_solver_model) :: this
      real(r8), intent(in) :: u(:), du(:)
    end function
    logical function is_converged(this, itr, u, du, f_lnorm, tol)
      import :: nlk_solver_model, r8
      class(nlk_solver_model) :: this
      integer, intent(in) :: itr
      real(r8), intent(in) :: tol
      real(r8), intent(in) :: u(:), du(:), f_lnorm(:)
    end function
  end interface

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(nlk_solver), intent(out) :: this
    class(nlk_solver_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    integer :: maxv
    real(r8) :: vtol
    procedure(pardp), pointer :: dp

    stat = 0
    this%model => model

    this%n = model%size()
    INSIST(this%n > 0)

    context = 'processing ' // params%path() // ': '
    call params%get('nlk-max-iter', this%mitr, stat, errmsg, default=100)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%mitr < 2) then
      stat = 1
      errmsg = context//'"nlk-max-iter" must be > 1'
      return
    end if

    call params%get('nlk-tol', this%ntol, stat, errmsg, default=1e-12_r8)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%ntol <= 0.0_r8 .or. this%ntol > 1.0_r8) then
      stat = 1
      errmsg = context//'"nlk-tol" must be > 0.0 and <= 1.0'
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
    maxv = min(maxv, this%mitr-1)

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

    !! Initialize the NKA structure.
    if (maxv > 0) then
      allocate(this%nka)
      call this%nka%init(this%n, maxv)
      call this%nka%set_vec_tol(vtol)
      dp => pardp ! NB: in F2008 can make pardp an internal sub and pass directly
      call this%nka%set_dot_prod(dp)
    end if

  end subroutine init

  function pardp(a, b) result(dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    ASSERT(size(a) == size(b))
    dp = global_dot_product(a, b)
  end function


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

  subroutine solve(this, u, errc)

    class(nlk_solver), intent(inout) :: this
    real(r8), intent(inout) :: u(:)
    integer,  intent(out)   :: errc

    character(256) :: msg
    integer :: itr
    real(r8) :: error, lnormi(3), du(size(u))

    errc = 0

    if (allocated(this%nka)) call this%nka%restart
    call this%model%compute_precon(u)

    do itr = 1, this%mitr
      call this%model%compute_f(u, du)
      lnormi = lnorm(du)
      call this%model%apply_precon(u, du) !TODO: why u as input?
      if (allocated(this%nka)) call this%nka%accel_update(du)
      u = u + du

      !! Error estimate.
      error = this%model%du_norm(u, du)
      if (this%verbose >= 2) then
        write(msg,fmt=3) itr, error, lnormi(3)
        call tls_info(trim(msg))
      end if

      !! Check for convergence.
      if (this%model%is_converged(itr, u, du, lnormi, this%ntol)) exit
    end do

    this%itr = itr ! expose the number of iterations performed

    if (itr > this%mitr) errc = 1 ! too many iterations -- failed to converge

    if (this%verbose >= 1) then
      if (errc == 0) then
        write(msg,fmt=2) itr, error, lnormi(3)
      else
        write(msg,fmt=1) itr, error, lnormi(3)
      end if
      call tls_info(trim(msg))
    end if

1   format(2x,'NLK BCE solve FAILED: ',i6,' iterations (max), error,lnorm_inf = ',2es13.3)
2   format(2x,'NLK BCE solve succeeded: ',i6,' iterations, error,lnorm_inf = ',2es13.3)
3   format(2x,i6,': error,lnorm_inf = ',2es13.3)

  end subroutine solve


  function lnorm(u)
    use parallel_communication, only: global_sum, global_dot_product, global_maxval
    real(r8), intent(in) :: u(:)
    real(r8) :: lnorm(3)
    lnorm(1) = global_sum(abs(u))
    lnorm(2) = sqrt(global_dot_product(u,u))
    lnorm(3) = global_maxval(abs(u))
  end function lnorm

end module nlk_solver_type
