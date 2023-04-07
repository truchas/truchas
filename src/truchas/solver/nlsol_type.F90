!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module nlsol_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nka_type
  use truchas_logging_services
  implicit none
  private

  type, public :: nlsol
    private
    class(nlsol_model), pointer :: model => null() ! unowned reference
    integer  :: n                   ! number of unknowns
    ! integer  :: seq = -1            ! number of steps taken
    ! real(r8) :: hlast               ! last step size
    ! real(r8) :: hpc                 ! step size built into the current preconditioner
    ! logical  :: usable_pc = .false. ! whether the current preconditioner is usable
    ! integer  :: freeze_count = 0    ! don't increase step size for this number of steps
    integer, public :: itr = 0
    integer  :: mitr            ! maximum number of nonlinear iterations
    real(r8) :: ntol            ! nonlinear solver error tolerance (relative to 1)
    type(nka), allocatable :: nka   ! nonlinear Krylov accelerator
    ! type(state_history)   :: uhist        ! solution history structure

    !! Perfomance counters
    integer :: pcfun_calls = 0      ! number of calls to PCFUN
    ! integer :: updpc_calls = 0      ! number of calls to UPDPC
    ! integer :: updpc_failed = 0     ! number of UPDPC calls returning an error
    ! integer :: retried_bce = 0      ! number of retried BCE steps
    ! integer :: failed_bce = 0       ! number of completely failed BCE steps
    ! integer :: rejected_steps = 0   ! number of steps rejected on error tolerance
    ! real(r8) :: hmin = huge(1.0_r8) ! minimum step size used on a successful step
    ! real(r8) :: hmax = tiny(1.0_r8) ! maximum step size used on a successful step

    !! Diagnostics
    integer :: verbose = 1

    ! !! For state save/restore
    ! type(nlsol), pointer :: cache => null()
  contains
    procedure :: init
    procedure :: solve
  end type nlsol

  type, abstract, public :: nlsol_model
    real(r8), allocatable :: rhs(:) ! made accessible to solvers
  contains
    procedure(model_size), deferred :: size
    procedure(compute_f), deferred :: compute_f
    procedure(apply_precon), deferred :: apply_precon
    procedure(compute_precon), deferred :: compute_precon
    procedure(du_norm), deferred :: du_norm
    procedure(is_converged), deferred :: is_converged
  end type nlsol_model

  abstract interface
    integer function model_size(this)
      import nlsol_model
      class(nlsol_model), intent(in) :: this
    end function
    subroutine compute_f(this, t, u, udot, f, ax)
      import nlsol_model, r8
      class(nlsol_model) :: this
      real(r8), intent(in) :: t
      real(r8), intent(in), contiguous, target :: u(:), udot(:)
      real(r8), intent(out), contiguous, target :: f(:)
      logical, intent(in), optional :: ax
    end subroutine
    subroutine apply_precon(this, t, u, f)
      import nlsol_model, r8
      class(nlsol_model) :: this
      real(r8), intent(in) :: t
      real(r8), intent(in), contiguous, target :: u(:)
      real(r8), intent(inout), contiguous, target :: f(:)
    end subroutine
    subroutine compute_precon(this, t, u, dt)
      import nlsol_model, r8
      class(nlsol_model) :: this
      real(r8), intent(in) :: t, dt
      real(r8), intent(in), contiguous, target :: u(:)
    end subroutine
    real(r8) function du_norm(this, t, u, du)
      import :: nlsol_model, r8
      class(nlsol_model) :: this
      real(r8), intent(in) :: t
      real(r8), intent(in), contiguous, target :: u(:), du(:)
    end function
    logical function is_converged(this, itr, t, u, du, f_lnorm, tol)
      import :: nlsol_model, r8
      class(nlsol_model) :: this
      integer, intent(in) :: itr
      real(r8), intent(in) :: t, tol
      real(r8), intent(in), contiguous, target :: u(:), du(:), f_lnorm(:)
    end function
  end interface

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(nlsol), intent(out) :: this
    class(nlsol_model), intent(in), target :: model
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
  subroutine solve(this, t, h, u0, u, errc)

    class(nlsol), intent(inout) :: this
    real(r8), intent(in)    :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out)   :: errc

    character(256) :: msg
    integer :: itr
    real(r8) :: error, lnormi(3), du(size(u)), udot(size(u))

    errc = 0

    if (allocated(this%nka)) call this%nka%restart
    call this%model%compute_precon(t, u0, h)

    do itr = 1, this%mitr
      !! Evaluate the nonlinear function and precondition it.
      this%pcfun_calls = this%pcfun_calls + 1
      if (h /= 0) then
        udot = (u-u0)/h
      else
        ! If passed zero step size, this is expected to be a system
        ! independent of udot.
        udot = 0
      end if
      call this%model%compute_f(t, u, udot, du)
      lnormi = lnorm(du)
      call this%model%apply_precon(t, u, du)

      !! NKA accelerated correction.
      if (allocated(this%nka)) call this%nka%accel_update(du)

      !! Next solution iterate.
      u = u + du

      !! Error estimate.
      error = this%model%du_norm(t, u, du)
      if (this%verbose >= 2) then
        write(msg,fmt=3) itr, error, lnormi(3)
        call tls_info(trim(msg))
      end if

      !! Check for convergence.
      if (this%model%is_converged(itr, t, u, du, lnormi, this%ntol)) exit
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

end module nlsol_type
