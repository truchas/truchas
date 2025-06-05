!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module cgsol_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  use nlsol_type, only: nlsol_model
  use truchas_logging_services
  implicit none
  private

  public :: nlsol_model ! re-export

  type, public :: cgsol
    private
    integer, public :: itr = 0
    class(nlsol_model), pointer :: model => null() ! unowned reference
    integer :: verbose
    integer :: mitr ! maximum number of iterations
    real(r8) :: tol ! CG stopping tolerance on the residual error
    real(r8) :: red ! CG stopping tolerance on the residual error reduction
    real(r8), allocatable :: r(:), p(:), q(:), udot(:)
  contains
    procedure :: init
    procedure :: solve
  end type cgsol

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(cgsol), intent(out) :: this
    class(nlsol_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context

    stat = 0
    this%model => model
    INSIST(model%size() > 0)

    context = 'processing ' // params%path() // ': '
    call params%get('nlk-max-iter', this%mitr, default=100, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%mitr < 2) then
      stat = 1
      errmsg = context//'"nlk-max-iter" must be > 1'
      return
    end if

    call params%get('nlk-tol', this%tol, default=1e-12_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%tol < 0 .or. this%tol > 1) then
      stat = 1
      errmsg = context//'"nlk-tol" must be >= 0.0 and <= 1.0'
      return
    end if

    call params%get('red-tol', this%red, default=1e-8_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%red <= 0 .or. this%red > 1) then
      stat = 1
      errmsg = context//'"red-tol" must be > 0.0 and <= 1.0'
      return
    end if

    call params%get('verbosity', this%verbose, default=1, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if

    allocate(this%r(model%size()), this%p(model%size()), this%q(model%size()), &
        this%udot(model%size()))
    this%udot = 0

  end subroutine init


  subroutine solve(this, t, h, u0, u, errc)

    use parallel_communication, only: global_dot_product

    class(cgsol), intent(inout) :: this
    real(r8), intent(in) :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out) :: errc

    character(256) :: msg
    integer :: itr
    real(r8) :: error, lnormi(3), rho, rho_cutoff, s

    ASSERT(size(u) == size(this%r))

    errc = 0
    call this%model%compute_precon(t, u, h)
    call this%model%compute_f(t, u, this%udot, this%r)
    this%p = this%r
    call this%model%apply_precon(t, u, this%p)
    rho = global_dot_product(this%r, this%p)
    rho_cutoff = max(this%tol**2, this%red**2 * rho)
    print *, rho, rho_cutoff
    !INSIST(rho >= 0)

    do itr = 1, this%mitr
      call this%model%compute_f(t, this%p, this%udot, this%q, ax=.true.)
      s = rho / global_dot_product(this%p, this%q)
      !print *, "s: ", s
      u = u + s * this%p
      this%r = this%r - s * this%q
      !this%r = this%r + s * this%q
      lnormi = lnorm(this%r)
      !INSIST(s >= 0)

      this%q = this%r
      call this%model%apply_precon(t, u, this%q)
      s = rho
      rho = global_dot_product(this%r, this%q)
      s = rho / s
      this%p = this%q + s * this%p
      ASSERT(.not.ieee_is_nan(rho))
      ASSERT(ieee_is_finite(rho))
      ASSERT(ieee_is_finite(s))
      print *, "r: ", rho, rho_cutoff
      ! print *
      !INSIST(rho >= 0)

      ! if (this%verbose >= 2) then
      !   error = this%model%du_norm(t, u, this%du)
      !   write(msg,fmt=3) itr, error, lnormi(3)
      !   call tls_info(trim(msg))
      ! end if
      ! if (this%model%is_converged(itr, t, u, this%du, lnormi, this%ntol)) exit
      if (abs(rho) < rho_cutoff) exit
    end do

    this%itr = itr ! expose the number of iterations performed

    if (itr > this%mitr) errc = 1 ! too many iterations -- failed to converge

    if (this%verbose >= 1) then
      !if (this%verbose < 2) error = this%model%du_norm(t, u, this%du)
      if (errc == 0) then
        write(msg,fmt=2) itr, sqrt(rho), lnormi(3)
      else
        write(msg,fmt=1) itr, sqrt(rho), lnormi(3)
      end if
      call tls_info(trim(msg))
    end if

1   format(2x,'CG solve FAILED: ',i6,' iterations (max), error, lnorm_inf = ',2es13.3)
2   format(2x,'CG solve succeeded: ',i6,' iterations, error, lnorm_inf = ',2es13.3)
3   format(2x,i6,': error, lnorm_inf = ',2es13.3)

  end subroutine solve


  function lnorm(u)
    use parallel_communication, only: global_sum, global_dot_product, global_maxval
    real(r8), intent(in) :: u(:)
    real(r8) :: lnorm(3)
    lnorm(1) = global_sum(abs(u))
    lnorm(2) = sqrt(global_dot_product(u,u))
    lnorm(3) = global_maxval(abs(u))
  end function lnorm

end module cgsol_type
