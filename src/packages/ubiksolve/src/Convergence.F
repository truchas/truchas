! -*- Mode: f90 -*-

module Ubik_Convergence_module
  use ubik_types_module
  use ubik_constants_module
  use ubik_control_type_module
  use ubik_vector_type_module
  use Ubik_Meta_module, only: metaMatVec
  use Ubik_Vector_module, only: Norm, spew
  use Ubik_WriteString_module
  implicit none

  private
  public :: CheckConvergence, &
            StoppingTest

contains

  !*********************************************************************
  !
  ! <para>
  ! Determine whether the iteration has converged and compute true
  ! residual and error estimate based on true residual, if required.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   x - x-vector
  !   b - source vector
  !   MatVec - routine to form product of A and x
  !
  !  Input (optional):
  !   rnorm - residual norm to use in test (if other than 
  !           control%rnorm(control%iter))
  ! 
  !  In/Out:
  !   control - structure containing control parameters
  ! 
  !  Output:
  !   status - return status
  ! 
  ! <SUBROUTINES REQUIRED>
  ! 
  !  WriteString
  !  WriteVector
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  Norm
  !
  !*********************************************************************
  subroutine StoppingTest (x_vec, b, control, MatVec, rnorm_arg)

#include "matvec_proto.F"

    real(real_type), intent(in), dimension(:) :: b
    real(real_type), intent(in), optional :: rnorm_arg
    type(control_type), intent(inout) :: control
    type(vector_type), intent(inout) :: x_vec
    ! </arguments>

    ! Local:
    logical(logical_type) :: stopping_relchg_, stopping_axb_, stopping_b_, &
         stopping_x_, stopping_true_, stopping_r_, stopping_r0_
    logical :: outmode_iterates_
    integer(int_type) :: luout_, norm_type_, sizeofx, status
    real(real_type) :: rnorm_, err_, errt_, xnorm
    real(real_type), allocatable, dimension(:) :: w

    ! Initialize return status.
    call set_status_converged (control)

    ! shorthand

    ! flags related to output mode
    outmode_iterates_ = outmode_iterates(control)

    ! flags related to stopping test
    stopping_relchg_ = stopping_relchg(control)
    stopping_axb_ = stopping_axb(control)
    stopping_b_ = stopping_b(control)
    stopping_r_ = stopping_r(control)
    stopping_r0_ = stopping_r0(control)
    stopping_x_ = stopping_x(control)
    stopping_true_ = stopping_true(control)

    norm_type_ = norm_type(control)
    luout_ = luout(control)

    if (PRESENT(rnorm_arg)) then
       rnorm_ = rnorm_arg
    else
       rnorm_ = rnorm(control)
    end if

    ! Calculate || x || if needed.
    if (stopping_relchg_ .or. stopping_axb_ .or. stopping_x_) then
       xnorm = Norm(norm_type_, x_vec)
    end if

    ! Calculate error estimate based on || b - Ax || as appropriate.
    if (residual_compute_only(control) .or. &
        residual_compute_update(control) .or. &
        stopping_true_) then

       if (rnormt(control) < TINY(0.0)) then

          ! || b - Ax || was not input, so compute.  Note that since this is a
          ! member of control_type we must reset it in control before returning.

          sizeofx = SIZE(x_vec)
          allocate (w(sizeofx), stat=status)
          if (status /= 0) then
             call set_status_alloc_failure (control)
             go to 9999
          end if

          call metaMatVec (MatVec, control, x_vec, w)
          if (matmul_error(control)) then
             go to 9999
          end if
          w = b - w
          call set_rnormt (control, Norm(norm_type_, w))

          ! Write out b-Ax, if appropriate.
          if (outmode_iterates_) then
             call spew (luout_, w, 'StoppingTest: b - Ax (w):')
          end if

          deallocate (w)
       end if

       ! Compute error estimate based on || b - Ax ||.
       if (stopping_relchg_) then

          ! Stopping test is || x - xold || / || x ||, so do not
          ! compute an error estimate based on rnormt.
          errt_ = 0.0
       else if (stopping_axb_) then

          ! Use || b - Ax || / ( || A ||*|| x || + || b || )
          errt_ = rnormt(control) / (anorm(control)*xnorm + bnorm(control))
       else if (stopping_b_) then

          ! Use || b - Ax || / || b ||
          errt_ = rnormt(control) / bnorm(control)
       else if (stopping_x_) then

          ! Use || b - Ax || / || x ||
          errt_ = rnormt(control)
          if (xnorm /= 0.0) then
             errt_ = errt_ / xnorm
          end if
       else if (stopping_r0_) then

          ! Use || b - Ax || / || b - Ax0 ||
          errt_ = rnormt(control) / rnormt(control, 0)
       else if (stopping_r_) then

          ! Use || b - Ax ||
          errt_ = rnormt(control)
       end if
       call set_errt (control, errt_)
    end if

    if (stopping_relchg_) then

       ! Use || x - xold || / || x || as stopping test.
       sizeofx = SIZE(x_vec)
       allocate (w(sizeofx), stat=status)
       if (status /= 0) then
          call set_status_alloc_failure (control)
          go to 9999
       end if
#ifdef dec
       w = vector_minus_vector_type(xold(control), x_vec)
#else
       w = xold(control) - x_vec
#endif
       err_ = Norm(norm_type_, w)
       if (xnorm /= 0.0) then
          err_ = err_/xnorm
       end if
       call set_xold (control, values(x_vec))
       deallocate (w)
    else if (stopping_axb_) then

       ! Use || r || / ( || A ||*|| x || + || b || )
       err_ = rnorm_ / (anorm(control)*xnorm + bnorm(control))
    else if (stopping_b_) then

       ! Use || r || / || b ||
       err_ = rnorm_ / bnorm(control)
    else if (stopping_x_) then

       ! Use || r || / || x ||
       err_ = rnorm_
       if (xnorm /= 0.0) then
          err_ = err_ / xnorm
       end if
    else if (stopping_r0_) then

       ! Use || r || / || b - Ax0 ||
       err_ = rnorm_ / rnorm(control, 0)
    else if (stopping_r_) then

       ! Use || r ||
       err_ = rnorm_
    end if
    call set_err (control, err_)

9999 continue
    return
  end subroutine StoppingTest

  !*********************************************************************
  !
  ! <para>
  ! Determine whether the iteration has converged and compute true
  ! residual and error estimate based on true residual, if required.
  ! </para>
  !
  ! Note that:
  !   - the arrays cpu, bnorm, rnorm, err, rnormt, and errt are members of
  !     control_type
  !   - rnormt(0) and rnorm(iter) are input
  !   - rnormt(iter) can either be input or output
  !   - err(iter) is output
  !   - errt(iter) is conditionally output
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   x - x-vector
  !   b - source vector
  !   MatVec - routine to form product of A and x
  !
  !  Input (optional):
  !   rnorm - residual norm to use in test (if other than 
  !           control%rnorm(control%iter))
  ! 
  !  In/Out:
  !   control - structure containing control parameters
  ! 
  !  Output:
  !   status - return status
  !
  ! <SUBROUTINES REQUIRED>
  ! 
  !  StoppingTest
  !  WriteString
  !  WriteVector
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  Norm
  !
  !*********************************************************************
  subroutine CheckConvergence (x_vec, b, control, MatVec, rnorm_arg)

#include "matvec_proto.F"

    real(real_type), intent(in), dimension(:) :: b
    real(real_type), intent(in), optional :: rnorm_arg
    type(control_type), intent(inout) :: control
    type(vector_type), intent(inout) :: x_vec
    ! </arguments>

    ! Local:
    logical(logical_type) :: stopping_relchg_, stopping_true_, tty_update
    logical(logical_type) :: outmode_errors_, outmode_system_, outmode_iterates_, outmode_summary_
    character(len=132) :: out_string
    integer(int_type) :: luerr_, luout_, lutty_, status
    real(real_type) :: err_test, rnorm_, err_, rnormt_, errt_
    real(real_type) :: relchg, rnorm_old

    ! shorthand
    outmode_errors_ = outmode_errors(control)
    outmode_system_ = outmode_system(control)
    outmode_iterates_ = outmode_iterates(control)
    outmode_summary_ = outmode_summary(control)
    stopping_relchg_ = stopping_relchg(control)
    stopping_true_ = stopping_true(control)
    luerr_ = luerr(control)
    luout_ = luout(control)
    lutty_ = lutty(control)

    if (PRESENT(rnorm_arg)) then
       rnorm_ = rnorm_arg
    else
       rnorm_ = rnorm(control)
    end if

    ! Write out current iterate, if appropriate.
    if (outmode_iterates_) then
       call spew (luout_, x_vec, 'CheckConvergence: CURRENT ITERATE (x):')
    end if

    ! Compute residual norm(s) and error estimate(s) as appropriate.
    call StoppingTest (x_vec, b, control, MatVec, rnorm_)

    ! shorthand
    err_ = err(control)
    rnormt_ = rnormt(control)
    errt_ = errt(control)

    ! write one-line status update to tty
    tty_update = .false.
    if (output_frequency(control) /= 0) then
       if (MOD(iter(control),output_frequency(control)) == 0) then
          tty_update = .true.
       end if
    end if
    if (tty_update) then
       if (iter(control) > 0) then
          rnorm_old = rnorm(control, iter(control)-1)
       else
          rnorm_old = 0.0
       end if
       if (rnorm_old > smallest(control)) then
          relchg = ABS( rnorm_ - rnorm_old ) / rnorm_old
       else
          relchg = 0.0
       end if
       write(out_string,100) iter(control), rnorm_, relchg, err_
       call spew (lutty_, out_string)
    end if
100 format (1p, 5x, &
         'iter =', i5, &
         ',  ||r|| =', e9.2, &
         ',  relative change in ||r|| =', e9.2, &
         ',  error estimate =', e9.2)

    ! write one-line status update to out
    if (outmode_summary_) then
       if (residual_compute_only(control) .or. &
           residual_compute_update(control) .or. &
           stopping_true_) then
          if (iter(control) == 0 .or. outmode_iterates_) then
             write(out_string,112)
             call spew (luout_, out_string)
          end if
          write(out_string,110) iter(control), rnorm_, err_, rnormt_, errt_
          call spew (luout_, out_string)
       else
          if (iter(control) == 0 .or. outmode_iterates_) then
             write(out_string,111)
             call spew (luout_, out_string)
          end if
          write(out_string,110) iter(control), rnorm_, err_
          call spew (luout_, out_string)
       end if
    end if
110 format (5x, i6, 2x, 2(3x,e12.5), 2x, 2(3x,e12.5))
111 format (t8, 'iter', t21, '||r||', t37, 'err')
112 format (t8, 'iter', t21, '||r||', t37, 'err', t52, '||b-Ax||', t68, 'errt')

    ! Test for convergence using appropriate estimate of error.
    ! If converged, set status, write final solution estimate if
    ! desired, and return.  If not, increment iteration counter 
    ! and check to see whether itmax has been exceeded.
    if (stopping_true_) then
       err_test = errt_
    else
       err_test = err_
    end if

    if (err_test <= eps(control)) then

       ! converged
       if (outmode_system_) then
          call spew (luout_, x_vec, 'CheckConvergence: AFTER CONVERGENCE (x):')
       end if
       call set_status_converged (control)
    else

       ! convergence not achieved, but no problems either, so continue (unless
       ! max number of iters has been exceeded)
       call set_status_next_iter (control)

       ! increment_iter will set status to indicate itmax exceeded if appropriate
       call increment_iter (control)

       if (itmax_exceeded(control)) then
          if (outmode_errors_) then
             call spew (luerr_, ' ')
             call spew (luerr_,    '*************************************************************************')
             call spew (luerr_,    'CheckConvergence:')
             call spew (luerr_,    '  iteration failed to converge in itmax iterations')
             write(out_string,905)            '  itmax = ', itmax(control)
             call spew (luerr_, out_string)
             write(out_string,940)            '  convergence criterion = ', eps(control)
             call spew (luerr_, out_string)
             call spew (luerr_,    '  parameters at termination:')
             if (.not.stopping_relchg_) then
                write(out_string,940)         '    norm of recursively-computed residual = ', rnorm_
                call spew (luerr_, out_string)
                write(out_string,940)         '    error estimate using recursively-computed residual = ', err_
                call spew (luerr_, out_string)
             else
                write(out_string,940)         '    relative change in successive iterations = ', err_
                call spew (luerr_, out_string)
             end if
             if (residual_compute_only(control) .or. &
                 residual_compute_update(control) .or. &
                 stopping_true_) then
                write(out_string,940)         '    norm of b-Ax = ', rnormt_
                call spew (luerr_, out_string)
                write(out_string,940)         '    error estimate using b-Ax = ', errt_
                call spew (luerr_, out_string)
             end if
             if (stopping_relchg_) then
                call spew (luerr_, '  convergence based on relative change in successive iterations')
             else if (stopping_true_) then
                call spew (luerr_, '  convergence based on error estimate using b-Ax')
             else
                call spew (luerr_, '  convergence based on error estimate using recursively-computed residual')
             end if
             call spew (luerr_,    '*************************************************************************')
          end if
       end if
    end if

9999 continue
    call set_cpu (control)

#include "formats.F"
    return
  end subroutine CheckConvergence

end module Ubik_Convergence_module
