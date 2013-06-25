! -*- Mode: f90 -*-

!**********************************************************************
! $Id: BCGS.F,v 1.30 2002/12/12 01:12:52 turner Exp $
! 
! <para>
! Compute the solution to a nonsymmetric system of linear equations
! of the form Ax=b by bi-conjugate gradients squared (BCGS; also
! known as conjugate gradients squared) or the Freund transpose-free
! quasi-minimal residuals method (TFQMR).
! 
! An error message is written out and control is returned to the calling
! routine if a breakdown in iterates is encountered or if the maximum
! number of iterations is exceeded.  At that time, err(control%iter) can
! be examined to determine if execution can continue.
! </para>
! 
! <PARAMETER LIST>
! 
!  Input:
!   b - source vector
!   MatVec - routine to form product of coefficient and vector
!
!  Input (optional):
!   ApplyPrecond - routine to perform preconditioning
! 
!  In/Out:
!   control - structure containing control parameters
!   x - solution vector (whatever is in x on entry is used as an
!       initial guess)
! 
! <SUBROUTINES REQUIRED>
! 
!  CheckConvergence
!  Swap
!  WriteIt
! 
! <FUNCTIONS REQUIRED>
! 
!  Dot_Product
!  Norm
!
!**********************************************************************
module Ubik_BCGS_module
  implicit none

  private
  public :: BCGS

contains

  recursive subroutine BCGS (x, b, control, MatVec, ApplyPrecond)
    use ubik_types_module
    use ubik_constants_module
    use ubik_control_type_module
    use ubik_vector_type_module
    use Ubik_Convergence_module, only: CheckConvergence
    use Ubik_Meta_module
    use Ubik_ParallelInfo_module, only: initialize_scope
    use Ubik_Vector_module, only: Ubik_Dot_Product, &
                             Norm, &
                             Swap, &
                             spew
    use Ubik_WriteString_module

#include "matvec_proto.F"
#include "precond_proto.F"

    optional :: ApplyPrecond
    real(real_type), intent(in), dimension(:) :: b
    type(control_type), intent(inout) :: control
    real(real_type), intent(inout), dimension(:), target :: x
    ! </arguments>

    ! Local:
    logical(logical_type) :: outmode_full_, precond
    character(len=132) :: out_string
    integer(int_type) :: luout_, luerr_, norm_type_, nunk, status
    real(real_type) :: &
         alpha, beta, rho, rhoold, sigma, eta, rnorm_qmr, tau, term, &
         theta, theta2, rnorm_iter
    real(real_type), allocatable, dimension(:) :: d, h, q, u, work
    real(real_type), allocatable, dimension(:), target :: g, v, r, w
    type(vector_type) :: g_vec, v_vec, x_vec, r_vec, w_vec

    ! Reset return status.
    call set_status_converged (control)

    ! local alias for number of unknowns
    nunk = SIZE(x)

    ! check validity of stopping test - modify if necessary
    call check_stopping_nonstationary (control)

    ! initialize control params
    call initialize (control)

    ! by default use Ubik_BCGS
    if (method(control) /= method_BCGS .and. method(control) /= method_TFQMR) then
       call set_method (control, method_TFQMR)
    end if

    ! shorthand
    luerr_ = luerr(control)
    luout_ = luout(control)
    norm_type_ = norm_type(control)
    outmode_full_ = outmode_full(control)

    ! check to see if b is essentially zero
    !
    ! this is dangerous since we're not scaling by anything, but
    ! the default value of smallest(control) is TINY(1.0), or
    ! typically about 1e-38, so the hope is that if the 1-norm of
    ! b (the sum of abs. values of the elements of the rhs) is less
    ! than that, it really is essentially zero, so returning the
    ! trivial solution is reasonable
    call set_bnorm (control, Norm(1, b))
    if (bnorm(control) < smallest(control)) then

       ! b is essentially zero, so declare victory
       x = 0.0

       if (outmode_warnings(control)) then
          call Spew (luerr_, ' ')
          call Spew (luerr_, '*****************************************')
          call Spew (luerr_, 'BCGS: WARNING')
          call Spew (luerr_, '  it appears that the right-hand side (rhs) of the system')
          call Spew (luerr_, '  is essentially zero:')
          call Spew (luerr_, ' ')
          write(out_string,*) '    ||b||_1 =', bnorm(control)
          call Spew (luerr_, out_string)
          call Spew (luerr_, ' ')
          call Spew (luerr_, '  where ||b||_1 is the 1-norm of the rhs, or the sum of')
          call Spew (luerr_, '  the absolute values of the elements of the rhs')
          call Spew (luerr_, ' ')
          call Spew (luerr_, '  if you believe the system should still be valid, try')
          call Spew (luerr_, '  scaling the system before attempting to solve')
          call Spew (luerr_, '*****************************************')
       end if
       go to 999
    else if (need_bnorm(control) .and. norm_type_ /= 1) then

       ! re-compute || b || if using Linf or L2 norm
       call set_bnorm (control, Norm(norm_type_, b))
    end if

    ! flag indicating whether preconditioning
    precond = PRESENT(ApplyPrecond)

    ! Allocate local arrays.
    allocate (d(nunk), h(nunk), q(nunk), r(nunk), u(nunk), &
              g(nunk), v(nunk), w(nunk), work(nunk), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if

    ! Initialize local arrays.
    d = 0.0
    h = 0.0
    q = 0.0
    r = 0.0
    u = 0.0
    g = 0.0
    v = 0.0
    w = 0.0
    work = 0.0

    ! Initialize vector derived types.
    call set_values_ptr (g_vec, g)
    call set_values_ptr (v_vec, v)
    call set_values_ptr (x_vec, x)
    call set_values_ptr (r_vec, r)
    call set_values_ptr (w_vec, w)

    ! Write initial system if desired.
    if (outmode_system(control)) then
       call spew (luout_, x, b, 'BCGS: ON ENTRY (x,b):')
    end if

    ! Compute || b || if needed.
    if (need_bnorm(control)) then
       call set_bnorm (control, Norm(norm_type_, b))
       if (bnorm(control) < smallest(control)) then
          call set_status_invalid_input (control)
          if (outmode_errors(control)) then
             call spew (luerr_, ' ')
             call spew (luerr_, '*******************************')
             call spew (luerr_, 'BCGS: ERROR')
             call spew (luerr_, '  stopping test requires norm of right-hand side, but')
             call spew (luerr_, '  apparently the right-hand side is effectively zero')
             call spew (luerr_, '*******************************')
          end if
          go to 888
       end if
    end if

    ! Initialize xold if needed.
    if (need_xold(control)) then
       call set_xold (control, x)
    end if

    ! Compute initial residual (w=b-Ax) and write out if desired.
    call metaMatVec (MatVec, control, x_vec, w)
    if (matmul_error(control)) then
       go to 888
    end if

    w = b - w
    if (outmode_full_) then
       call spew (luout_, w, 'BCGS: AFTER CALCULATING INITIAL RESIDUAL (w):')
    end if

    ! Compute || r0 ||.  Return if zero.
    call set_rnorm (control, Norm(norm_type_, w))
    call set_rnormt (control, rnorm(control))

    ! Check to see if initial guess meets the convergence criterion.
    if (rnorm(control) <= smallest(control)) then
       call set_status_converged (control)
       go to 888
    end if

    ! Initial convergence check.  Note that on entry control%iter = 0 and
    ! CheckConvergence increments control%iter.
    call CheckConvergence (x_vec, b, control, MatVec)
    if (.not.next_iter(control)) then
       go to 888
    end if

    ! Initialize TFQMR parameters if necessary.
    theta = 0.0
    theta2 = 0.0
    eta = 0.0
    select case (method(control))
    case (method_TFQMR)
       tau = rnorm(control, iteration=0)
       d = 0.0
    case (method_BCGS)
       tau = 0.0
    case default
       call set_status_invalid_input (control)
       go to 888
    end select

    ! Initialize r.
    if (precond) then

       ! Call user routine to perform preconditioning.
       call metaPrecond (ApplyPrecond, control, w, r_vec)
       if (precond_error(control)) then
          go to 888
       end if
    else
       r = w
    end if

    ! Initialize q.
    q = r

    ! Compute rhoold.
    rhoold = Ubik_Dot_Product(r, q) ! rhoold = (r,q)

    ! First part of first iteration is special since beta is zero
    ! and there is no h vector yet.
    u = r
    g = r
    if (outmode_full_) then
       call spew (luout_, r, 'BCGS: BEFORE START OF MAIN LOOP (r):')
    end if

    ! ************* START MAIN LOOP *************

    do
       if (outmode_full_) then
          call spew (luout_, ' ')
          write(out_string,905) '*** ITERATION NUMBER: ', iter(control)
          call spew (luout_, out_string)
          call spew (luout_, u, 'BCGS: AFTER START OF MAIN LOOP (u):')
          call spew (luout_, g, 'BCGS: AFTER START OF MAIN LOOP (g):')
       end if

       ! Call user routine to compute matrix-vector product.
       call metaMatVec (MatVec, control, g_vec, work)
       if (matmul_error(control)) then
          exit
       end if

       if (precond) then

          ! Call user routine to perform preconditioning.
          call metaPrecond (ApplyPrecond, control, work, w_vec)
          if (precond_error(control)) then
             exit
          end if
       else
          work = w
       end if

       if (outmode_full_) then
          call spew (luout_, w, 'BCGS: AFTER FIRST PRECONDITIONING (w):')
       end if
       sigma = Ubik_Dot_Product(q, w)
       if (ABS(sigma) < smallest(control)) then
          if (outmode_errors(control)) then
             call spew (luerr_, ' ')
             call spew (luerr_, '*******************************')
             call spew (luerr_, 'BCGS: ERROR')
             call spew (luerr_, '  Breakdown (sigma = 0.0)')
             write(out_string,905)         '  current iteration number: ', iter(control)
             call spew (luerr_, out_string)
             write(out_string,940)         '  residual norm at breakdown: ', rnorm(control)
             call spew (luerr_, out_string)
             call spew (luerr_, '*******************************')
          end if
          call set_status_breakdown (control)
          exit
       end if
       alpha = rhoold/sigma

       h = u - alpha*w
       v = alpha*(u + h)

       if (outmode_full_) then
          call spew (luout_, h, 'BCGS: BEFORE CALCULATION OF RESIDUAL (h):')
          call spew (luout_, v, 'BCGS: BEFORE CALCULATION OF RESIDUAL (v):')
       end if

       ! Call user routine to compute matrix-vector product.
       call metaMatVec (MatVec, control, v_vec, work)
       if (matmul_error(control)) then
          exit
       end if

       if (precond) then

          ! Call user routine to perform preconditioning.
          call metaPrecond (ApplyPrecond, control, work, w_vec)
          if (precond_error(control)) then
             exit
          end if
       else
          work = w
       end if

       ! Compute BCGS residual.
       r = r - w
       if (outmode_full_) then
          call spew (luout_, r, 'BCGS: AFTER CALCULATION OF RESIDUAL (r):')
       end if

       ! Compute || r ||.
       call set_rnorm (control, Norm(norm_type_, r))

       select case (method(control))
       case (method_BCGS)

          ! Compute new BCGS iterate.
          x = x + v
          rnorm_iter = rnorm(control)
       case (method_TFQMR)

          ! Compute 1st TFQMR iterate and check convergence.
          term = theta2*eta/alpha
          d = term*d + u
          theta = SQRT(rnorm(control, MAX(0,iter(control)-2))*rnorm(control)) / tau
          theta2 = theta*theta
          sigma = 1.0 / SQRT(1.0 + theta2)
          eta = sigma*sigma*alpha
          x = x + eta*d
          rnorm_qmr = tau*SQRT(REAL(iter(control)+1))
          call CheckConvergence (x_vec, b, control, MatVec, rnorm_arg=rnorm_qmr)
          call set_rnorm (control, rnorm(control, iter(control)-1))
          if (.not.next_iter(control) .or. theta == 0.0) then
             exit
          end if

          ! Compute 2nd TFQMR iterate and check convergence.
          tau = tau*theta*sigma
          term = theta2*eta/alpha
          d = term*d + h
          theta = rnorm(control, iter(control)-1) / tau
          theta2 = theta*theta
          sigma = 1.0 / SQRT(1.0 + theta2)
          tau = tau*theta*sigma
          eta = sigma*sigma*alpha
          x = x + eta*d
          rnorm_qmr = tau*SQRT(REAL(iter(control)))
          rnorm_iter = rnorm_qmr
       case default
          call set_status_internal_error (control)
          exit
       end select

       ! Check convergence.
       call CheckConvergence (x_vec, b, control, MatVec, rnorm_arg=rnorm_iter)
       if (method(control) == method_TFQMR) then
          call set_rnorm (control, rnorm(control, iter(control)-1))
       end if
       if (.not.next_iter(control)) then
          exit
       end if

       ! Begin new iteration.
       rho = Ubik_Dot_Product(q, r)
       if (ABS(rho) < smallest(control)) then
          if (outmode_errors(control)) then
             call spew (luerr_, ' ')
             call spew (luerr_, '*******************************')
             call spew (luerr_, 'BCGS: ERROR')
             call spew (luerr_, '  Breakdown (rho = 0.0)')
             write(out_string,905)         '  current iteration number: ', iter(control)
             call spew (luerr_, out_string)
             write(out_string,940)         '  residual norm at breakdown: ', rnorm(control)
             call spew (luerr_, out_string)
             call spew (luerr_, '*******************************')
          end if
          call set_status_breakdown (control)
          exit
       end if
       beta = rho/rhoold
       rhoold = rho

       u = beta*h + r
       g = u + beta*(beta*g + h)

    end do

    ! common clean up and exit point
888 continue

    ! Deallocate local arrays in reverse order.
    call destroy (w_vec, overlap_only=.true.)
    call destroy (r_vec, overlap_only=.true.)
    call destroy (x_vec, overlap_only=.true.)
    call destroy (v_vec, overlap_only=.true.)
    call destroy (g_vec, overlap_only=.true.)
    if (allocated(work)) then
       deallocate(work)
    end if
    if (allocated(w)) then
       deallocate(w)
    end if
    if (allocated(v)) then
       deallocate(v)
    end if
    if (allocated(g)) then
       deallocate(g)
    end if
    if (allocated(u)) then
       deallocate(u)
    end if
    if (allocated(r)) then
       deallocate(r)
    end if
    if (allocated(q)) then
       deallocate(q)
    end if
    if (allocated(h)) then
       deallocate(h)
    end if
    if (allocated(d)) then
       deallocate(d)
    end if

    ! exit w/o cleanup
999 continue

    ! Subtract time at start of solution from each element of cpu array
    ! to obtain true CPU time for each iteration.
    call normalize_cpu (control)

#include "formats.F"
    return
  end subroutine BCGS

end module Ubik_BCGS_module
