! -*- Mode: f90 -*-

!*******************************************************************
! $Id: BCGSTAB.F,v 1.30 2002/12/12 01:07:52 turner Exp $
! 
! <para>
! Compute the solution to a system of linear equations of the form
! Ax=b by Bi-CGSTAB.
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
!  spew
! 
! <FUNCTIONS REQUIRED>
! 
!  Dot_Product
!  Norm
!
!*******************************************************************
module Ubik_BCGSTAB_module
  implicit none

  private
  public :: BCGSTAB

contains

!*******************************************************************

  recursive subroutine BCGSTAB (x, b, control, MatVec, ApplyPrecond)
    use ubik_types_module
    use ubik_constants_module
    use ubik_control_type_module
    use ubik_vector_type_module
    use Ubik_Convergence_module, only: CheckConvergence
    use Ubik_Meta_module
    use Ubik_ParallelInfo_module, only: initialize_scope
    use Ubik_Vector_module, only: Ubik_Dot_Product, &
                             Norm, &
                             Ubik_Sum, &
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
    real(real_type) :: alpha, beta, delta, rho, rhoold
    real(real_type), allocatable, dimension(:) :: p, r, r0, s, t, v
    real(real_type), allocatable, dimension(:), target :: w, z
    type(vector_type) :: x_vec, w_vec, z_vec

    ! reset return status
    call set_status_converged (control)

    ! local alias for number of unknowns
    nunk = SIZE(x)

    ! check validity of stopping test - modify if necessary
    call check_stopping_nonstationary (control)

    ! initialize control params
    call initialize (control)

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
          call Spew (luerr_, 'BCGSTAB: WARNING')
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

    ! Nullify aux vector pointers in vector derived types.
    call create (z_vec, overlap_only=.true.)
    call create (w_vec, overlap_only=.true.)
    call create (x_vec, overlap_only=.true.)

    ! Allocate local arrays.
    allocate (p(nunk), r(nunk), r0(nunk), s(nunk), t(nunk),  &
         v(nunk), w(nunk), z(nunk), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if

    ! Initialize local arrays.
    p = 0.0
    r = 0.0
    r0 = 0.0
    s = 0.0
    t = 0.0
    v = 0.0
    w = 0.0
    z = 0.0

    ! Initialize vector derived types.
    call set_values_ptr (z_vec, z)
    call set_values_ptr (w_vec, w)
    call set_values_ptr (x_vec, x)

    ! Write initial system if desired.
    if (outmode_system(control)) then
       call spew (luout_, x, b, 'BCGSTAB on entry: Initial guess and source (x,b):')
    end if

    ! Compute || b || if needed.
    if (need_bnorm(control)) then
       call set_bnorm (control, Norm(norm_type_, b))
       if (bnorm(control) < smallest(control)) then
          call set_status_invalid_input (control)
          if (outmode_errors(control)) then
             call spew (luerr_, ' ')
             call spew (luerr_, '*******************************')
             call spew (luerr_, 'BCGSTAB: ERROR')
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

    ! ************* START MAIN LOOP *************

    do
       if (outmode_full_) then
          call spew (luout_, ' ')
          write(out_string,*) '*** ITERATION NUMBER: ', iter(control)
          call spew (luout_, out_string)
       end if

       ! Compute residual and || r ||.
       if (update_residual(control)) then

          ! Use r = b - Ax.  In this case also set rnormt = rnorm since it
          ! will be input to CheckConvergence rather than output.
          if (outmode_iterates(control)) then
             call spew (luout_, ' ')
             write(out_string,905) 'Computing true residual at iteration: ', iter(control)
             call spew (luout_, out_string)
             call spew (luout_, ' ')
          end if

          ! Call user routine to compute matrix-vector product.
          call metaMatVec (MatVec, control, x_vec, r)
          if (matmul_error(control)) then
             exit
          end if
          r = b - r
          call set_rnorm (control, Norm(norm_type_, r))
          call set_rnormt (control, rnorm(control))
       else

          ! Use recursive formula
          r = s - delta*t
          call set_rnorm (control, Norm(norm_type_, r))
       end if
       if (outmode_full_) then
          call spew (luout_, r, 'BCGSTAB: RESIDUAL (r):')
       end if

       ! If the residual norm is "small enough", simply declare victory.
       if (rnorm(control) <= smallest(control)) then
          call set_status_converged (control)
          exit
       end if

       ! Check for convergence.  Note that on entry iter = 0 and CheckConvergence increments.
       call CheckConvergence (x_vec, b, control, MatVec)
       if (.not.next_iter(control)) then
          exit
       end if

       if (iter(control) == 1) then
          r0 = r
       else
          rhoold = rho
       end if
       rho = Ubik_Dot_Product(r0, r)
       if (ABS(rho) < smallest(control)) then
          if (outmode_errors(control)) then
             call spew (luerr_, ' ')
             call spew (luerr_, '*******************************')
             call spew (luerr_, 'BCGSTAB: ERROR')
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

       if (iter(control) == 1) then
          p = r
       else
          beta = (rho/rhoold) * (alpha/delta)
          p = r + beta*(p - delta*v)
       end if

       ! Precondition.
       if (precond) then

          ! Call user routine to perform preconditioning.
          call metaPrecond (ApplyPrecond, control, p, w_vec)
          if (precond_error(control)) then
             exit
          end if
          if (outmode_full_) then
             call spew (luout_, w, 'BCGSTAB: AFTER 1st PRECONDITIONING (w):')
          end if
       else

          ! No preconditioning.
          w = p
       end if

       call metaMatVec (MatVec, control, w_vec, v)
       if (matmul_error(control)) then
          exit
       end if
       if (outmode_full_) then
          call spew (luout_, t, 'BCGSTAB: AFTER 1st MATVEC (t):')
       end if

       alpha = rho / Ubik_Dot_Product(r0, v)
       s = r - alpha*v

       ! should test for convergence based on ||s|| here

       ! Precondition.
       if (precond) then

          ! Call user routine to perform preconditioning.
          call metaPrecond (ApplyPrecond, control, s, z_vec)
          if (precond_error(control)) then
             exit
          end if
          if (outmode_full_) then
             call spew (luout_, z, 'BCGSTAB: AFTER 2nd PRECONDITIONING (z):')
          end if
       else

          ! No preconditioning.
          z = s
       end if

       call metaMatVec (MatVec, control, z_vec, t)
       if (matmul_error(control)) then
          exit
       end if
       if (outmode_full_) then
          call spew (luout_, t, 'BCGSTAB: AFTER 2nd MATVEC (t):')
       end if

       ! Compute new iterate.
       delta = Ubik_Dot_Product(s,t) / Ubik_Dot_Product(t,t)
       x = x + alpha*w + delta*z
    end do

    ! common clean up and exit point
888 continue

    ! Deallocate local arrays.
    call destroy (w_vec, overlap_only=.true.)
    call destroy (x_vec, overlap_only=.true.)
    call destroy (z_vec, overlap_only=.true.)
    if (allocated(p)) then
       deallocate(p)
    end if
    if (allocated(r)) then
       deallocate(r)
    end if
    if (allocated(r0)) then
       deallocate(r0)
    end if
    if (allocated(s)) then
       deallocate(s)
    end if
    if (allocated(t)) then
       deallocate(t)
    end if
    if (allocated(v)) then
       deallocate(v)
    end if
    if (allocated(w)) then
       deallocate(w)
    end if

    ! exit w/o cleanup
999 continue

    ! Subtract time at start of solution from each element of cpu array
    ! to obtain true CPU time for each iteration.
    call normalize_cpu (control)

#include "formats.F"
    return
  end subroutine BCGSTAB

end module Ubik_BCGSTAB_module
