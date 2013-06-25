! -*- Mode: f90 -*-

!**********************************************************************
! $Id: GMRES.F,v 1.32 2002/12/05 00:39:06 turner Exp $
!**********************************************************************
module Ubik_GMRES_module
  use ubik_types_module
  implicit none

  interface GMRES
     module procedure GMRES_Left_Radicati
  end interface

  private
  public :: &
       GMRES, &
       GMRES_Left_Radicati, &
       GMRES_Left_Saad, &
!       GMRES_Right_Saad, &
       FGMRES

contains

  ! <para>
  ! Compute Givens rotation.
  ! see Golub and Van Loan, 2nd Ed., p. 202 (Algorithm 5.1.5)
  ! </para>
  subroutine Givens (a, b, c, s)

    ! <arguments>
    real(real_type), intent(in) :: a, b
    real(real_type), intent(out) :: c, s
    ! </arguments>
    real(real_type) :: tau

    if (ABS(b) < TINY(0.0)) then
       c = 1.0
       s = 0.0
    else
       if (ABS(b) > ABS(a)) then
          tau = -a/b
          s = 1.0 / (SQRT(1.0 + tau*tau))
          c = s*tau
       else
          tau = -b/a
          c = 1.0 / (SQRT(1.0 + tau*tau))
          s = c*tau
       end if
    end if

    return
  end subroutine Givens

  !**********************************************************************
  ! 
  ! <para>
  ! Compute the solution to a nonsymmetric system of linear equations of the
  ! form Ax=b by the restarted Generalized Minimal Residual method, GMRES(k),
  ! where k is the frequency of restart, with left preconditioning.
  ! 
  ! An error message is written out and control is returned to the calling
  ! routine if a breakdown in iterates is encountered or if the maximum number
  ! of iterations is exceeded.  At that time, err(control%iter) can be 
  ! examined to determine if execution can continue.
  !
  ! This implementation is a slight modificiation of the one given in paper
  ! by G. Radicati di Brozolo and Y. Robert.
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
  !  SolveUTriang
  !  spew
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  Dot_Product
  !  Norm
  !
  !**********************************************************************
  recursive subroutine GMRES_Left_Radicati (x, b, control, MatVec, ApplyPrecond)

    use ubik_constants_module
    use ubik_control_type_module
    use ubik_vector_type_module
    use Ubik_Convergence_module, only: CheckConvergence
    use Ubik_FlushUnit_module
    use Ubik_Full_module, only: SolveUTriang
    use Ubik_Meta_module
    use Ubik_ParallelInfo_module, only: initialize_scope
    use Ubik_Vector_module, only: Ubik_Dot_Product, &
                             Norm, &
                             Swap, &
                             spew
    use Ubik_WriteString_module

#include "GMRES-decl.F"

    ! Local:
    real(real_type) :: gamma
    real(real_type) :: hijp, hjjinv
    real(real_type), allocatable, dimension(:,:), target :: v
    type(vector_type), allocatable, dimension(:) :: v_vec

    ! string containing name of this routine
    routine = 'GMRES_Left_Radicati'

#include "GMRES-init.F"

    ! allocate and initialize vector type corresponding to v
    allocate (v_vec(m+1), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if
    do j=1,m+1
       call set_values_ptr (v_vec(j), v(:,j))
    end do

    OUTER: do
       ! Compute residual (w=b-Ax) and write out if desired.

       ! Call user routine to compute matrix-vector product.
       call metaMatVec (MatVec, control, x_vec, w)
       if (matmul_error(control)) then
          go to 888
       end if

       ! compute unpreconditioned residual
       w = b - w

       if (outmode_full_) then
          call spew (luout_, w, TRIM(routine)//': AT START OF OUTER LOOP (w):')
       end if

       ! Initialize directions and write out if desired (formerly passed v(:,1)).
       if (precond) then

          ! Call user routine to perform preconditioning.
          call metaPrecond (ApplyPrecond, control, w, v_vec(1))
          if (precond_error(control)) then
             exit
          end if
       else
          v(:,1) = w
       end if

       ! v(:,1) now contains the preconditioned residual
       if (outmode_full_) then
          call spew (luout_, v(:,1), TRIM(routine)//': JUST BEFORE INNER LOOP (v):')
       end if

       ! clear vector that will hold residual norm estimates
       r = 0.0

       ! set initial residual
       r(1) = Norm(2, v(:,1))
       call set_rnorm (control, r(1))
       h(1,1) = r(1)

       ! if this is first iteration, check for convergence
       if (iter(control) == 0) then

          ! initialize true residual norm if necessary
          if (stopping_true(control)) then
             call set_rnormt (control, rnorm(control))
          end if

          if (rnorm(control) <= smallest(control)) then

             ! If the residual norm is smaller than control%tiny, simply declare victory.
             call set_status_converged (control)
          else

             ! Check to see if initial guess meets convergence criterion.
             ! Note that on entry iter = 0 and CheckConvergence increments iter.  
             call CheckConvergence (x_vec, b, control, MatVec)
          end if
          if (.not.next_iter(control)) then
             exit
          end if
       end if

       INNER: do j=1,m

          if (outmode_full_) then
             call spew (luout_, ' ')
             write(out_string,905) '*** ITERATION NUMBER: ', iter(control)
             call spew (luout_, out_string)
          end if

          if (h(j,j) < smallest(control)) then

             ! Not quite sure what should be done here.  This really shouldn't happen,
             ! since h(j,j) = 0 should indicate that the solution has been found, so
             ! convergence should have been detected.  So call it a breakdown, since
             ! it probably indicates a problem.
             hjjinv = 1.0
          else
             hjjinv = 1.0 / h(j,j)
          end if

          ! Step 1 (Kernel 2).

          ! Call user routine to compute matrix-vector product.
          call metaMatVec (MatVec, control, v_vec(j), w)
          if (matmul_error(control)) then
             exit OUTER
          end if

          w = hjjinv*w
          if (outmode_full_) then
             call spew (luout_, w, TRIM(routine)//': AFTER START OF INNER LOOP (w):')
          end if

          ! Apply preconditioner (formerly passed v(:,j+1)).
          if (precond) then

             ! Call user routine to perform preconditioning.
             call metaPrecond (ApplyPrecond, control, w, v_vec(j+1))
             if (precond_error(control)) then
                exit OUTER
             end if
          else

             ! No preconditioning.
             v(:,j+1) = w
          end if

          if (outmode_full_) then
             call spew (luout_, v(:,j+1), TRIM(routine)//': AFTER PRECONDITIONING (v):')
          end if

          h(1,j+1) = Ubik_Dot_Product(v(:,1), v(:,j+1))

          ! Step 2 (Kernel 3).
          do i=1,j-1
             v(:,j+1) = v(:,j+1) - h(i,j+1)*v(:,i)
             h(i+1,j+1) = Ubik_Dot_Product(v(:,i+1), v(:,j+1))
          end do
          if (outmode_full_) then
             call spew (luout_, v(:,j+1), TRIM(routine)//': AFTER STEP 2 (v):')
          end if

          ! Step 3.
          h(j,j+1) = h(j,j+1)*hjjinv

          ! Steps 4 (Kernel 4) and 5.
          v(:,j) = hjjinv*v(:,j)
          v(:,j+1) = v(:,j+1) - h(j,j+1)*v(:,j)
          h(j+1,j+1) = Norm(2, v(:,j+1))
          if (outmode_full_) then
             call spew (luout_, v(:,j+1), TRIM(routine)//': AFTER STEP 4 (v):')
          end if

          ! Step 6.
          do i=1,j-1

             ! 2x2 matrix multiplication.
             hijp = h(i,j+1)
             h(i,j+1) = c(i)*hijp + s(i)*h(i+1,j+1)
             h(i+1,j+1) = -s(i)*hijp + c(i)*h(i+1,j+1)
          end do

          ! Compute new Gram-Schmidt parameters via Givens rotation.
          gamma = SQRT( h(j,j+1)*h(j,j+1) + h(j+1,j+1)*h(j+1,j+1) )
          if (gamma < smallest(control)) then
             call set_status_breakdown (control)
             exit OUTER
          end if
          c(j) = h(j,j+1) / gamma
          s(j) = h(j+1,j+1) / gamma
!!$          call givens (h(j,j+1), h(j+1,j+1), c(j), s(j))
!!$          s(j) = -s(j)

          ! Compute residual norm estimate.
          r(j+1) = -s(j)*r(j)
          call set_rnorm (control, ABS(r(j+1)))
          r(j) = c(j)*r(j)

          h(j,j+1) = c(j)*h(j,j+1) + s(j)*h(j+1,j+1)

          ! Solve upper triangular system and compute new iterate *only*
          ! if needed for stopping test or to calculate the true residual
          ! and error estimate.  Note that this is costly.
          if (update_iterate) then
             rt(1:j) = r(1:j)
             call SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., &
                  a=h(1:j,2:j+1), x=rt(1:j))
             xt = x
             do i=1,j
                xt = xt + rt(i)*v(:,i)
             end do
          end if

          ! Check for convergence.
          call CheckConvergence (xt_vec, b, control, MatVec)
          if (.not.next_iter(control)) then
             exit INNER
          end if

       end do INNER

       ! Solve upper triangular system.

       ! Update iterate.
       if (update_iterate) then

          ! May have already done it.
          x = xt
       else

          ! If not, first solve upper triangular system.
          n = MIN(j,m)
          call SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., &
               a=h(1:n,2:n+1), x=r(1:n))

          ! Now update iterate.
          do i=1,n
             x = x + r(i)*v(:,i)
          end do
          if (outmode_iterates_) then
             call spew (luout_, x, TRIM(routine)//': ITERATE AT END OF OUTER LOOP (x):')
          end if

          ! Copy iterate into xt if needed for stopping test.
          if (stopping_uses_x) then
             xt = x
          end if
       end if

       ! End of outer loop.
       if (.NOT.next_iter(control)) then
          exit
       end if
    end do OUTER

    ! common clean up and exit point
888 continue

    ! deallocate local arrays
    do j=1,m+1
       call destroy (v_vec(j), overlap_only=.true.)
    end do
    if (allocated(v_vec)) then
       deallocate(v_vec)
    end if

#include "GMRES-cleanup.F"

    return
  end subroutine GMRES_Left_Radicati

  !**********************************************************************
  ! 
  ! <para>
  ! Compute the solution to a nonsymmetric system of linear equations of the
  ! form Ax=b by the restarted Generalized Minimal Residual method, GMRES(k),
  ! where k is the frequency of restart, with left preconditioning.
  ! 
  ! An error message is written out and control is returned to the calling
  ! routine if a breakdown in iterates is encountered or if the maximum number
  ! of iterations is exceeded.  At that time, err(control%iter) can be 
  ! examined to determine if execution can continue.
  !
  ! This implementation is based on the one given in Saad's original paper.
  ! </para>
  !
  !**********************************************************************
  recursive subroutine GMRES_Left_Saad (x, b, control, MatVec, ApplyPrecond)

    use ubik_constants_module
    use ubik_control_type_module
    use ubik_vector_type_module
    use Ubik_Convergence_module, only: CheckConvergence
    use Ubik_FlushUnit_module
    use Ubik_Full_module, only: SolveUTriang
    use Ubik_Meta_module
    use Ubik_ParallelInfo_module, only: initialize_scope
    use Ubik_Vector_module, only: Ubik_Dot_Product, &
                             Norm, &
                             Swap, &
                             spew
    use Ubik_WriteString_module

#include "GMRES-decl.F"

    ! Local:
    real(real_type) :: beta, delta, himj
    real(real_type), allocatable, dimension(:,:), target :: v
    type(vector_type), allocatable, dimension(:) :: v_vec

    ! string containing name of this routine
    routine = 'GMRES_Left_Saad'

#include "GMRES-init.F"

    ! allocate and initialize vector type corresponding to v
    allocate (v_vec(m+1), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if
    do j=1,m+1
       call set_values_ptr (v_vec(j), v(:,j))
    end do

    ! Compute residual (w=b-Ax) and write out if desired.                                    Step 1.
    OUTER: do

       ! Call user routine to compute matrix-vector product.
       call metaMatVec (MatVec, control, x_vec, w)
       if (matmul_error(control)) then
          exit OUTER
       end if

       ! compute unpreconditioned residual
       w = b - w

       if (outmode_full_) then
          call spew (luout_, w, TRIM(routine)//': AT START OF OUTER LOOP (w = b-Ax):')
       end if

       ! Initialize directions and write out if desired (formerly passed v(:,1)).
       if (precond) then

          ! Call user routine to perform preconditioning.
          call metaPrecond (ApplyPrecond, control, w, v_vec(1))
          if (precond_error(control)) then
             exit
          end if
       else
          v(:,1) = w
       end if

       ! v(:,1) now contains the preconditioned residual
       call set_rnorm (control, Norm(norm_type_, v(:,1)))

       ! If this is first iteration, compute || r0 ||.
       if (iter(control) == 0) then

          ! hereafter rnormt will be 
          call set_rnormt (control, rnorm(control))

          if (rnorm(control) <= smallest(control)) then

             ! If the residual norm is smaller than control%tiny, simply declare victory.
             call set_status_converged (control)
          else

             ! Check to see if initial guess meets convergence criterion.
             ! Note that on entry iter = 0 and CheckConvergence increments iter.  
             call CheckConvergence (x_vec, b, control, MatVec, rnorm_arg=rnormt(control))
          end if
          if (.not.next_iter(control)) then
             exit
          end if
       end if

       ! clear vector that will hold residual norm estimates
       r = 0.0

       ! set initial residual - recall that CheckConvergence incremented iteration
       ! count
       r(1) = rnorm(control, iter(control)-1)
       v(:,1) = v(:,1) / r(1)

       if (outmode_full_) then
          call spew (luout_, v(:,1), TRIM(routine)//': JUST BEFORE INNER LOOP (v):')
       end if

       INNER: do j=1,m                                                                     ! Step 2.

          if (outmode_full_) then
             call spew (luout_, ' ')
             write(out_string,905) '*** ITERATION NUMBER: ', iter(control)
             call spew (luout_, out_string)
          end if

          ! Step 3.

          ! Call user routine to compute matrix-vector product.
          call metaMatVec (MatVec, control, v_vec(j), w)
          if (matmul_error(control)) then
             exit OUTER
          end if

          if (outmode_full_) then
             call spew (luout_, w, TRIM(routine)//': AFTER START OF INNER LOOP (w = A*v_j):')
          end if

          ! Apply preconditioner.
          if (precond) then

             ! Call user routine to perform preconditioning.
             ! CHECK: should the following be j or j+1???
             call metaPrecond (ApplyPrecond, control, w, v_vec(j))
             if (precond_error(control)) then
                exit OUTER
             end if
          else

             ! No preconditioning.
             v(:,j) = w
          end if

          if (outmode_full_) then
             call spew (luout_, v(:,j), TRIM(routine)//': AFTER PRECONDITIONING (v):')
          end if

          do i=1,j                                                                         ! Step 4.
             h(i,j) = Ubik_Dot_Product(w, v(:,i))                                          ! Step 5.
             w = w - h(i,j)*v(:,i)                                                         ! Step 6.
          end do                                                                           ! Step 7.
          if (outmode_full_) then
             call spew (luout_, w, TRIM(routine)//': AFTER STEP 7 (w):')
          end if

          ! Step 8.
          h(j+1,j) = Norm(2, w)
          v(:,j+1) = w / h(j+1,j)
          if (outmode_full_) then
             call spew (luout_, v(:,j+1), TRIM(routine)//': AFTER STEP 8 (v):')
          end if

          ! Step 10.
          do i=2,j

             ! 2x2 matrix multiplication.
             himj = h(i-1,j)
             h(i-1,j) = c(i-1)*himj + s(i-1)*h(i,j)
             h(i,j) = -s(i-1)*himj + c(i-1)*h(i,j)
          end do

          ! Compute new Gram-Schmidt parameters.
          delta = SQRT( h(j,j)*h(j,j) + h(j+1,j)*h(j+1,j) )
          c(j) = h(j,j) / delta
          s(j) = h(j+1,j) / delta

          ! Compute residual norm estimate.
          h(j,j) = c(j)*h(j,j) + s(j)*h(j+1,j)
          call set_rnorm (control, ABS(h(j,j)))

          r(j+1) = -s(j)*r(j)
          r(j) = c(j)*r(j)

          ! Solve upper triangular system and compute new iterate *only*
          ! if needed for stopping test or to calculate the true residual
          ! and error estimate.  Note that this is costly.
          if (update_iterate) then
             rt(1:j) = r(1:j)
             call SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., &
                  a=h(1:j,2:j+1), x=rt(1:j))
             xt = x
             do i=1,j
                xt = xt + rt(i)*v(:,i)
             end do
          end if

          ! Check for convergence.
          call CheckConvergence (xt_vec, b, control, MatVec)
          if (.not.next_iter(control)) then
             exit INNER
          end if

       end do INNER

       ! Solve upper triangular system.

       ! Update iterate.
       if (update_iterate) then

          ! May have already done it.
          x = xt
       else

          ! If not, first solve upper triangular system.
          n = MIN(j,m)
          call SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., &
               a=h(1:n,2:n+1), x=r(1:n))

          ! Now update iterate.
          do i=1,n
             x = x + r(i)*v(:,i)
          end do
          if (outmode_iterates_) then
             call spew (luout_, x, TRIM(routine)//': ITERATE AT END OF OUTER LOOP (x):')
          end if

          ! Copy iterate into xt if needed for stopping test.
          if (stopping_uses_x) then
             xt = x
          end if
       end if

       ! End of outer loop.
       if (.NOT.next_iter(control)) then
          exit
       end if
    end do OUTER

    ! common clean up and exit point
888 continue

    ! deallocate local arrays
    do j=1,m+1
       call destroy (v_vec(j), overlap_only=.true.)
    end do
    if (allocated(v_vec)) then
       deallocate(v_vec)
    end if

#include "GMRES-cleanup.F"

    return
  end subroutine GMRES_Left_Saad

  !**********************************************************************
  ! 
  ! <para>
  ! Compute the solution to a nonsymmetric system of linear equations of the
  ! form Ax=b by the restarted Flexible Generalized Minimal Residual method, 
  ! FGMRES(k), where k is the frequency of restart.  FGMRES(k) allows a
  ! different preconditioner to be applied at each iteration, but if a constant
  ! one is used it is equivalent to GMRES(k) with right preconditioning.
  ! 
  ! An error message is written out and control is returned to the calling
  ! routine if a breakdown in iterates is encountered or if the maximum number
  ! of iterations is exceeded.  At that time, err(control%iter) can be 
  ! examined to determine if execution can continue.
  !
  ! This implementation is a corrected version of the one given in Chapter 8
  ! of the book Numerical Linear Algebra on High-Performance Computers
  ! (Fig. 8.7, p. 177).
  ! </para>
  !
  !**********************************************************************
  recursive subroutine FGMRES (x, b, control, MatVec, ApplyPrecond)

    use ubik_constants_module
    use ubik_control_type_module
    use ubik_vector_type_module
    use Ubik_Convergence_module, only: CheckConvergence
    use Ubik_FlushUnit_module
    use Ubik_Full_module, only: SolveUTriang
    use Ubik_Meta_module
    use Ubik_ParallelInfo_module, only: initialize_scope
    use Ubik_Vector_module, only: Ubik_Dot_Product, &
                             Norm, &
                             Swap, &
                             spew
    use Ubik_WriteString_module

#include "GMRES-decl.F"

    ! Local:
    real(real_type) :: delta, hkmi, rho
    real(real_type), allocatable, dimension(:,:) :: v
    real(real_type), allocatable, dimension(:,:), target :: z
    type(vector_type), allocatable, dimension(:) :: z_vec

    ! string containing name of this routine
    routine = 'FGMRES'

#include "GMRES-init.F"

    ! allocate and initialize z and corresponding vector type
    allocate (z(nunk,m+1), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if
    z = 0.0
    allocate (z_vec(m+1), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if
    do i=1,m+1
       call set_values_ptr (z_vec(i), z(:,i))
    end do

    OUTER: do

       ! compute residual

       ! call user routine to compute matrix-vector product
       call metaMatVec (MatVec, control, x_vec, w)
       if (matmul_error(control)) then
          go to 888
       end if
       w = b - w

       if (outmode_full_) then
          call spew (luout_, w, TRIM(routine)//': residual (r = b-Ax):')
       end if

       ! clear vector that will hold residual norm estimates
       r = 0.0

       ! set initial residual
       r(1) = Norm(2, w)
       call set_rnorm (control, r(1))
       v(:,1) = w / r(1)

       INNER: do i=1,m

          if (outmode_full_) then
             call spew (luout_, ' ')
             write(out_string,905) '*** ITERATION NUMBER: ', iter(control)
             call spew (luout_, out_string)
          end if

          ! apply preconditioner
          if (precond) then

             ! call user routine to perform preconditioning
             call metaPrecond (ApplyPrecond, control, v(:,i), z_vec(i))
             if (precond_error(control)) then
                exit OUTER
             end if
          else

             ! no preconditioning
             z(:,i) = v(:,i)
          end if

          if (outmode_full_) then
             call spew (luout_, z(:,i), TRIM(routine)//': AFTER PRECONDITIONING (z_i):')
          end if

          ! call user routine to compute matrix-vector product
          call metaMatVec (MatVec, control, z_vec(i), w)
          if (matmul_error(control)) then
             exit OUTER
          end if

          if (outmode_full_) then
             call spew (luout_, w, TRIM(routine)//': AFTER MATVEC (w = A*z_i):')
          end if

          do k=1,i
             h(k,i) = Ubik_Dot_Product(v(:,k), w)
             w = w - h(k,i)*v(:,k)
          end do
          if (outmode_full_) then
             call spew (luout_, w, TRIM(routine)//': (w):')
          end if

          h(i+1,i) = Norm(2, w)
          if (h(i+1,i) < smallest(control)) then

             ! Not quite sure what should be done here.  This really shouldn't happen,
             ! but h(i+1,i) = 0 should indicate that the solution has been found, so
             ! just prevent the divide by zero that would otherwise occur here and
             ! convergence should be detected below.
             v(:,i+1) = w
          else
             v(:,i+1) = w / h(i+1,i)
          end if
          if (outmode_full_) then
             call spew (luout_, v(:,i+1), TRIM(routine)//': (v_i+1):')
          end if

          do k=2,i

             ! 2x2 matrix multiplication
             hkmi = h(k-1,i)
             h(k-1,i) = c(k-1)*hkmi + s(k-1)*h(k,i)
             h(k,i) = -s(k-1)*hkmi + c(k-1)*h(k,i)
          end do

          ! compute new Gram-Schmidt parameters
          delta = SQRT( h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i) )
          if (delta < smallest(control)) then
             call set_status_breakdown (control)
             exit OUTER
          end if
          c(i) = h(i,i) / delta
          s(i) = h(i+1,i) / delta
          h(i,i) = c(i)*h(i,i) + s(i)*h(i+1,i)

          ! compute residual norm estimate
          r(i+1) = -s(i)*r(i)
          r(i) = c(i)*r(i)

          rho = ABS(r(i+1))
          call set_rnorm (control, rho)

          ! Solve upper triangular system and compute new iterate *only*
          ! if needed for stopping test or to calculate the true residual
          ! and error estimate.  Note that this is costly.
          if (update_iterate) then
             rt(1:i) = r(1:i)
             call SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., &
                  a=h(1:i,1:i), x=rt(1:i))
             xt = x
             do k=1,i
                xt = xt + rt(k)*z(:,k)
             end do
          end if

          ! Check for convergence.
          call CheckConvergence (xt_vec, b, control, MatVec)
          if (.not.next_iter(control)) then
             n = i
             goto 777
          end if

       end do INNER
       n = m

777    continue

       ! Update iterate.
       if (update_iterate) then

          ! May have already done it.
          x = xt
       else

          ! If not, first solve upper triangular system.
          call SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., &
               a=h(1:n,1:n), x=r(1:n))

          ! Now update iterate.
          do i=1,n
             x = x + r(i)*z(:,i)
          end do
          if (outmode_iterates_) then
             call spew (luout_, x, TRIM(routine)//': ITERATE AT END OF OUTER LOOP (x):')
          end if

          ! Copy iterate into xt if needed for stopping test.
          if (stopping_uses_x) then
             xt = x
          end if
       end if

       ! End of outer loop.
       if (.NOT.next_iter(control)) then
          exit
       end if
    end do OUTER

    ! common clean up and exit point
888 continue

    ! deallocate local arrays
    do i=1,m+1
       call destroy (z_vec(i), overlap_only=.true.)
    end do
    if (allocated(z_vec)) then
       deallocate(z_vec)
    end if
    if (allocated(z)) then
       deallocate(z)
    end if

#include "GMRES-cleanup.F"

    return
  end subroutine FGMRES

end module Ubik_GMRES_module
