! -*- Mode: f90 -*-

!**********************************************************************
! $Id: Stationary.F,v 1.16 2002/12/03 04:58:11 turner Exp $
! 
! <para>
! Computes the solution to a linear system of equations of the
! form Ax=b by a stationary iterative method (Jacobi, Jacobi with
! relaxation, Gauss-Seidel, SOR, Symmetric Gauss-Seidel, or SSOR).
! 
! An error message is written out and control is returned to the calling
! routine if the maximum number of iterations is exceeded.  At that time,
! err(control%iter) can be examined to determine if execution can continue.
! </para>
! 
! <PARAMETER LIST>
! 
!  Input:
!   b - source vector
! 
!  In/Out:
!   control - structure containing control parameters
!   x - solution vector (whatever is in x on entry is used as an
!       initial guess)
! 
! <SUBROUTINES REQUIRED>
! 
!  CheckConvergence
!  WriteSystem
!  WriteVector
! 
! <FUNCTIONS REQUIRED>
! 
!  Norm
!
!**********************************************************************
module Ubik_Stationary_module
  implicit none

  private
  public :: Stationary

contains

  subroutine Stationary (x, b, control, MatVec, Iterate)
    use ubik_types_module
    use ubik_constants_module
    use ubik_control_type_module
    use ubik_vector_type_module
    use Ubik_Convergence_module, only: CheckConvergence
    use Ubik_Meta_module, only: metaMatVec
    use Ubik_ParallelInfo_module, only: initialize_scope
    use Ubik_Vector_module, only: Norm, spew
    use Ubik_WriteString_module

#include "matvec_proto.F"
#include "iterate_proto.F"

    real(real_type), intent(in), dimension(:) :: b
    type(control_type), intent(inout) :: control
    real(real_type), intent(inout), dimension(:), target :: x
    ! </arguments>

    ! Local:
    logical(logical_type) :: outmode_full_, outmode_system_
    logical(logical_type) :: scope_global_
    integer(int_type) :: nunk, luout_, norm_type_, status
    real(real_type), allocatable, dimension(:) :: w
    type(vector_type) :: x_vec

    ! reset return status
    call set_status_converged (control)

    ! local alias for number of unknowns
    nunk = SIZE(x)

    ! check validity of stopping test - modify if necessary
    ! we know that after this stopping_true(control) will be true
    ! unless stopping_relchg(control) is true
    call check_stopping_stationary (control)

    ! initialize control params
    call initialize (control)
    call initialize_scope (scope(control))

    ! shorthand
    luout_ = luout(control)
    norm_type_ = norm_type(control)
    outmode_system_ = outmode_system(control)
    outmode_full_ = outmode_full(control)
    scope_global_ = scope_global(control)

    ! Allocate local arrays.
    allocate (w(nunk), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 9999
    end if

    ! Initialize vector derived types.
    call set_values_ptr (x_vec, x)

    ! Write initial system if desired.
    if (outmode_system_) then
       call spew (luout_, x, b, 'Stationary on entry: Initial guess and source (x,b):')
    end if

    ! Initialization dependent on the stopping test.
    if (need_bnorm(control)) then

       ! Compute || b || if needed.
       call set_bnorm (control, Norm(norm_type_, b))
    end if

    ! Initialize xold if needed.
    if (need_xold(control)) then

       call set_xold (control, x)

    else if (stopping_r0(control)) then

       ! Compute initial residual if using || r || / || r0 || stopping
       ! test and write out if desired.

       call metaMatVec (MatVec, control, x_vec, w)
       if (matmul_error(control)) then
          go to 9999
       end if
       w = b - w
       if (outmode_full_) then
          call spew (luout_, w, 'Stationary: INITIAL RESIDUAL (w):')
       end if
       call set_rnormt (control, Norm(norm_type_, w))
       if (rnormt(control) < TINY(0.0)) then
          call set_status_converged (control)
          go to 9999
       end if
    end if

    ! Main loop.
    do
       ! Deallocate aux vector pointers in vector derived types if solving globally.
       if (scope_global_) then
          call destroy (x_vec, .true.)
       end if

       call Iterate (b, x, status)

       call CheckConvergence (x_vec, b, control, MatVec)
       if (.not.next_iter(control)) then
          exit
       end if
    end do

9999 continue

    ! Deallocate local arrays.
    if (allocated(w)) deallocate (w)
    call destroy (x_vec, overlap_only=.true.)

    ! Subtract time at start of solution from each element of cpu array
    ! to obtain true CPU time for each iteration.
    call normalize_cpu (control)

    return 
  end subroutine Stationary

end module Ubik_Stationary_module
