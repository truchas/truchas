! -*- Mode: f90 -*-

!**********************************************************************
!
! <para>
! Compute the solution to a nonsymmetric system of linear equations
! of the form Ax=b by Feund's transpose-free quasi-minimal residuals 
! method (TFQMR).
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
!**********************************************************************
module Ubik_TFQMR_module
  implicit none

  private
  public :: TFQMR

contains

  recursive subroutine TFQMR (x, b, control, MatVec, ApplyPrecond)
    use ubik_types_module
    use ubik_constants_module
    use ubik_control_type_module
    use Ubik_BCGS_module

#include "matvec_proto.F"
#include "precond_proto.F"

    optional :: ApplyPrecond
    real(real_type), intent(in), dimension(:) :: b
    type(control_type), intent(inout) :: control
    real(real_type), intent(inout), dimension(:), target :: x
    ! </arguments>

    call set_method (control, method_TFQMR)
    call BCGS (x, b, control, MatVec, ApplyPrecond)

    return
  end subroutine TFQMR

end module Ubik_TFQMR_module
