! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: Meta.F,v 1.6 2002/12/03 04:58:11 turner Exp $
! 
!*********************************************************************
module Ubik_Meta_module
  use ubik_types_module
  use ubik_control_type_module
  use ubik_vector_type_module
  implicit none

  private
  public :: &
       metaMatVec, &
       metaPrecond

contains

  !*****************************************************************
  !
  ! <para>
  ! Wrapper for matrix-vector multiplication, y = Ax.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   x - x-vector
  ! 
  !  Output:
  !   y - y-vector
  !
  !*****************************************************************
  recursive subroutine metaMatVec(MatVec, control, x, y)

#include "matvec_proto.F"

    ! <arguments>
    type(control_type), intent(inout) :: control
    type(vector_type), intent(inout) :: x
    real(real_type), intent(inout), dimension(:), target :: y
    ! </arguments>

    ! Local:
    integer(int_type) :: status

    ! deallocate aux vectors if solving globally
    if (scope_global(control)) then
       call destroy (x, overlap_only=.true.)
    end if

    ! increment number of matvecs performed
    call increment_matvecs (control)

    ! call user routine to compute matrix-vector product
    call MatVec (x, y, status)

    ! check return status
    if (status /= 0) then
       call set_status_matmul_error (control)
    end if

    return
  end subroutine metaMatVec

  !*****************************************************************
  !
  ! <para>
  ! Wrapper for preconditioning (solve Ax = b).
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   b - rhs
  ! 
  !  Output:
  !   x - unknowns
  !
  !*****************************************************************
  subroutine metaPrecond(ApplyPrecond, control, b, x)

#include "precond_proto.F"

    ! <arguments>
    real(real_type), intent(in), dimension(:) :: b
    type(control_type), intent(inout) :: control
    type(vector_type), intent(inout) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: status

    ! deallocate aux vectors if solving globally
    if (scope_global(control)) then
       call destroy (x, overlap_only=.true.)
    end if

    ! increment number of preconditioning applications performed
    call increment_preconds (control)

    ! call user routine to apply preconditioner
    x = 0.0d0
    call ApplyPrecond (b, x, status)

    ! check return status
    if (status /= 0) then
       call set_status_precond_error (control)
    end if

    return
  end subroutine metaPrecond

end module Ubik_Meta_module
