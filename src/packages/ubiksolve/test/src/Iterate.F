! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: Iterate.F,v 1.5 2002/07/25 17:28:58 turner Exp $
! 
!*********************************************************************
module Iterate_module
  use Ubik_module
  implicit none

  private
  public :: &
       iterate_control_type, &
       iterate_control, &
       SetIterateDefaults, &
       Iterate

  type iterate_control_type

     ! method - controls iteration method [DEFAULT = Ubik_method_SOR]
     integer(Ubik_int_type) :: method

     ! omega - relaxation parameter [DEFAULT = 1.0]
     real(Ubik_real_type) :: omega

  end type iterate_control_type

  type(iterate_control_type) :: iterate_control

contains

  !*********************************************************************
  ! Sets reasonable default values for control parameters.
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In/Out:
  !   iterate_control - structure containing iteration control parameters
  !*********************************************************************
  subroutine SetIterateDefaults (iterate_control)
    type(iterate_control_type), intent(out) :: iterate_control

    iterate_control%method = Ubik_method_SOR ! SOR
    iterate_control%omega = 1.0              ! relaxation parameter

    return
  end subroutine SetIterateDefaults

  subroutine Iterate_Full (coeff, b, x, status)
    real(Ubik_real_type), intent(in), dimension(:,:) :: coeff

#include "Iterate-guts.F"

  end subroutine Iterate_Full

  !*********************************************************************

  subroutine Iterate_ELL (coeff, b, x, status)
    type(Ubik_ELL_matrix), intent(in) :: coeff

#include "Iterate-guts.F"

  end subroutine Iterate_ELL

  !*********************************************************************
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   n - number of unknowns
  !   b - b-vector (source)
  ! 
  !  Input/Output:
  !   x - x-vector (unknown)
  ! 
  !*********************************************************************
  subroutine Iterate (b, x, status)
    use Full_data_module
    use ELL_data_module
    use test_data_module

    real(Ubik_real_type), intent(in), dimension(:) :: b
    real(Ubik_real_type), intent(inout), dimension(:) :: x
    integer(Ubik_int_type), intent(out) :: status

    select case (prob)
    case (Saad_Full)
       call Iterate_Full (coeff_full, b, x, status)
    case (Saad_ELL)
       call Iterate_ELL (coeff_ell, b, x, status)
    end select

    return
  end subroutine Iterate

end module Iterate_module
