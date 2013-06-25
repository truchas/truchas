!!
!! BDF2_CONTROLLER
!!

#include "f90_assert.fpp"

module bdf2_controller

  use bdf2_kinds
  implicit none
  private
  
  public :: bdf2_set_param, destroy, error_norm
  
  type, public :: bdf2_control
    real(kind=rk) :: hmin = 0.0_rk        ! step size lower bound
    real(kind=rk) :: hmax = huge(1.0_rk)  ! step size upper bound
    integer       :: mtry = 10            ! Maximum number of attempts at a step.
    integer       :: mitr = 5             ! Maximum number of BCE step iterations
    real(kind=rk) :: ntol = 0.1_rk        ! BCE step error tolerance (relative to 1)
    real(kind=rk) :: vtol = 1.0e-3_rk     ! Vector drop tolerance (FPI accelerator)
    integer       :: mvec = 0             ! Maximum number of vectors (FPI accelerator)
    real(kind=rk) :: rtol = 0.0_rk              ! relative error tolerance
    real(kind=rk), pointer :: atol(:) => null() ! absolute error tolerance
  end type bdf2_control
  
  interface destroy
    module procedure destroy_bdf2_control
  end interface

contains

  subroutine destroy_bdf2_control (control)
    type(bdf2_control), intent(inout) :: control
    type(bdf2_control) :: default
    if (associated(control%atol)) deallocate(control%atol)
    control = default ! default initialize
  end subroutine destroy_bdf2_control

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF2_SET_PARAM
 !!

  subroutine bdf2_set_param (control, hmin, hmax, mtry, atol, rtol, mitr, ntol, mvec, vtol)
  
    type(bdf2_control), intent(inout) :: control
    real(kind=rk), intent(in), optional :: hmin, hmax, atol(:), rtol, ntol, vtol
    integer, intent(in), optional :: mtry, mitr, mvec
  
    if (present(hmin)) control%hmin = hmin
    if (present(hmax)) control%hmax = hmax
    if (present(mtry)) control%mtry = mtry
    
    if (present(atol)) then
      if (.not.associated(control%atol)) allocate(control%atol(size(atol)))
      control%atol = atol
    end if
    if (present(rtol)) control%rtol = rtol
    
    if (present(mitr)) control%mitr = mitr
    if (present(ntol)) control%ntol = ntol
    if (present(mvec)) control%mvec = mvec
    if (present(vtol)) control%vtol = vtol
    
  end subroutine bdf2_set_param
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERROR_NORM
 !!
 
  function error_norm (control, u, du) result (norm)
  
    type(bdf2_control), intent(in) :: control
    real(kind=rk),      intent(in) :: u(:), du(:)
    real(kind=rk) :: norm
    
    norm = maxval(abs(du)/(control%atol + control%rtol*abs(u)))
    
  end function error_norm
  
end module bdf2_controller
