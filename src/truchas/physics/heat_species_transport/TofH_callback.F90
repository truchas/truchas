module TofH_callback

  use kinds, only: r8
  use property_mesh_function
  implicit none
  private
  
  public :: set_context_HofT, set_context_cell_H, f
  
  type :: context
    type(prop_mf), pointer :: HofT => null()
    integer  :: cell
    real(r8) :: H
  end type context
  type(context), save :: this
  
contains

  subroutine set_context_HofT (HofT)
    type(prop_mf), intent(in), target :: HofT
    this%HofT => HofT
  end subroutine set_context_HofT
  
  subroutine set_context_cell_H (cell, H)
    integer, intent(in) :: cell
    real(r8), intent(in) :: H
    this%cell = cell
    this%H = H
  end subroutine set_context_cell_H

  function f(x) result (fx)
    real(r8), intent(in) :: x
    real(r8) :: fx
    call pmf_eval (this%HofT, this%cell, (/x/), fx)
    fx = this%H - fx
  end function

end module TofH_callback
