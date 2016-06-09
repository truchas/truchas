#include "f90_assert.fpp"

module ded_head_type

  use toolpath_type
  use gb_laser_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private
  
  type, public :: ded_head
    type(toolpath), pointer :: tp => null()
    type(gb_laser) :: laser
    real(r8) :: absorp
  contains
    procedure :: init
    procedure :: laser_irrad
  end type ded_head

contains

  subroutine init(this, params)
    use parameter_list_type
    use toolpath_table, only: toolpath_ptr
    class(ded_head), intent(out) :: this
    type(parameter_list) :: params
    character(:), allocatable :: tp_name
    type(parameter_list), pointer :: plist
    call params%get('toolpath', tp_name)
    this%tp => toolpath_ptr(tp_name)
    INSIST(associated(this%tp))
    call params%get('laser-absorp', this%absorp)
    plist => params%sublist('laser')
    call this%laser%init(plist)
  end subroutine init

  function laser_irrad(this, t, r)
    class(ded_head), intent(in) :: this
    real(r8), intent(in) :: t, r(:)
    real(r8) :: laser_irrad, dr(3)
    if (this%tp%is_flag_set(0)) then
      call this%tp%get_position(t, dr)
      dr = r - dr
      laser_irrad = this%absorp * this%laser%irrad(dr(1), dr(2), dr(3))
    else
      laser_irrad = 0.0_r8
    end if
  end function laser_irrad

end module ded_head_type
