#include "f90_assert.fpp"

module ded_head_type

  use toolpath_type
  use laser_irrad_class
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: ded_head
    type(toolpath), pointer :: tp => null()
    class(laser_irrad), allocatable :: laser
    real(r8) :: absorp
  contains
    procedure :: init
    procedure :: laser_irrad => irrad
  end type ded_head

contains

  subroutine init(this, params)
    use parameter_list_type
    use toolpath_table, only: toolpath_ptr
    use laser_irrad_factory
    class(ded_head), intent(out) :: this
    type(parameter_list) :: params
    character(:), allocatable :: tp_name
    type(parameter_list), pointer :: plist
    call params%get('toolpath', tp_name)
    this%tp => toolpath_ptr(tp_name)
    INSIST(associated(this%tp))
    call params%get('laser-absorp', this%absorp)
    plist => params%sublist('laser')
    call alloc_laser_irrad(this%laser, plist)
  end subroutine init

  function irrad(this, t, r)
    class(ded_head), intent(in) :: this
    real(r8), intent(in) :: t, r(:)
    real(r8) :: irrad, dr(3)
    if (this%tp%is_flag_set(0)) then
      call this%tp%get_position(t, dr)
      dr = r - dr
      irrad = this%absorp * this%laser%irrad(dr(1), dr(2), dr(3))
    else
      irrad = 0.0_r8
    end if
  end function irrad

end module ded_head_type
