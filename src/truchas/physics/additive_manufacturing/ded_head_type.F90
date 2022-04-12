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
    logical :: laser_is_on, fade_is_on
    real(r8):: tau, a0, t0
  contains
    procedure :: init
    procedure :: next_tp_segment
    procedure :: laser_irrad => irrad
    procedure, private :: fade_factor
  end type ded_head

contains

  subroutine init(this, params, t)
    use parameter_list_type
    use toolpath_table, only: toolpath_ptr
    use laser_irrad_factory
    class(ded_head), intent(out) :: this
    type(parameter_list) :: params
    real(r8), intent(in) :: t
    character(:), allocatable :: tp_name
    type(parameter_list), pointer :: plist
    call params%get('toolpath', tp_name)
    this%tp => toolpath_ptr(tp_name)
    INSIST(associated(this%tp))
    call this%tp%set_segment(t)
    plist => params%sublist('laser')
    call alloc_laser_irrad(this%laser, plist)
    call params%get('laser-absorp', this%absorp)
    call params%get('laser-time-constant', this%tau)
    this%laser_is_on = this%tp%is_flag_set(0)
    this%t0 = this%tp%prec_flag_flip_time(0)
    this%a0 = merge(0.0_r8, 1.0_r8, this%laser_is_on)
    this%fade_is_on = ((this%tau > 0.0_r8) .and. (this%t0 > -huge(1.0_r8)))
  end subroutine init

  subroutine next_tp_segment(this, t)
    class(ded_head), intent(inout) :: this
    real(r8), intent(in) :: t
    if (t /= this%tp%final_time()) return ! do not advance to the next segment
    call this%tp%next_segment()
    if (this%tp%is_flag_set(0) .neqv. this%laser_is_on) then ! laser flipped state
      this%fade_is_on = (this%tau > 0.0_r8)
      if (this%fade_is_on) then
        this%a0 = fade_factor(this, t)  ! new starting value
        this%t0 = t ! new starting time
      end if
      this%laser_is_on = .not.this%laser_is_on  ! new laser state
    end if
  end subroutine next_tp_segment

  pure function fade_factor(this, t) result(a)
    class(ded_head), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8) :: a
    if (this%laser_is_on) then  ! fade to 1
      if (this%fade_is_on) then
        a = 1 - (1 - this%a0)*exp((this%t0 - t)/this%tau)
      else
        a = 1
      end if
    else  ! fade to 0
      if (this%fade_is_on) then
        a = this%a0*exp((this%t0 - t)/this%tau)
      else
        a = 0
      end if
    end if
  end function fade_factor

  function irrad(this, t, r)
    class(ded_head), intent(in) :: this
    real(r8), intent(in) :: t, r(:)
    real(r8) :: irrad, dr(3)
    call this%tp%get_position(t, dr)
    dr = r - dr
    irrad = this%absorp * this%fade_factor(t) * this%laser%irrad(t, dr(1), dr(2), dr(3))
  end function irrad

end module ded_head_type
