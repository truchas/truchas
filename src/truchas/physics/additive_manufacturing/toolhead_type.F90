!!
!! TOOLHEAD_TYPE
!!
!! This module defines the TOOLHEAD derived type that describes a toolhead and
!! its attached energy and mass sources that are moving with respect to the
!! computational domain. An object of this type couples a TOOLPATH object that
!! describes the motion of the toolhead with reference configuration source
!! functions to yield moving sources. It is in essence a factory for allocating
!! standard abstract function objects that implement the moving sources for use
!! by physics models. Objects also handle advancing the toolpath between path
!! segments and updating the source functions as necessary.
!!
!! Currently only laser irradiation energy sources are implemented.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module toolhead_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use toolpath_type
  use laser_irrad_class
  use vector_func_class
  implicit none
  private

  type, public :: toolhead
    type(toolpath), allocatable :: tp
    class(laser_irrad), allocatable :: laser
    real(r8) :: absorp, dir(3)
    logical :: laser_is_on, fade_is_on
    real(r8):: tau, a0, t0
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: next_tp_segment
    procedure :: alloc_laser_func
    procedure :: laser_irrad => irrad
    procedure, private :: fade_factor
  end type toolhead

  type, extends(vector_func) :: laser_vector_func
    type(toolhead), pointer :: th => null() ! reference only
  contains
    procedure :: eval => laser_func_eval
    procedure :: eval_comp => laser_func_eval_comp
  end type laser_vector_func

contains

  !TODO: Add error checking of parameter list; currently assumed valid
  !      because the namelist reader checks and creates the parameter list.

  subroutine init(this, params)

    use parameter_list_type
    use toolpath_driver, only: alloc_toolpath
    use laser_irrad_factory

    class(toolhead), intent(out) :: this
    type(parameter_list) :: params

    character(:), allocatable :: tp_name
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: array(:)
    integer :: stat
    character(:), allocatable :: errmsg

    call params%get('toolpath', tp_name)
    call alloc_toolpath(this%tp, tp_name, stat, errmsg)
    INSIST(stat == 0)

    plist => params%sublist('laser')
    call alloc_laser_irrad(this%laser, plist)
    call params%get('laser-absorp', this%absorp)
    call params%get('laser-time-constant', this%tau)
    call params%get('laser-direction', array)
    this%dir = array / norm2(array) ! unit normal direction

  end subroutine init

  subroutine set_initial_state(this, t)
    class(toolhead), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%tp%set_segment(t)
    this%laser_is_on = this%tp%is_flag_set(0)
    this%t0 = this%tp%prec_flag_flip_time(0)
    this%a0 = merge(0.0_r8, 1.0_r8, this%laser_is_on)
    this%fade_is_on = ((this%tau > 0.0_r8) .and. (this%t0 > -huge(1.0_r8)))
  end subroutine

  !! Advance the toolhead to the next segment of its toolpath.
  subroutine next_tp_segment(this, t)
    class(toolhead), intent(inout) :: this
    real(r8), intent(in) :: t
    if (t /= this%tp%final_time()) return
    call this%tp%next_segment
    if (this%tp%is_flag_set(0) .neqv. this%laser_is_on) then ! laser flipped state
      this%fade_is_on = (this%tau > 0.0_r8)
      if (this%fade_is_on) then
        this%a0 = fade_factor(this, t)  ! new starting value
        this%t0 = t ! new starting time
      end if
      this%laser_is_on = .not.this%laser_is_on  ! new laser state
    end if
  end subroutine

  pure function fade_factor(this, t) result(a)
    class(toolhead), intent(in) :: this
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
    class(toolhead), intent(in) :: this
    real(r8), intent(in) :: t, r(:)
    real(r8) :: irrad, dr(3)
    call this%tp%get_position(t, dr)
    dr = r - dr
    !NB: This exploits a feature of the current implementations of laser%irrad where
    !    the profile is radially symmetric in the x and y coordinates, and needs to
    !    be refactored. FIXME
    dr = dr - dot_product(dr,this%dir)*this%dir
    dr(1) = norm2(dr)
    dr(2) = 0
    dr(3) = 0
    irrad = this%absorp * this%fade_factor(t) * this%laser%irrad(t, dr(1), dr(2), dr(3))
  end function irrad

  subroutine alloc_laser_func(this, f)
    class(toolhead), intent(in), target :: this
    class(vector_func), allocatable, intent(out) :: f
#ifdef GNU_PR103394
    block
      type(laser_vector_func), allocatable :: tmp
      allocate(tmp)
      tmp%dim = 3
      tmp%th => this
      call move_alloc(tmp, f)
    end block
#else
    allocate(f, source=laser_vector_func(dim=3,th=this))
#endif
  end subroutine

  !! Type bound EVAL procedure; delegates to the modules TOOLHEAD object.
  function laser_func_eval(this, x) result(fx)
    class(laser_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%dim)
    integer :: i
    fx = this%th%laser_irrad(x(1),x(2:4)) * this%th%dir
  end function

  function laser_func_eval_comp(this, i, x) result(fx)
    class(laser_vector_func), intent(in) :: this
    integer, intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    real(r8) :: f(this%dim)
    f = this%eval(x)
    fx = f(i)
  end function

end module toolhead_type
