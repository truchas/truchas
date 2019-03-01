!!
!! GB_LASER_TYPE
!!
!! Laser irradiance function based on the gaussian beam equation.
!!
!! Neil Carlson <nnc@lanl.gov.>
!! April 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gb_laser_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: gb_laser
    private
    real(r8) :: c1, c2, c3
  contains
    procedure :: init
    procedure :: irrad
  end type gb_laser
  
contains

  subroutine init(this, params)
    use parameter_list_type
    class(gb_laser), intent(out) :: this
    type(parameter_list) :: params
    real(r8) :: lambda, w0r, msq, p
    real(r8), parameter :: PI = 3.141592653589793_r8
    call params%get('wave-length', lambda)
    call params%get('waist-radius', w0r)
    call params%get('beam-param', msq)
    call params%get('power', p, default=1.0_r8)
    this%c1 = 0.5_r8 * w0r**2   ! length^2
    this%c2 = (pi * w0r**2) / (lambda * msq)  ! length
    this%c3 = p / PI
  end subroutine init

  pure function irrad(this, dx, dy, dz)
    class(gb_laser), intent(in) :: this
    real(r8), intent(in) :: dx, dy, dz
    real(r8) :: irrad, c4
    c4 = this%c1 * (1.0_r8 + (dz/this%c2)**2)
    irrad = this%c3 * exp(-(dx**2 + dy**2)/c4) / c4
  end function irrad

end module gb_laser_type
