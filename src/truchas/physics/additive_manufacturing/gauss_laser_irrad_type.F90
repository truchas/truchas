!!
!! GAUSS_LASER_IRRAD_TYPE
!!
!! Laser irradiance function with a constant Gaussian profile along the beam axis.
!!
!! Neil Carlson <nnc@lanl.gov>
!! March 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gauss_laser_irrad_type

  use laser_irrad_class
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, extends(laser_irrad), public :: gauss_laser_irrad
    private
    real(r8) :: c1, c2
  contains
    procedure :: init
    procedure :: irrad
  end type gauss_laser_irrad

contains

  subroutine init(this, params)
    use parameter_list_type
    class(gauss_laser_irrad), intent(out) :: this
    type(parameter_list) :: params
    real(r8) :: sigma, p
    real(r8), parameter :: PI = 3.141592653589793_r8
    call params%get('sigma', sigma)
    call params%get('power', p, default=1.0_r8)
    this%c1 = 0.5_r8 * sigma**2   ! length^2
    this%c2 = p / (PI * this%c1)
  end subroutine init

  pure function irrad(this, dx, dy, dz)
    class(gauss_laser_irrad), intent(in) :: this
    real(r8), intent(in) :: dx, dy, dz
    real(r8) :: irrad
    irrad = this%c2 * exp(-(dx**2 + dy**2)/this%c1)
  end function irrad

end module gauss_laser_irrad_type
