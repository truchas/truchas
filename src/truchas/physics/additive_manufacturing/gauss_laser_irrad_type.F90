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
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "f90_assert.fpp"
module gauss_laser_irrad_type

  use laser_irrad_class
  use scalar_func_class
  use scalar_func_factories
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, extends(laser_irrad), public :: gauss_laser_irrad
    private
    real(r8) :: c1
    class(scalar_func), allocatable :: power
  contains
    procedure :: init
    procedure :: irrad
  end type gauss_laser_irrad

contains

  subroutine init(this, params)
    use parameter_list_type
    class(gauss_laser_irrad), intent(out) :: this
    type(parameter_list) :: params
    real(r8) :: sigma
    integer  :: stat
    character(:), allocatable :: errmsg

    call params%get('sigma', sigma)
    call alloc_scalar_func(params, 'power', this%power, stat, errmsg)
    INSIST(stat == 0)
    this%c1 = 2.0_r8 * sigma**2   ! length^2
  end subroutine init

  function irrad(this, t, dx, dy, dz)
    class(gauss_laser_irrad), intent(in) :: this
    real(r8), intent(in) :: t, dx, dy, dz
    real(r8), parameter :: PI = 3.141592653589793_r8
    real(r8) :: irrad
    irrad = this%power%eval([t]) / (PI * this%c1) * exp(-(dx**2 + dy**2)/this%c1)
  end function irrad

end module gauss_laser_irrad_type
