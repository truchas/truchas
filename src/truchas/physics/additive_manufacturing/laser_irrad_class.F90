!!
!! LASER_IRRAD_CLASS
!!
!! This module provides the abstract base class LASER_IRRAD that defines a
!! common interface to models of laser irradiance.
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

module laser_irrad_class

  use parameter_list_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: laser_irrad
  contains
    procedure(init),  deferred :: init
    procedure(irrad), deferred :: irrad
  end type laser_irrad

  abstract interface
    subroutine init(this, params)
      import laser_irrad, parameter_list
      class(laser_irrad), intent(out) :: this
      type(parameter_list) :: params
    end subroutine
    pure function irrad(this, dx, dy, dz)
      import laser_irrad, r8
      class(laser_irrad), intent(in) :: this
      real(r8), intent(in) :: dx, dy, dz
      real(r8) :: irrad, c4
    end function
  end interface

end module laser_irrad_class
