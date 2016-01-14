!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module defines the abstract type AM_COORD_SYSTEM, which
!! which provides an interface for representing and manipulating the AM
!! laser-powder system.
!!
!! PROGRAMMING INTERFACE
!!
!!  The class AM_COORD_SYSTEM has the following type bound procedures:
!!
!!     INIT (THIS, PARAMS) 
!!      Initializes the AM_COORD_SYSTEM object THIS from parameters read in from the input file
!!      (encapsulated in the object parmas of type PARAMETER_LIST).
!!
!!   ADVANCE_AM_BLOCK_NOZZLE (THIS, DT)
!!      Updates the private data this%nozzle_x, this%nozzle_y, and this%nozzle_z of the 
!!      laser-powder system coordinates by advancing the coords.; here dt is the time increment
!!      used for advancing the laser-powder system coordinates.
!! 

module am_coord_system_class

  use kinds, only: r8
  use parameter_list_type
  implicit none
  private

  type, abstract, public :: am_coord_system
    real(r8) :: nozzle_x, nozzle_y, nozzle_z, &
                nozzle_speed
    real(r8) :: am_coords(3), travel_dir(3)
    contains
       procedure(init), deferred :: init
       procedure(advance_am_nozzle), deferred :: advance_am_nozzle
       procedure(print_obj), deferred :: print_obj
  end type am_coord_system


  abstract interface

    subroutine init (this, params)
      import am_coord_system, parameter_list
      class(am_coord_system), intent(inout) :: this
      type(parameter_list), intent(inout) :: params
    end subroutine

    subroutine advance_am_nozzle (this, dt)
      use kinds, only: r8
      import am_coord_system
      real(r8), intent(in) :: dt
      class(am_coord_system), intent(inout) :: this
    end subroutine

    subroutine print_obj (this)
      import am_coord_system
      class(am_coord_system), intent(inout) :: this
    end subroutine

  end interface


end module am_coord_system_class

