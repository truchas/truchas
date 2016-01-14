!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module defines the abstract type INTERFACE_ENERGY, which
!! which provides an interface for the energy deposition process
!!
!! PROGRAMMING INTERFACE
!!
!!  The class INTERFACE_ENERGY has the following type bound procedures:
!!
!!     ALLOC (THIS) 
!!      Allocates whatever arrays are needed for the INTERFACE_ENERGY object THIS
!!
!!     INIT (THIS, PARAMS) 
!!      Initializes the INTERFACE_ENERGY object THIS from parameters read in from the input file
!!      (encapsulated in the object parmas of type PARAMETER_LIST).
!!
!!   DEPOSIT_ENERGY (THIS, STATE, AM_GEOM, DT)
!!      Updates the state data (e.g. temperatures and enthalpies) stored in the object STATE 
!!      (of type AM_STATE). This procedure uses as input the object AM_GEOM
!!       (of some type derived from class AM_COORD_SYSTEM) and the time step DT



module interface_energy_class

  use kinds, only: r8
  use parameter_list_type
  use am_state_type
  use am_coord_system_class
  implicit none
  private

  type, abstract, public :: interface_energy
    contains
      procedure(init), deferred :: init
      procedure(alloc), deferred :: alloc
      procedure(deposit_energy), deferred :: deposit_energy
      procedure(print_obj), deferred :: print_obj
    end type interface_energy


    abstract interface

      subroutine init (this, params)
        import interface_energy, parameter_list
        class(interface_energy), intent(inout) :: this
        type(parameter_list), intent(inout) :: params
      end subroutine

      subroutine alloc (this)
        import interface_energy
        class(interface_energy), intent(inout) :: this
      end subroutine

      subroutine deposit_energy(this, state, am_geom, dt)
        use kinds, only: r8
        import interface_energy, am_state, am_coord_system
        class(interface_energy), intent(inout) :: this
        type(am_state), intent(inout) :: state
        class(am_coord_system), intent(in) :: am_geom
        real(r8), intent(in) :: dt
      end subroutine deposit_energy

      subroutine print_obj(this)
        import interface_energy
        class(interface_energy), intent(inout) :: this
      end subroutine print_obj

    end interface


end module interface_energy_class

