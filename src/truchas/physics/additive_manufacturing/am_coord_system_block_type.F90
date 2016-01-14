!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module defines the derived type AM_COORD_SYSTEM_BLOCK, which
!! inherits from the abstract type AM_COORD_SYSTEM. The type
!! AM_COORD_SYSTEM_BLOCK is for representing and manipulating the AM
!! laser-powder system, assuming that the laser-powder system moves
!! long a straight line in the x-coordinate direction and between 
!! user-specified bounds min_x and max_x; upon reaching min_x or max_x,
!! the direction of the laser-powder system is reversed
!! 
!! 
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type AM_COORD_SYSTEM_BLOCK has the following type bound procedures:
!!
!!     INIT (THIS, PARAMS) 
!!      Initializes the AM_COORD_SYSTEM_BLOCK object THIS from parameters read in from the input file 
!!      (encapsulated in the object parmas of type PARAMETER_LIST).
!!
!!   ADVANCE_AM_BLOCK_NOZZLE (THIS, DT)
!!      Updates the private data this%nozzle_x, this%nozzle_y, and this%nozzle_z of the 
!!      laser-powder system coordinates by advancing the coords. by an increment dt*(1,0,0)
!! 

module am_coord_system_block_type

  use kinds, only: r8
  use am_coord_system_class
  implicit none
  private

  type, extends(am_coord_system), public :: am_coord_system_block
    real(r8) :: min_x, max_x
    contains
      procedure :: init => init_block
      procedure :: advance_am_nozzle => advance_am_block_nozzle
      procedure :: print_obj => print_obj
  end type am_coord_system_block


  contains


    subroutine init_block(this, params)

      use parameter_list_type

      class(am_coord_system_block), intent(inout) :: this
      type(parameter_list), intent(inout) :: params

      call params%get ('nozzle_x', this%nozzle_x)
      call params%get ('nozzle_y', this%nozzle_y)
      call params%get ('nozzle_z', this%nozzle_z)
      call params%get ('nozzle_speed', this%nozzle_speed)
      call params%get ('min_x', this%min_x)
      call params%get ('max_x', this%max_x)
      call params%get ('direction_x', this%travel_dir(1))
      call params%get ('direction_y', this%travel_dir(2))
      call params%get ('direction_z', this%travel_dir(3))

      ! Set the normalized velocity of the nozzle
    !  this%travel_dir(1:3) = 0
    !  this%travel_dir(1) = 1
      !!!! REMOVE, UNCOMMENT ABOVE !!!!
    !  this%travel_dir(2) = 1
      !!!! REMOVE, UNCOMMENT ABOVE !!!!

      ! Set current nozzle coordinates
      this%am_coords = (/ this%nozzle_x, this%nozzle_y, this%nozzle_z /)
      

    end subroutine init_block


    subroutine advance_am_block_nozzle(this, dt)

      real(r8), intent(in) :: dt
      class(am_coord_system_block), intent(inout) :: this


      ! Advance the laser-powder nozzle
   !   if ( this%am_coords(1) < this%min_x ) then
   !     this%travel_dir = -this%travel_dir
   !     this%am_coords(1) = this%min_x + .00000001_r8
   !   end if
   !   if ( this%am_coords(1) > this%max_x ) then
   !     this%travel_dir = -this%travel_dir
   !     this%am_coords(1) = this%max_x - .00000001_r8
   !   end if
   !   this%am_coords(1:3) = this%am_coords(1:3) + &
   !                          dt * this%nozzle_speed * this%travel_dir(1:3)
      !!!! REMOVE, UNCOMMENT ABOVE !!!!
      ! Advance the laser-powder nozzle
   !   if ( this%am_coords(2) < this%min_x ) then
   !     this%travel_dir = -this%travel_dir
   !     this%am_coords(2) = this%min_x + .00000001_r8
   !   end if
   !   if ( this%am_coords(2) > this%max_x ) then
   !     this%travel_dir = -this%travel_dir
   !     this%am_coords(2) = this%max_x - .00000001_r8
   !   end if
      this%am_coords(1:3) = this%am_coords(1:3) + &
                             dt * this%nozzle_speed * this%travel_dir(1:3)
      !!!! REMOVE, UNCOMMENT ABOVE !!!!

    end subroutine advance_am_block_nozzle



    subroutine print_obj(this)

      class(am_coord_system_block), intent(inout) :: this

      print *, "track bounds for am geom.: ", this%min_x, this%max_x
      print *, "init. nozzle coords for am geom.: ", this%nozzle_x, this%nozzle_y, this%nozzle_z  
      print *, "current nozzle coords for am geom.: ", this%am_coords(1:3)
      print *, "travel dir. for am geom.: ", this%travel_dir(1:3)

    end subroutine print_obj


end module am_coord_system_block_type
