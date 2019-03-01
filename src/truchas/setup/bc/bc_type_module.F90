!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BC_TYPE_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define derived types (structures) associated with boundary conditions.
  !
  ! Contains: None
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Bryan R. Lally, LANL ESA-EPE (lally@lanl.gov)
  !
  !=======================================================================
  use legacy_mesh_api, only: nfc
  use kinds, only: r8

  implicit none
  private

  ! Public Variables
  public :: BC_C, BC_P, BC_T, BC_V, Conc, Prs, Vel, &
            BC_STRUCTURE, BIT_LOCATION, BOUNDARY_CONDITION

  ! Public Subroutines

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define BC_STRUCTURE structure ..........................
  Type BC_STRUCTURE
     character :: ID

     integer :: Exist
     integer :: IFlag
     integer :: Shift
     integer, pointer, dimension(:) :: Flag => null()

     integer, pointer, dimension(:,:)   :: BCID => null()
     integer, pointer, dimension(:,:)   :: BCOperator => null()

     logical, pointer, dimension(:,:) :: UseFunc => null()

     real(r8), pointer, dimension(:,:)   :: Value1 => null()
     real(r8), pointer, dimension(:,:)   :: Value2 => null()
     real(r8), pointer, dimension(:,:)   :: Value3 => null()
     real(r8), pointer, dimension(:,:,:) :: Value_Array => null()
  End Type BC_STRUCTURE

  ! Declare BC_STRUCTURE structures
  type(BC_STRUCTURE), target, save :: BC_C ! Concentration
  type(BC_STRUCTURE), target, save :: BC_P ! Pressure
  type(BC_STRUCTURE), target, save :: BC_T ! Temperature
  type(BC_STRUCTURE), target, save :: BC_V ! Velocity

  ! Define BIT_LOCATION structure ..........................
  type BIT_LOCATION
     integer, dimension(nfc) :: Face_bit
  end type BIT_LOCATION

  ! Declare BIT_LOCATION structures
  type(BIT_LOCATION) :: Prs  ! Pressure
  type(BIT_LOCATION) :: Vel  ! Velocity
  type(BIT_LOCATION) :: Conc ! Concentration

  ! Define BOUNDARY_CONDITION structure ....................
  type BOUNDARY_CONDITION
     integer :: Flag
     integer :: Internal
  end type BOUNDARY_CONDITION

END MODULE BC_TYPE_MODULE
