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
  use scalars_module

  implicit none

  ! Private Module
  private

  ! Public Variables
  public :: BC_C, BC_P, BC_T, BC_V, Conc, Prs, Vel, &
            BC_STRUCTURE, BIT_LOCATION, BOUNDARY_CONDITION

  ! Public Subroutines

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define BC_STRUCTURE structure ..........................
  Type BC_STRUCTURE
     character(LEN = 1)                              :: ID

     integer(KIND = int_kind)                        :: Exist
     integer(KIND = int_kind)                        :: IFlag
     integer(KIND = int_kind)                        :: Shift
     integer(KIND = int_kind), pointer, dimension(:) :: Flag => null()

     integer,                  pointer,      dimension(:,:)   :: BCID => null()
     integer,                  pointer,      dimension(:,:)   :: BCOperator => null()

     logical(KIND = log_kind), pointer, dimension(:,:) :: UseFunc => null()

     real(KIND = real_kind), pointer, dimension(:,:)   :: Value1 => null()
     real(KIND = real_kind), pointer, dimension(:,:)   :: Value2 => null()
     real(KIND = real_kind), pointer, dimension(:,:)   :: Value3 => null()
     real(KIND = real_kind), pointer, dimension(:,:,:) :: Value_Array => null()
  End Type BC_STRUCTURE

  ! Declare BC_STRUCTURE structures
  type(BC_STRUCTURE), target, save :: BC_C ! Concentration
  type(BC_STRUCTURE), target, save :: BC_P ! Pressure
  type(BC_STRUCTURE), target, save :: BC_T ! Temperature
  type(BC_STRUCTURE), target, save :: BC_V ! Velocity

  ! Define BIT_LOCATION structure ..........................
  type BIT_LOCATION
     integer(KIND = int_kind), dimension(nfc) :: Face_bit
  end type BIT_LOCATION

  ! Declare BIT_LOCATION structures
  type(BIT_LOCATION) :: Prs  ! Pressure
  type(BIT_LOCATION) :: Vel  ! Velocity
  type(BIT_LOCATION) :: Conc ! Concentration

  ! Define BOUNDARY_CONDITION structure ....................
  type BOUNDARY_CONDITION
     integer(KIND = int_kind) :: Flag
     integer(KIND = int_kind) :: Internal
  end type BOUNDARY_CONDITION

END MODULE BC_TYPE_MODULE
