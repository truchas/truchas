MODULE FLUID_TYPE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define and declare derived types for the fluid flow solution.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kind_module,      only: int_kind, real_kind
  use parameter_module, only: nfc

  implicit none

  ! Private Module
  private

  ! Public Variables
  public :: NORMS, DIV_NORMS, Div_c, Div_f

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define NORMS Structure
  type NORMS
     real(KIND = real_kind)   :: Linf, L1, L2
     integer(KIND = int_kind) :: Linf_Location
  end type NORMS

  ! Define DIV_NORMS Structure
  type DIV_NORMS
     type(NORMS) :: V_f
  end type DIV_NORMS

  ! Declare DIV_NORMS Structure
  type(DIV_NORMS), save                 :: Div_c
  type(DIV_NORMS), save, dimension(nfc) :: Div_f
  
END MODULE FLUID_TYPE_MODULE
