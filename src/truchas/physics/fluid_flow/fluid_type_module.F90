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
  use kinds, only: r8
  use parameter_module, only: nfc
  implicit none
  private

  ! Public Variables
  public :: NORMS, DIV_NORMS, Div_c, Div_f

  ! Define NORMS Structure
  type NORMS
     real(r8) :: Linf, L1, L2
     integer  :: Linf_Location
  end type NORMS

  ! Define DIV_NORMS Structure
  type DIV_NORMS
     type(NORMS) :: V_f
  end type DIV_NORMS

  ! Declare DIV_NORMS Structure
  type(DIV_NORMS), save                 :: Div_c
  type(DIV_NORMS), save, dimension(nfc) :: Div_f
  
END MODULE FLUID_TYPE_MODULE
