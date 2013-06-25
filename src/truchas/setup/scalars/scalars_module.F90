MODULE SCALARS_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Primary repository module for most of the principal code
  !   scalars. This module does nothing more than serve as a
  !   principal repository for other scalar modules, simply by
  !   "using" certain scalar modules. This allows only one use
  !   statement ("use scalars_module") for scalars per routine.
  !
  ! Contains: None
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use code_module
  use constants_module
  use cutoffs_module
  use kind_module
  use parameter_module

  implicit none
  save

  ! Public Module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

END MODULE SCALARS_MODULE
