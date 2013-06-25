MODULE DISCRETE_OPS_DATA
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables needed for discrete operators options
  !
  ! Contains: None
  !
  !
  !=======================================================================
  use kind_module, only: log_kind

  implicit none

  ! Private Module
  private


  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! NUMERICS namelist variable - set to 'default', 'ortho', 'nonortho'
  character(LEN=80), public, save :: discrete_ops_type

  ! Flag for enabling/disabling the use of ortho face gradients
  logical(KIND = log_kind), public, save :: use_ortho_face_gradient

END MODULE DISCRETE_OPS_DATA
