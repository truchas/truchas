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
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! NUMERICS namelist variable - set to 'default', 'ortho', 'nonortho'
  character(80), public, save :: discrete_ops_type

  ! Flag for enabling/disabling the use of ortho face gradients
  logical, public, save :: use_ortho_face_gradient

END MODULE DISCRETE_OPS_DATA
