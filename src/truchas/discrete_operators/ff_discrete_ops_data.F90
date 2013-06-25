MODULE FF_DISCRETE_OPS_DATA
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables needed for discrete operators options
  !
  ! Contains: None
  !
  !
  !=======================================================================
  use kind_module,            only: log_kind
  use support_operators ,     only: SO_Control

  implicit none

  ! Private Module
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! NUMERICS namelist variable - set to 'default', 'ortho', 'nonortho', 'SO'
  character(LEN=80), public, save :: ff_discrete_ops_type

  ! Flag for enabling/disabling the use of ortho face gradients
  logical(KIND = log_kind), public, save :: use_ff_ortho_face_gradient, use_ff_support_operators
  
  ! This structure contains the mass matrix for the SO solve
  TYPE(SO_Control), public, POINTER, SAVE   :: FF_SO_Control_Data => Null()

END MODULE FF_DISCRETE_OPS_DATA
