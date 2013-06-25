MODULE CUTOFFS_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define cutoff quantities.
  !
  ! Contains: None
  !
  ! Author(s): The Telluridians (telluride-info@lanl.gov)
  !
  !=======================================================================
  use kind_module, only: real_kind

  implicit none
  save

  ! Public Module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define cutoff quantities
  real(KIND = real_kind) :: alittle = EPSILON(1.0_REAL_KIND) ! Relative zero for REAL_KIND
                                                             ! Can be overwritten in Numerics namelist
  real(KIND = real_kind) :: cutvof  ! Volume fraction cutoff

END MODULE CUTOFFS_MODULE
