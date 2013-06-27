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
  use kinds, only: r8
  implicit none
  private

  ! Define cutoff quantities
  real(r8), public, save :: alittle = EPSILON(1.0_r8) ! Relative zero for real values
                                                      ! Can be overwritten in Numerics namelist
  real(r8), public, save :: cutvof  ! Volume fraction cutoff

END MODULE CUTOFFS_MODULE
