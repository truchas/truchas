MODULE LIMITER_DATA
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables used for slope or flux limiting.
  !
  ! Contains: None
  !
  ! Author(s): Ed Dendy     (dendy@lanl.gov)
  !            
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  real(r8), public,  dimension(:), allocatable :: phiUpMin   
  real(r8), public,  dimension(:), allocatable :: phiUpMax   
  real(r8), public,  dimension(:), allocatable :: phiMp1Min  
  real(r8), public,  dimension(:), allocatable :: phiMp1Max  
  real(r8), public,  dimension(:), allocatable :: sumVolIn   
  real(r8), public,  dimension(:), allocatable :: sumVolOut  
  real(r8), public,  dimension(:), allocatable :: sumPhiVolInMin
  real(r8), public,  dimension(:), allocatable :: sumPhiVolInMax
  real(r8), public,  dimension(:), allocatable :: PhiOutMin
  real(r8), public,  dimension(:), allocatable :: PhiOutMax

END MODULE LIMITER_DATA

