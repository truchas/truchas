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
  use kind_module,        only: real_kind

  implicit none

  ! Private Module
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    real(real_kind), public,  dimension(:),         allocatable :: phiUpMin   
    real(real_kind), public,  dimension(:),         allocatable :: phiUpMax   
    real(real_kind), public,  dimension(:),         allocatable :: phiMp1Min  
    real(real_kind), public,  dimension(:),         allocatable :: phiMp1Max  
    real(real_kind), public,  dimension(:),         allocatable :: sumVolIn   
    real(real_kind), public,  dimension(:),         allocatable :: sumVolOut  
    real(real_kind), public,  dimension(:),         allocatable :: sumPhiVolInMin
    real(real_kind), public,  dimension(:),         allocatable :: sumPhiVolInMax
    real(real_kind), public,  dimension(:),         allocatable :: PhiOutMin
    real(real_kind), public,  dimension(:),         allocatable :: PhiOutMax

END MODULE LIMITER_DATA

