MODULE PROBE_DATA_MODULE    
!=======================================================================   
! Purpose(s):   !   Define PROBE variables to be read in from input file.   
!                   Note: This is a public module.   
!   
! Contains: None   
!   
! Author(s): Sharen Cummins (scummins@lanl.gov)   
!   
!=======================================================================   
use kind_module,      only: int_kind, log_kind, real_kind   
use parameter_module, only: ndim, nprobes, MAX_PROBES
    
implicit none    
! Public Module   
public    

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>    
! Probe Namelist Variables    
character(LEN = 80), dimension(0:MAX_PROBES), save :: &
                        probe_name,                   &
                        probe_description

real(KIND = real_kind), dimension(ndim, 0:MAX_PROBES), save :: &
                        probe_coords

real(KIND = real_kind), dimension(0:MAX_PROBES), save :: &
                        probe_coords_scale


END MODULE PROBE_DATA_MODULE
