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
  use kinds, only: r8
  use parameter_module, only: ndim, nprobes, MAX_PROBES
  implicit none    
  public    

  ! Probe Namelist Variables    
  character(80), dimension(0:MAX_PROBES), save :: probe_name, probe_description
  real(r8), dimension(ndim, 0:MAX_PROBES), save :: probe_coords
  real(r8), dimension(0:MAX_PROBES), save :: probe_coords_scale

END MODULE PROBE_DATA_MODULE
