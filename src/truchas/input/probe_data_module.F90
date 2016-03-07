!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use parameter_module, only: nprobes, MAX_PROBES
  use legacy_mesh_api, only: ndim
  implicit none    
  public    

  ! Probe Namelist Variables    
  character(80), dimension(0:MAX_PROBES), save :: probe_name, probe_description
  real(r8), dimension(ndim, 0:MAX_PROBES), save :: probe_coords
  real(r8), dimension(0:MAX_PROBES), save :: probe_coords_scale

END MODULE PROBE_DATA_MODULE
