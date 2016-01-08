!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
