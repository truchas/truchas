!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
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
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  ! Define cutoff quantities
  real(r8), public, save :: cutvof  ! Volume fraction cutoff

END MODULE CUTOFFS_MODULE
