!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! NUMERICS namelist variable - set to 'default', 'ortho', 'nonortho'
  character(80), public, save :: ff_discrete_ops_type

  ! Flag for enabling/disabling the use of ortho face gradients
  logical, public, save :: use_ff_ortho_face_gradient
  
END MODULE FF_DISCRETE_OPS_DATA
