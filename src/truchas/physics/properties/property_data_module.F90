!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PROPERTY_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define material property variables and control parameters.
  !
  ! Contains: None
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Anand V. Reddy, LANL T-3 (jsbrock@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: maxcon, maxmat
  use legacy_mesh_api, only: ndim
  implicit none
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Material Namelist Variables
  integer, save :: background_material

  character(80), dimension(0:maxmat), save :: Material_Name, Material_Feature
  
  real(r8), dimension(0:maxmat), save :: density, Sound_Speed, Void_Temperature

  integer, dimension(0:maxmat), save :: Priority

  real(r8), dimension(ndim,0:maxmat), save :: Permeability_Constant

  ! Interface-tracking advection priority control parameters
  integer, dimension(0:maxmat), save :: Matpri

  ! Is the material not a fluid
  logical, save :: isImmobile(maxmat)

END MODULE PROPERTY_DATA_MODULE
