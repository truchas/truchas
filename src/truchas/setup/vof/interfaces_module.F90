!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE INTERFACES_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables necessary for the generation of initial interface
  !   locations.
  !
  ! Public Interface(s): None
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: mbody, mphi
  use legacy_mesh_api, only: ndim
  use scalar_func_containers, only: scalar_func_box
  implicit none
  private

  integer, public, save :: Matnum(mbody), nbody
  real(r8), public, save :: Body_Mass(mbody), Body_Vel(ndim,mbody), Body_Phi(mbody,mphi)
  type(scalar_func_box), public, save :: Body_Temp(mbody)

END MODULE INTERFACES_MODULE
