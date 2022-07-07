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
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_containers, only: scalar_func_box
  implicit none
  private

  integer, parameter, public :: mbody = 50
  integer, parameter, public :: mphi  =  5

  integer, public, save :: Matnum(mbody), nbody
  real(r8), public, save :: Body_Mass(mbody), Body_Vel(3,mbody), Body_Phi(mbody,mphi)
  type(scalar_func_box), public, save :: Body_Temp(mbody)

END MODULE INTERFACES_MODULE
