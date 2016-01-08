!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BODY_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Body force modeling variables and procedures
  !
  ! Author(s):  Jim Sicilian (sicilian@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: ndim

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Body Force Vector
  real(r8), dimension(ndim), public, save :: Body_Force

  ! Time-weighting for treatment of body forces
  real(r8), public, save :: body_force_implicitness

  ! Turn on face-based treatment of body-forces
  logical, public, save :: body_force_face_method

  ! potential + kinetic energy bounds
  real(r8), public, save :: mechanical_energy_bound

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

END MODULE BODY_DATA_MODULE
