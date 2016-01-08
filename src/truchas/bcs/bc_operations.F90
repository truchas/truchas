!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_Operations
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the operations to support boundary condition operations.
  !   This module provides the application visible routines and
  !   the application visible data.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@lanl.gov)
  !-----------------------------------------------------------------------------
  use bc_data_types
  Implicit None

  ! Data instances provided by this module.
  type (BC_Specifier), SAVE, target :: Temperature_BC
  type (BC_Specifier), SAVE, target :: Pressure_BC

END Module BC_Operations
