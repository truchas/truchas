!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CONSTANTS_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define commonly-used constants.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  integer,  parameter, public :: KILO = 1024
  integer,  parameter, public :: FLOATBYTES = 8

  real(r8), parameter, public :: pi  = 3.1415926535
  real(r8), parameter, public :: big = 1e10_r8

END MODULE CONSTANTS_MODULE
