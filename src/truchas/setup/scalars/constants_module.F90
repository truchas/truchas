!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
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

  real(r8), parameter, public :: pi  = 3.14159265358979323846_r8
  real(r8), parameter, public :: big = 1e10_r8

END MODULE CONSTANTS_MODULE
