!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PARAMETER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all integer and scalar parameters.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
  !
  !=======================================================================
  implicit none
  public
  
  ! Maximum number of outputs
  integer, parameter :: mops = 21

  ! Maximum output character string array size
  integer, parameter :: string_len = 256

  ! Interface (body) initialization parameters
  integer, parameter :: mbody = 50
  integer, parameter :: mphi  =  5

  ! Maximum number of material constants          
  integer, parameter :: max_slots = 10
  integer, parameter :: maxmat = 64
  integer :: nmat, mat_slot = 0, mat_slot_new = 0

  ! Maximum time intervals for electromagnetics
  integer, parameter, public :: MAXSV = 32

END MODULE PARAMETER_MODULE
