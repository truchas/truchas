!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vof_velocity_overwrite
  !=======================================================================
  ! Purpose(s):
  !
  !   Storage of information needed to overwrite velocity to run
  !   time varying canonical interface and advection tests.
  !
  ! Public Interface(s): None
  !
  ! Contains: None
  !
  ! Author(s): Robert Chiodi (robertchiodi@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only : string_len
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Logical telling if this was requested in the name list
  logical, public, save :: velocity_overwrite_requested

  ! Case used to identify velocity set in vtrack_driver.f90
  character(string_len), public, save :: velocity_overwrite_case

end module vof_velocity_overwrite
