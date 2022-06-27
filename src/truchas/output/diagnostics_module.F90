!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DIAGNOSTICS_MODULE
  !=======================================================================
  ! Purpose:
  !
  !    Data and procedures for diagnostic quantities for output
  !
  ! Contains: DIAGNOSTICS
  !           DIVERGENCE
  !
  ! Author(s): Kin Lam (klam@lanl.gov)
  !            Sharen Cummins (scummins@lanl.gov)
  !=======================================================================
  implicit none
  private

  ! public procedures
  public :: DIAGNOSTICS!, DIVERGENCE

  logical, save, public :: solidification_diagnostics = .false., &
                           alloy_solid_diagnostics    = .false.
  
CONTAINS

 !!
 !! NNC, 2 May 2012 -- As part of removing old HT, DIAGNOSTICS, which dealt with
 !! solidification times, has been gutted leaving just a stub.  We need to add
 !! back the appropriate analysis for the new heat transfer solver.
 !!

  SUBROUTINE DIAGNOSTICS ()
  
    logical, save :: first_time = .true.

    ! Decide whether to do various diagnostics calculations
    ! Create (and initialize) or nullify pointer arrays
    if (first_time) then
       ! If problem does not involve phase change, no solidification diagnostics
       solidification_diagnostics = .false.
       alloy_solid_diagnostics    = .false.
       if (solidification_diagnostics) then
          !TODO: Do something here based on the new phase change.
       end if
       first_time = .false.
       return   ! don't do diagnostics calculations the first time
    end if

    if (solidification_diagnostics) then
      !TODO: do something for the new phase change
    end if

  END SUBROUTINE DIAGNOSTICS

END MODULE DIAGNOSTICS_MODULE
