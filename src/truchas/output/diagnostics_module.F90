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
  use kinds, only: r8
  implicit none
  private

  ! public procedures
  public :: DIAGNOSTICS, DIVERGENCE

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


  SUBROUTINE DIVERGENCE (D)

    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell-centered divergence (D) for diagnostic purposes
    !
    !=======================================================================
    use fluid_data_module, only: fluidRho
    use fluid_data_module, only: Fluxing_Velocity
    use legacy_mesh_api,   only: nfc, ncells, Cell 
    use time_step_module,  only: dt

    ! Arguments
    real(r8), dimension(ncells), intent(OUT) :: D

    ! Local Variables
    integer :: f 

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over faces and accumulate the divergence
    ! from these face-centered velocities
    D = 0.0_r8
    do f = 1,nfc
       D = D + Fluxing_Velocity(f,:)*Cell%Face_Area(f)
    end do

    ! Normalize by cell volume
    D = D*dt/Cell%Volume

    ! Zero out divergences in void cells
    where (FluidRho == 0.0_r8) D = 0.0_r8

    ! Eliminate noise
!   D = MERGE(0.0_r8, D, ABS(D) <= alittle)

  END SUBROUTINE DIVERGENCE

END MODULE DIAGNOSTICS_MODULE
