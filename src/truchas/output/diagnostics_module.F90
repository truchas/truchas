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

  ! The get_global_cell procedure is not used any more in Truchas and 
  ! is a candidate for removal.

  subroutine get_global_cell (coordinates, icell_o, multiplicity)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the global cell number that is closest to a given coordinate set
    !
    !=======================================================================

    use legacy_mesh_api,      only: ndim, ncells, Cell
    use PGSLIB_module,        only: PGSLib_Global_SUM, pgslib_global_minval

    ! Argument List

    integer, intent(out)  :: icell_o
    real(r8), intent(out) :: multiplicity
    real(r8), intent(in)  :: coordinates(ndim)

    ! Local Variables

    real(r8) :: tmp(ncells)
    real(r8) :: Minimum_distance 
    real(r8) :: Minimum_distance_global
    real(r8) :: Minimum_distance_global_tot
    integer  :: icell(1)
    
    ! Find cell at minimum distance from the point

    tmp(:) = (coordinates(1) - Cell(:)%Centroid(1))**2 +        &
         (coordinates(2) - Cell(:)%Centroid(2))**2 +            &
         (coordinates(3) - Cell(:)%Centroid(3))**2 

    Minimum_distance = minval(tmp)
    icell            = minloc(tmp)

    Minimum_distance_global = pgslib_global_minval(tmp)

    if ( Minimum_distance_global == Minimum_distance) then
       icell_o          = icell(1)
    else
       icell_o          = -1
       Minimum_distance = 0.0_r8
    endif

    if (Minimum_distance_global == 0) then
       multiplicity                = 1.0_r8
    else
       Minimum_distance_global_tot = pgslib_global_sum(Minimum_distance)
       multiplicity                = Minimum_distance_global_tot / Minimum_distance_global 
    endif

  END SUBROUTINE get_global_cell 

  ! The get_global_node procedure is not used any more in Truchas and 
  ! is a candidate for removal.

  subroutine get_global_node (coordinates, inode_o, multiplicity)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the global node number that is closest to a given coordinate set
    !
    !=======================================================================

    use legacy_mesh_api,      only: ndim, nnodes, Vertex
    use PGSLIB_module,        only: PGSLib_Global_SUM, pgslib_global_minval

    ! Argument List

    integer, intent(out)  :: inode_o
    real(r8), intent(out) :: multiplicity
    real(r8), intent(in)  :: coordinates(ndim)

    ! Local Variables

    real(r8) :: tmp(nnodes)
    real(r8) :: Minimum_distance 
    real(r8) :: Minimum_distance_global
    real(r8) :: Minimum_distance_global_tot
    integer  :: inode(1)
    
    ! Find node at minimum distance from the point

    tmp(:) = (coordinates(1) - Vertex(:)%Coord(1))**2 +        &
         (coordinates(2) - Vertex(:)%Coord(2))**2 +            &
         (coordinates(3) - Vertex(:)%Coord(3))**2 

    Minimum_distance = minval(tmp)
    inode            = minloc(tmp)

    Minimum_distance_global = pgslib_global_minval(tmp)

    if ( Minimum_distance_global == Minimum_distance) then
       inode_o          = inode(1)
    else
       inode_o          = -1
       Minimum_distance = 0.0_r8
    endif

    if (Minimum_distance_global == 0) then
       multiplicity     = 1.0_r8
    else
       Minimum_distance_global_tot = pgslib_global_sum(Minimum_distance)
       multiplicity                = Minimum_distance_global_tot / Minimum_distance_global 
    endif

  END SUBROUTINE get_global_node

END MODULE DIAGNOSTICS_MODULE
