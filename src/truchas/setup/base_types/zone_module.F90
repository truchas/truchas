!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ZONE_MODULE
  !=======================================================================
  ! Purpose:
  !    Define the CELL_AVG derived type
  !
  ! Contains:
  !    STREAM_ZONE_OUT (Stream, Zone)
  !    ZONE_COLLATE (Zone)
  !    STREAM_ZONE_IN (Stream, Zone)
  !    ZONE_DISTRIBUTE (Zone_Tot)
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !=======================================================================
  use kinds, only: r8
  use legacy_mesh_api, only: ndim
  implicit none
  private

  ! public variables and types
  public :: CELL_AVG, Zone, read_zone_data

  !-----------------------------------------------------------------------------

  ! CELL_AVG structure - hold physical state variables
  type CELL_AVG
     real(r8)                      :: Rho       ! current cell average density
     real(r8)                      :: Rho_Old   ! past cell average density
     real(r8)                      :: Temp      ! current temperature
     real(r8)                      :: Temp_Old  ! past temperature
     real(r8)                      :: Enthalpy      ! current enthalpy
     real(r8)                      :: Enthalpy_Old  ! past enthalpy
     real(r8)                      :: P         ! current pressure
     real(r8), dimension(ndim)     :: Vc        ! current cell-centered velocity
     real(r8), dimension(ndim)     :: Vc_old    ! past cell-centered velocity
  end type CELL_AVG

  type(CELL_AVG), dimension(:), pointer :: Zone

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_ZONE_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 6 Jan 2005, revised 21 Apr 2005
 !!
 !! This subroutine reads the zone data from a restart file opened (and pre-
 !! positioned) on UNIT, and initializes the module structure ZONE with this
 !! data (properly distributed and permuted).  VERSION is the version number
 !! of the restart file format.
 !!
 !! NB: We have no data to initialize the *_OLD components of the zone
 !! structure, so we simply set those to be the current values.  We wouldn't
 !! normally need to do this (they are overwritten at the start of a time
 !! step) except the output routine wants to write them at the start :-(
 !!

  subroutine read_zone_data (unit, version)

    use legacy_mesh_api, only: pcell => unpermute_mesh_vector
    use restart_utilities, only: read_dist_array

    integer, intent(in) :: unit, version

    integer :: n

    call read_dist_array (unit, zone%rho,      pcell, 'READ_ZONE_DATA: error reading RHO record')
    call read_dist_array (unit, zone%temp,     pcell, 'READ_ZONE_DATA: error reading TEMP record')
    call read_dist_array (unit, zone%enthalpy, pcell, 'READ_ZONE_DATA: error reading ENTHALPY record')
    call read_dist_array (unit, zone%p,        pcell, 'READ_ZONE_DATA: error reading PRESSURE record')

    zone%rho_old      = zone%rho
    zone%temp_old     = zone%temp
    zone%enthalpy_old = zone%enthalpy

    do n = 1, ndim
      call read_dist_array (unit, zone%vc(n), pcell, 'READ_ZONE_DATA: error reading VC records')
      zone%vc_old(n) = zone%vc(n)
    end do

  end subroutine read_zone_data

END MODULE ZONE_MODULE
