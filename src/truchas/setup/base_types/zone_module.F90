!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module zone_module

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: zone_init, zone_free, read_zone_data

  ! physical state variables on a cell
  type, public :: cell_avg
     real(r8) :: rho       ! current cell average density
     real(r8) :: rho_old   ! past cell average density
     real(r8) :: temp      ! current temperature
     real(r8) :: temp_old  ! past temperature
     real(r8) :: enthalpy      ! current enthalpy
     real(r8) :: enthalpy_old  ! past enthalpy
     real(r8) :: p         ! current pressure
     real(r8) :: vc(3)     ! current cell-centered velocity
     real(r8) :: vc_old(3) ! past cell-centered velocity
  end type cell_avg

  type(cell_avg), allocatable, target, public :: zone(:)

CONTAINS

  subroutine zone_init(ncell)
    integer, intent(in) :: ncell
    integer :: n
    allocate(zone(ncell))
    zone%rho          = 0.0_r8
    zone%rho_old      = 0.0_r8
    zone%temp         = 0.0_r8
    zone%temp_old     = 0.0_r8
    zone%enthalpy     = 0.0_r8
    zone%enthalpy_old = 0.0_r8
    zone%p            = 0.0_r8
    do n = 1, 3
      zone%vc(n)      = 0.0_r8
      zone%vc_old(n)  = 0.0_r8
    end do
  end subroutine zone_init

  subroutine zone_free
    if (allocated(zone)) deallocate(zone)
  end subroutine

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

    do n = 1, 3
      call read_dist_array (unit, zone%vc(n), pcell, 'READ_ZONE_DATA: error reading VC records')
      zone%vc_old(n) = zone%vc(n)
    end do

  end subroutine read_zone_data

end module zone_module
