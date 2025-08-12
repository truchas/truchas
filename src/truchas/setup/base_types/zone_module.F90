!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module zone_module

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: zone_init, read_zone_data

  type :: zone_type
     real(r8), allocatable :: rho(:)       ! current cell average density
     real(r8), allocatable :: rho_old(:)   ! past cell average density
     real(r8), allocatable :: temp(:)      ! current temperature
     real(r8), allocatable :: temp_old(:)  ! past temperature
     real(r8), allocatable :: enthalpy(:)      ! current enthalpy
     real(r8), allocatable :: enthalpy_old(:)  ! past enthalpy
     real(r8), allocatable :: p(:)         ! current pressure
     real(r8), allocatable :: vc(:,:)      ! current cell-centered velocity
     real(r8), allocatable :: vc_old(:,:)  ! past cell-centered velocity
  end type
  type(zone_type), target, public :: zone

CONTAINS

  subroutine zone_init(ncell)
    integer, intent(in) :: ncell
    allocate(zone%rho(ncell), source=0.0_r8)
    allocate(zone%rho_old(ncell), source=0.0_r8)
    allocate(zone%temp(ncell), source=0.0_r8)
    allocate(zone%temp_old(ncell), source=0.0_r8)
    allocate(zone%enthalpy(ncell), source=0.0_r8)
    allocate(zone%enthalpy_old(ncell), source=0.0_r8)
    allocate(zone%p(ncell), source=0.0_r8)
    allocate(zone%vc(3,ncell), source=0.0_r8)
    allocate(zone%vc_old(3,ncell), source=0.0_r8)
  end subroutine zone_init

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

    use restart_utilities, only: read_dist_array
    use unstr_mesh_type
    use mesh_manager, only: unstr_mesh_ptr

    integer, intent(in) :: unit, version
    type(unstr_mesh), pointer :: mesh

    integer :: n

    mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(mesh))

    associate (pcell => mesh%xcell(1:mesh%ncell_onP))
      call read_dist_array (unit, zone%rho,      pcell, 'READ_ZONE_DATA: error reading RHO record')
      call read_dist_array (unit, zone%temp,     pcell, 'READ_ZONE_DATA: error reading TEMP record')
      call read_dist_array (unit, zone%enthalpy, pcell, 'READ_ZONE_DATA: error reading ENTHALPY record')
      call read_dist_array (unit, zone%p,        pcell, 'READ_ZONE_DATA: error reading PRESSURE record')

      zone%rho_old      = zone%rho
      zone%temp_old     = zone%temp
      zone%enthalpy_old = zone%enthalpy

      do n = 1, 3
        call read_dist_array (unit, zone%vc(n,:), pcell, 'READ_ZONE_DATA: error reading VC records')
        zone%vc_old(n,:) = zone%vc(n,:)
      end do
    end associate

  end subroutine read_zone_data

end module zone_module
