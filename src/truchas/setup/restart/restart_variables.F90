!!
!! RESTART_VARIABLES
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 29 Apr 2005; revised 28 July 2006.
!!
!! This module provides the variables associated with restarts, and a procedure
!! for reading those that belong to the RESTART namelist.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module restart_variables

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private
  save

  public :: read_restart_namelist

  !! These are defined externally by PROCESS_COMMAND_LINE.
  logical, public :: restart
  character(256), public :: restart_file

  !! The namelist variables.
  logical, public :: ignore_t  = .false.
  logical, public :: ignore_dt = .false.
  logical, public :: ignore_legacy_solid_mechanics = .false.
  logical, public :: ignore_joule_heat = .false.

  !! Common data read by OPEN_RESTART_FILE.
  integer, public :: restart_cycle_number, restart_ncells, restart_nnodes
  real(r8), public :: restart_t, restart_dt

  !! Optional restart file data segments; redefined by OPEN_RESTART_FILE.
  logical, public :: have_fluid_flow_data = .false.
  logical, public :: have_joule_heat_data = .false.
  logical, public :: have_solid_mechanics_data = .false.
  logical, public :: have_legacy_solid_mechanics_data = .false.
  logical, public :: have_species_data = .false.
  logical, public :: have_microstructure_data = .false.

contains

  subroutine read_restart_namelist (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_communication, only: is_IOP, broadcast

    integer, intent(in) :: lun

    namelist /restart/ ignore_t, ignore_dt, ignore_legacy_solid_mechanics, ignore_joule_heat

    logical :: found
    integer :: ios
    character(128) :: iom

    !! Locate the restart namelist (it's optional).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'RESTART', found)
    end if
    call broadcast (found)

    if (.not.found) return  ! fine, we'll just use the default values.

    !! Read the namelist.
    call input_info ('Reading RESTART namelist ...')
    if (is_IOP) read(lun,nml=restart,iostat=ios,iomsg=iom)
    call broadcast (ios)
    if (ios /= 0) call input_error ('error reading the RESTART namelist: ' // trim(iom))

    !! Broadcast the namelist variables to the other processors.
    call broadcast (ignore_t)
    call broadcast (ignore_dt)
    call broadcast (ignore_legacy_solid_mechanics)
    call broadcast (ignore_joule_heat)

  contains

    subroutine input_info (message)
      use truchas_logging_services
      character(len=*), intent(in) :: message
      call TLS_info (message)
    end subroutine input_info

    subroutine input_error (message)
      use truchas_logging_services
      character(len=*), intent(in) :: message
      call TLS_fatal ('READ_RESTART_NAMELIST: ' // message)
    end subroutine input_error

  end subroutine read_restart_namelist

end module restart_variables
