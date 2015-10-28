!!
!! RESTART_VARIABLES
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 29 Apr 2005; revised 28 July 2006.
!!
!! This module provides the variables associated with restarts, and a procedure
!! for reading those that belong to the RESTART namelist.
!!

#include "f90_assert.fpp"

module restart_variables

  use kinds, only: r8
  use parameter_module, only: string_len
  implicit none
  private
  save

  public :: read_restart_namelist

  !! These are defined externally by PROCESS_COMMAND_LINE.
  logical, public :: restart
  character(len=string_len), public :: restart_file

  !! The namelist variables.
  logical, public :: ignore_t  = .false.
  logical, public :: ignore_dt = .false.
  logical, public :: ignore_solid_mechanics = .false.
  logical, public :: ignore_joule_heat = .false.

  !! Common data read by OPEN_RESTART_FILE.
  integer, public :: restart_cycle_number, restart_ncells, restart_nnodes
  real(r8), public :: restart_t, restart_dt

  !! Optional restart file data segments; redefined by OPEN_RESTART_FILE.
  logical, public :: have_fluid_flow_data = .false.
  logical, public :: have_joule_heat_data = .false.
  logical, public :: have_solid_mechanics_data = .false.
  logical, public :: have_species_data = .false.
  logical, public :: have_microstructure_data = .false.

contains

  subroutine read_restart_namelist (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_info_module, only: p_info
    use pgslib_module, only: pgslib_bcast
    
    integer, intent(in) :: lun

    namelist /restart/ ignore_t, ignore_dt, ignore_solid_mechanics, ignore_joule_heat

    logical :: found
    integer :: ios

    !! Locate the restart namelist (it's optional).
    if (p_info%IOP) then
      rewind lun
      call seek_to_namelist (lun, 'RESTART', found)
    end if
    call pgslib_bcast (found)

    if (.not.found) return  ! fine, we'll just use the default values.

    !! Read the namelist.
    call input_info ('Reading the RESTART namelist ...')
    if (p_info%IOP) read(lun,nml=restart,iostat=ios)
    call pgslib_bcast (ios)
    if (ios /= 0) call input_error ('Error reading the RESTART namelist')

    !! Broadcast the namelist variables to the other processors.
    call pgslib_bcast (ignore_t)
    call pgslib_bcast (ignore_dt)
    call pgslib_bcast (ignore_solid_mechanics)
    call pgslib_bcast (ignore_joule_heat)

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
