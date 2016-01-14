!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_energy_data


  use kinds
  use truchas_logging_services
  use parameter_list_type
  implicit none
  private
  save

  public :: read_interface_energy_namelist

  !! parameter list
  type(parameter_list), public :: params_interface_energy

  real(r8) :: laser_power, laser_absorp, laser_sigma, &
              stefan_boltzmann, abszero


contains


  subroutine read_interface_energy_namelist (lun, interface_energy_type)
    
    integer, intent(in) :: lun
    character(len=32), intent(in) :: interface_energy_type

    if (interface_energy_type == "gaussian") then
      call read_interface_energy_gaussian_namelist (lun)
    else 
      call TLS_fatal ('interface_energy_type not supported')
    end if 


  end subroutine read_interface_energy_namelist


  subroutine read_interface_energy_gaussian_namelist (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use string_utilities
    
    integer, intent(in) :: lun

    integer :: ios
    logical :: found

    !! Magic values used to detect variables not initialized by input
  !  character, parameter :: NULL_C = char(0)
  !  integer,   parameter :: NULL_I = huge(1)
    real(r8),  parameter :: NULL_R = huge(1.0_r8)

    namelist /interface_energy/ laser_power, laser_absorp, &
                                laser_sigma, stefan_boltzmann, abszero

    !! Check that the manually-set module variable has a good value before reading.
  !  INSIST(am_enabled)
    
    call TLS_info ('')
    call TLS_info ('Reading INTERFACE_ENERGY namelist ...')

    !! Locate the INTERFACE_ENERGY namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'INTERFACE_ENERGY', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    call broadcast (found)
    if (.not.found) call TLS_fatal ('INTERFACE_ENERGY namelist not found')

    !! Read the namelist.
    if (is_IOP) then
      laser_power                = NULL_R 
      laser_absorp               = NULL_R 
      laser_sigma                = NULL_R
      stefan_boltzmann           = NULL_R
      abszero                    = NULL_R
      read(lun,nml=interface_energy,iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading INTERFACE_ENERGY namelist')

    !! Broadcast the namelist variables.
    call broadcast (laser_power)
    call broadcast (laser_absorp)
    call broadcast (laser_sigma)
    call broadcast (stefan_boltzmann)
    call broadcast (abszero)

    if (laser_power == NULL_R) then
        call TLS_fatal ('Need to specify the laser power: laser_power')
    end if


    if (laser_power == NULL_R) then
        call TLS_fatal ('Need to specify the absorption: laser_absorp')
    end if

    if (laser_sigma == NULL_R) then
        call TLS_fatal ('Need to specify the laser radius: laser_sigma')
    end if

    if ((stefan_boltzmann == NULL_R) .or. (abszero == NULL_R)) then
         stefan_boltzmann = 0
         abszero = 0
         call TLS_info ('  using default choice of no radiation: ' )
    end if

    ! Create parameter list
    params_interface_energy = &
      plist_interface_energy_gaussian(laser_power, laser_absorp, laser_sigma, &
                                      stefan_boltzmann, abszero)


  end subroutine read_interface_energy_gaussian_namelist




  function plist_interface_energy_gaussian(laser_power, laser_absorp, laser_sigma, &
                                           stefan_boltzmann, abszero)

    real(r8), intent(in) :: laser_power, laser_absorp, laser_sigma, &
                            stefan_boltzmann, abszero
    type(parameter_list) :: plist_interface_energy_gaussian

    call plist_interface_energy_gaussian%set ('interface_energy_type', 'gaussian')
    call plist_interface_energy_gaussian%set ('laser_power', laser_power)
    call plist_interface_energy_gaussian%set ('laser_absorp', laser_absorp)
    call plist_interface_energy_gaussian%set ('laser_sigma', laser_sigma)
    call plist_interface_energy_gaussian%set ('stefan_boltzmann', stefan_boltzmann)
    call plist_interface_energy_gaussian%set ('abszero', abszero)


  end function plist_interface_energy_gaussian





end module interface_energy_data
