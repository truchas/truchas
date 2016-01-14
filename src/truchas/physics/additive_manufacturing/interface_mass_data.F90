!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_mass_data



  use kinds
  use truchas_logging_services
  use parameter_list_type
  implicit none
  private
  save

  public :: read_interface_mass_namelist

  !! parameter list
  type(parameter_list), public :: params_interface_mass

  real(r8) :: powder_rate, powder_sigma



contains


  subroutine read_interface_mass_namelist (lun, interface_mass_type)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use string_utilities
    
    integer, intent(in) :: lun
    character(len=32), intent(in) :: interface_mass_type

    integer :: ios
    logical :: found, exists

    !! Magic values used to detect variables not initialized by input
    character, parameter :: NULL_C = char(0)

    if (interface_mass_type == "gaussian") then
      call read_interface_mass_gaussian_namelist (lun)
    else 
      call TLS_fatal ('interface_mass_type not supported')
    end if 


  end subroutine read_interface_mass_namelist



  subroutine read_interface_mass_gaussian_namelist (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use string_utilities
    
    integer, intent(in) :: lun

    integer :: ios
    logical :: found, exists

    !! Magic values used to detect variables not initialized by input
    character, parameter :: NULL_C = char(0)
    integer,   parameter :: NULL_I = huge(1)
    real(r8),  parameter :: NULL_R = huge(1.0_r8)

    namelist /interface_mass/ powder_rate, powder_sigma
    
    call TLS_info ('')
    call TLS_info ('Reading INTERFACE_MASS namelist ...')

    !! Locate the ADDITIVE_MANUFACTURING namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'INTERFACE_MASS', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    call broadcast (found)
    if (.not.found) call TLS_fatal ('INTERFACE_MASS namelist not found')

    !! Read the namelist.
    if (is_IOP) then
      powder_rate         = NULL_R 
      powder_sigma        = NULL_R 
      read(lun,nml=interface_mass,iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading INTERFACE_MASS namelist')

    !! Broadcast the namelist variables.
    call broadcast (powder_rate)
    call broadcast (powder_sigma)

    if (powder_rate == NULL_R) then
        call TLS_fatal ('Need to specify the powder deposition rate: powder_rate')
    end if

    if (powder_sigma == NULL_R) then
        call TLS_fatal ('Need to specify the powder deposition radius: powder_sigma')
    end if

    ! Create parameter list
    params_interface_mass = &
      plist_interface_mass_gaussian(powder_rate, powder_sigma)


  end subroutine read_interface_mass_gaussian_namelist



  function plist_interface_mass_gaussian(powder_rate, powder_sigma)

    real(r8), intent(in) :: powder_rate, powder_sigma
    type(parameter_list) :: plist_interface_mass_gaussian

    call plist_interface_mass_gaussian%set ('interface_mass_type', 'gaussian')
    call plist_interface_mass_gaussian%set ('powder_rate', powder_rate)
    call plist_interface_mass_gaussian%set ('powder_sigma', powder_sigma)


  end function plist_interface_mass_gaussian




end module interface_mass_data
