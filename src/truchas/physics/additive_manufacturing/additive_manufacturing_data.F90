!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module additive_manufacturing_data

  use kinds
  use truchas_logging_services
  use parameter_list_type 
  use interface_energy_data, only: read_interface_energy_namelist, &
                                   params_interface_energy
  use interface_mass_data, only: read_interface_mass_namelist, &
                                 params_interface_mass
  use am_coord_system_data, only: read_am_coord_system_namelist, &
                                  params_am_coord_system                                 
  implicit none
  private
  save

  public :: read_am_namelists

  !! This variable is defined manually prior to calling READ_DS_NAMELIST
  logical, public :: am_enabled = .false.

 ! character(len=32), public :: mesh_name


  !! parameter lists
  character(len=32), public :: interface_energy_type, &
                               am_coord_system_type, interface_mass_type
  type(parameter_list), public :: params_am 
  type(parameter_list), public, pointer :: p_params_interface_energy, &
                                           p_params_material_interface, &
                                           p_params_interface_mass, &
                                           p_params_am_coord_system



contains


  subroutine read_am_namelists (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use string_utilities
    
    integer, intent(in) :: lun

    integer :: ios
    logical :: found

    !! Magic values used to detect variables not initialized by input
    character, parameter :: NULL_C = char(0)

    namelist /additive_manufacturing/ interface_energy_type, &
                                      am_coord_system_type, interface_mass_type

    !! Check that the manually-set module variable has a good value before reading.
  !  INSIST(am_enabled)
    
    call TLS_info ('')
    call TLS_info ('Reading ADDITIVE_MANUFACTURING namelist ...')

    !! Locate the ADDITIVE_MANUFACTURING namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'ADDITIVE_MANUFACTURING', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    call broadcast (found)
    if (.not.found) call TLS_fatal ('ADDITIVE_MANUFACTURING namelist not found')

    !! Read the namelist.
    if (is_IOP) then
      interface_energy_type   = NULL_C
      am_coord_system_type    = NULL_C 
      interface_mass_type     = NULL_C
      read(lun,nml=additive_manufacturing,iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading ADDITIVE_MANUFACTURING namelist')

    !! Broadcast the namelist variables.
    call broadcast (interface_energy_type)
    call broadcast (am_coord_system_type)
    call broadcast (interface_mass_type)

    if (interface_energy_type == NULL_C) then
        interface_energy_type = "gaussian"
        call TLS_info ('  using default interface_energy_type: gaussian' )
    end if
    if (interface_mass_type == NULL_C) then
        interface_mass_type = "gaussian"
        call TLS_info ('  using default interface_mass_type: gaussian' )
    end if
    if (am_coord_system_type == NULL_C) then
        am_coord_system_type = "block_3D"
        call TLS_info ('  using default am_coord_system_type: block_3D' )
    end if

    ! Read in namelists for the different am components

    call read_am_coord_system_namelist (lun, am_coord_system_type)
    call read_interface_energy_namelist (lun, interface_energy_type)
    call read_interface_mass_namelist (lun, interface_mass_type)

   ! Create parameter lists

    p_params_am_coord_system => &
          params_am%sublist('am_coord_system')
    p_params_interface_energy => &
          params_am%sublist('interface_energy')
    p_params_interface_mass => &
          params_am%sublist('interface_mass')

    p_params_am_coord_system = params_am_coord_system
    p_params_interface_energy = params_interface_energy
    p_params_interface_mass = params_interface_mass


    !! REMOVE !!
   ! call test_am_paramater_list()
    !! REMOVE !!

  end subroutine read_am_namelists


  ! Tests that the am namelists were constructed correctly
  subroutine test_am_paramater_list()

    real(r8) :: nozzle_x, nozzle_y, nozzle_z, &
                laser_power, laser_absorp, &
                nozzle_speed, powder_rate, min_x, max_x, &
                laser_sigma, powder_sigma, laser_depth
    type(parameter_list), pointer :: p_params_test

    p_params_test => params_am%sublist('am_coord_system')
    call p_params_test%get('nozzle_x', nozzle_x)
    call p_params_test%get('nozzle_y', nozzle_y)
    call p_params_test%get('nozzle_z', nozzle_z)
    call p_params_test%get('nozzle_speed', nozzle_speed)
    call p_params_test%get('min_x', min_x)
    call p_params_test%get('max_x', max_x)

    print *, "nozzle_x: ", nozzle_x
    print *, "nozzle_y: ", nozzle_y
    print *, "nozzle_z: ", nozzle_z
    print *, "nozzle_speed: ", nozzle_speed
    print *, "min_x: ", min_x
    print *, "max_x: ", max_x

    p_params_test => params_am%sublist('interface_energy')
    call p_params_test%get('laser_absorp', laser_absorp)
    call p_params_test%get('laser_power', laser_power)
    call p_params_test%get('laser_sigma', laser_sigma)
    call p_params_test%get('laser_depth', laser_depth)
    print *, "laser_absorp: ", laser_absorp
    print *, "laser_sigma: ", laser_sigma
    print *, "laser_power: ", laser_power
    print *, "laser_depth: ", laser_depth

    p_params_test => params_am%sublist('interface_mass')
    call p_params_test%get('powder_rate', powder_rate)
    call p_params_test%get('powder_sigma', powder_sigma)
    print *, "powder_rate: ", powder_rate
    print *, "powder_sigma: ", powder_sigma


  end subroutine test_am_paramater_list


end module additive_manufacturing_data
