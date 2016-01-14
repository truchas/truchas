!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module am_coord_system_data


  use kinds
  use truchas_logging_services
  use parameter_list_type
  implicit none
  private
  save


  public :: read_am_coord_system_namelist

  !! parameter list
  type(parameter_list), public :: params_am_coord_system

  real(r8) :: nozzle_x, nozzle_y, nozzle_z, nozzle_speed, min_x, max_x, &
              coordinate_scale_factor


contains


  subroutine read_am_coord_system_namelist (lun, am_coord_system_type)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use string_utilities
    
    integer, intent(in) :: lun
    character(len=32), intent(in) :: am_coord_system_type

    if (am_coord_system_type=="block_3D") then
      call read_am_coord_system_block_namelist (lun)
    else 
      call TLS_fatal ('am_coord_system_type not supported')
    end if 


  end subroutine read_am_coord_system_namelist


  subroutine read_am_coord_system_block_namelist (lun)

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

    namelist /am_coord_system/ nozzle_x, nozzle_y, &
                               nozzle_z, nozzle_speed, min_x, max_x, &
                               coordinate_scale_factor


    !! Check that the manually-set module variable has a good value before reading.
  !  INSIST(am_enabled)
    
    call TLS_info ('')
    call TLS_info ('Reading AM_COORD_SYSTEM namelist ...')

    !! Locate the ADDITIVE_MANUFACTURING namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'AM_COORD_SYSTEM', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    call broadcast (found)
    if (.not.found) call TLS_fatal ('AM_COORD_SYSTEM namelist not found')


    !! Read the namelist.
    if (is_IOP) then
      nozzle_x                = NULL_R
      nozzle_y                = NULL_R
      nozzle_z                = NULL_R
      nozzle_speed            = NULL_R
      min_x                   = NULL_R
      max_x                   = NULL_R
      coordinate_scale_factor = NULL_R
      read(lun,nml=am_coord_system,iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading AM_COORD_SYSTEM namelist')

    !! Broadcast the namelist variables.
    call broadcast (nozzle_x)
    call broadcast (nozzle_y)
    call broadcast (nozzle_z)
    call broadcast (nozzle_speed)
    call broadcast (min_x)
    call broadcast (max_x)
    call broadcast (coordinate_scale_factor)

    if (nozzle_x == NULL_R .or. nozzle_y == NULL_R .or. nozzle_z == NULL_R) then
        call TLS_fatal ('Need to specify the initial nozzle coordinates: nozzle_x, nozzle_y, nozzle_z')
    end if

    if (nozzle_speed == NULL_R) then
        call TLS_fatal ('Need to specify the nozzle speed: nozzle_speed')
    end if

    if (min_x == NULL_R) then
         min_x = nozzle_x
         call TLS_info ('  using default min. x value: ' )
    end if

    if (max_x == NULL_R) then
         max_x = -nozzle_x
         call TLS_info ('  using default max. x value: ' )
    end if

    if (coordinate_scale_factor == NULL_R) then
         coordinate_scale_factor = 1.0_r8
         call TLS_info ('  using default coordinate_scale_factor: 1.0' )
    end if

    ! scale spatial coordinates by coordinate_scale_factor
    nozzle_x = nozzle_x * coordinate_scale_factor
    nozzle_y = nozzle_y * coordinate_scale_factor
    nozzle_z = nozzle_z * coordinate_scale_factor
   ! nozzle_speed = nozzle_speed * coordinate_scale_factor
    min_x = min_x * coordinate_scale_factor
    max_x = max_x * coordinate_scale_factor

    ! create paramater list
    params_am_coord_system = &
        plist_am_coord_system_block(nozzle_x, nozzle_y, &
                                    nozzle_z, nozzle_speed, &
                                    min_x, max_x)


  end subroutine read_am_coord_system_block_namelist



  function plist_am_coord_system_block(nozzle_x, nozzle_y, &
                                       nozzle_z, nozzle_speed, min_x, max_x)

    real(r8), intent(in) :: nozzle_x, nozzle_y, nozzle_z, &
                            nozzle_speed, min_x, max_x
    type(parameter_list) :: plist_am_coord_system_block

    call plist_am_coord_system_block%set ('am_coord_system_type', 'block_3D')
    call plist_am_coord_system_block%set ('nozzle_x', nozzle_x)
    call plist_am_coord_system_block%set ('nozzle_y', nozzle_y)
    call plist_am_coord_system_block%set ('nozzle_z', nozzle_z)
    call plist_am_coord_system_block%set ('nozzle_speed', nozzle_speed)
    call plist_am_coord_system_block%set ('min_x', min_x)
    call plist_am_coord_system_block%set ('max_x', max_x)


  end function plist_am_coord_system_block



end module am_coord_system_data

