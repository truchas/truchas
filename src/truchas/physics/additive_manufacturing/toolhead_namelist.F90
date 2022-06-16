!!
!! TOOLHEAD_NAMELIST
!!
!! Reads the TOOLHEAD namelist whose variables define the action of the
!! toolhead of additive manufacturing processes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module toolhead_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  public :: read_toolhead_namelists

contains

  subroutine read_toolhead_namelists(lun)

    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use toolpath_driver, only: known_toolpath
    use toolhead_driver, only: table
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: n, ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom
    type(parameter_list), pointer :: plist

    character(31) :: name, toolpath, laser_type, laser_power_func
    real(r8) :: laser_direction(3), laser_absorp, laser_time_const, laser_power
    real(r8) :: laser_sigma, laser_wave_length, laser_waist_radius, laser_beam_quality_factor
    namelist /toolhead/ name, toolpath, laser_direction, laser_absorp, laser_time_const, &
        laser_power, laser_power_func, laser_type, laser_sigma, laser_wave_length, & 
        laser_waist_radius, laser_beam_quality_factor

    call TLS_info('Reading TOOLHEAD namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do  ! until all TOOLHEAD namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'toolhead', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'TOOLHEAD[' // i_to_c(n) // ']'

      name = NULL_C
      toolpath = NULL_C
      laser_direction = NULL_R
      laser_absorp = NULL_R
      laser_time_const = NULL_R
      laser_power = NULL_R
      laser_power_func = NULL_C
      laser_type = NULL_C
      laser_sigma = NULL_R
      laser_wave_length = NULL_R
      laser_waist_radius = NULL_R
      laser_beam_quality_factor = NULL_R

      if (is_IOP) read(lun,nml=toolhead,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(toolpath)
      call broadcast(laser_direction)
      call broadcast(laser_absorp)
      call broadcast(laser_time_const)
      call broadcast(laser_power)
      call broadcast(laser_power_func)
      call broadcast(laser_type)
      call broadcast(laser_sigma)
      call broadcast(laser_wave_length)
      call broadcast(laser_waist_radius)
      call broadcast(laser_beam_quality_factor)

      !! A unique NAME is required; becomes the toolheads sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (table%known_toolhead(name)) then
        call TLS_fatal(label // ': another TOOLHEAD namelist has this NAME: ' // trim(name))
      end if

      !! Check the variables (FIXME: need more checks)
      if (toolpath == NULL_C) call TLS_fatal(label // ': TOOLPATH not specified')
      if (.not.known_toolpath(toolpath)) call TLS_fatal(label // ': unknown TOOLPATH name: ' // trim(toolpath))
      if (all(laser_direction == NULL_R)) call TLS_fatal(label // ': LASER_DIRECTION not defined')
      if (any(laser_direction == NULL_R)) call TLS_fatal(label // ': LASER_DIRECTION not fully defined')
      if (laser_absorp == NULL_R) call TLS_fatal(label // ': LASER_ABSORP not defined')
      if (laser_time_const == NULL_R) call TLS_fatal(label // ': LASER_TIME_CONST not defined')
      if (laser_power == NULL_R .and. laser_power_func == NULL_C) call TLS_fatal('LASER_POWER not defined')
      if (laser_type == NULL_C) call TLS_fatal(label // ': LASER_TYPE not defined')
      select case (laser_type)
      case ('gaussian')
        if (laser_sigma == NULL_R) call TLS_fatal(label // ': LASER_SIGMA not defined')
      case ('gaussian beam')
        if (laser_wave_length == NULL_R) call TLS_fatal(label // ': LASER_WAVE_LENGTH not defined')
        if (laser_beam_quality_factor == NULL_R) call TLS_fatal(label // ': LASER_BEAM_QUALITY_FACTOR not defined')
        if (laser_waist_radius == NULL_R) call TLS_fatal(label // ': LASER_WAIST_RADIUS not defined')
      case default
        call TLS_fatal(label // ': unknown LASER_TYPE: ' // trim(laser_type))
      end select

      call TLS_info('  read namelist "' // trim(name) // '"')

      call make_toolhead
    end do

    if (n == 0) call TLS_info('  none found')

  contains

    subroutine make_toolhead

      use parameter_list_type
      use toolhead_type

      type(parameter_list), pointer :: params, plist
      type(toolhead), allocatable :: head

      !! Stuff them into a parameter list for later use
      allocate(params)
      call params%set('toolpath', toolpath)
      call params%set('laser-direction', laser_direction)
      call params%set('laser-absorp', laser_absorp)
      call params%set('laser-time-constant', laser_time_const)
      plist => params%sublist('laser')
      call process2(plist, laser_power, laser_power_func, 'LASER_POWER', 'power', label)
      call plist%set('type', laser_type)
      select case (laser_type)
      case ('gaussian')
        call plist%set('sigma', laser_sigma)
      case ('gaussian beam')
        call plist%set('wave-length', laser_wave_length)
        call plist%set('waist-radius', laser_waist_radius)
        call plist%set('beam-quality-factor', laser_beam_quality_factor)
      end select

      allocate(head)
      call head%init(params)
      call table%add_toolhead(name, head)

    end subroutine

  end subroutine read_toolhead_namelists

  subroutine process2(plist, const, fname, bname, pname, label)
    use input_utilities, only: NULL_C, NULL_R
    use scalar_func_table
    type(parameter_list), intent(inout) :: plist
    real(r8),     intent(in) :: const ! possible constant value
    character(*), intent(in) :: fname ! possible function name value
    character(*), intent(in) :: bname ! namelist variable base name
    character(*), intent(in) :: pname ! parameter list variable name
    character(*), intent(in) :: label
    if (const /= NULL_R .and. fname /= NULL_C) then
      call TLS_fatal(label // ': cannot specify both' // bname // ' and ' // bname // '_FUNC')
    else if (const /= NULL_R) then
      call plist%set(pname, const)
    else if (fname /= NULL_C) then
      if (known_func(fname)) then
        call plist%set(pname, trim(fname))
      else
        call TLS_fatal(label // ': unknown function for ' // bname // '_FUNC: ' // trim(fname))
      end if
    end if
  end subroutine


end module toolhead_namelist
