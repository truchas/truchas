!!
!! DED_HEAD_NAMELIST
!!
!! Reads the DED_HEAD namelist whose variables define the action of the
!! tool head of a LENS-like directed energy deposition (DED) AM process.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ded_head_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  implicit none
  private

  public :: read_ded_head_namelist

  type(parameter_list), pointer, public :: ded_params => null()

contains

  subroutine read_ded_head_namelist (lun)

    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use toolpath_table, only: known_toolpath
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios
    logical :: found
    type(parameter_list), pointer :: plist

    character(31) :: toolpath, laser_type
    real(r8) :: laser_absorp, laser_time_const, laser_power
    real(r8) :: laser_sigma, laser_wave_length, laser_waist_radius, laser_beam_param
    namelist /ded_head/ toolpath, laser_absorp, laser_time_const, laser_power, laser_type, &
        laser_sigma, laser_wave_length, laser_waist_radius, laser_beam_param

    !! Locate the optional DED_HEAD namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'DED_HEAD', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) return

    call TLS_info('')
    call TLS_info('Reading DED_HEAD namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      toolpath = NULL_C
      laser_absorp = NULL_R
      laser_time_const = NULL_R
      laser_power = NULL_R
      laser_type = NULL_C
      laser_sigma = NULL_R
      laser_wave_length = NULL_R
      laser_waist_radius = NULL_R
      laser_beam_param = NULL_R
      read(lun,nml=ded_head,iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading DED_HEAD namelist')

    !! Broadcast the namelist variables.
    call broadcast(toolpath)
    call broadcast(laser_absorp)
    call broadcast(laser_time_const)
    call broadcast(laser_power)
    call broadcast(laser_type)
    call broadcast(laser_sigma)
    call broadcast(laser_wave_length)
    call broadcast(laser_waist_radius)
    call broadcast(laser_beam_param)

    !! Check the variables (FIXME: need more checks)
    if (toolpath == NULL_C) call TLS_fatal('TOOLPATH not specified')
    if (.not.known_toolpath(toolpath)) call TLS_fatal('unknown TOOLPATH name: ' // trim(toolpath))
    if (laser_absorp == NULL_R) call TLS_fatal('LASER_ABSORP not defined')
    if (laser_time_const == NULL_R) call TLS_fatal('LASER_TIME_CONST not defined')
    if (laser_power == NULL_R) call TLS_fatal('LASER_POWER not defined')
    if (laser_type == NULL_C) call TLS_fatal('LASER_TYPE not defined')
    select case (laser_type)
    case ('gaussian')
      if (laser_sigma == NULL_R) call TLS_fatal('LASER_SIGMA not defined')
    case ('gaussian beam')
      if (laser_wave_length == NULL_R) call TLS_fatal('LASER_WAVE_LENGTH not defined')
      if (laser_beam_param == NULL_R) call TLS_fatal('LASER_BEAM_PARAM not defined')
      if (laser_waist_radius == NULL_R) call TLS_fatal('LASER_WAIST_RADIUS not defined')
    case default
      call TLS_fatal('unknown LASER_TYPE: ' // trim(laser_type))
    end select

    !! Stuff them into a parameter list for later use
    allocate(ded_params)
    call ded_params%set('toolpath', toolpath)
    call ded_params%set('laser-absorp', laser_absorp)
    call ded_params%set('laser-time-constant', laser_time_const)
    plist => ded_params%sublist('laser')
    call plist%set('power', laser_power)
    call plist%set('type', laser_type)
    select case (laser_type)
    case ('gaussian')
      call plist%set('sigma', laser_sigma)
    case ('gaussian beam')
      call plist%set('wave-length', laser_wave_length)
      call plist%set('waist-radius', laser_waist_radius)
      call plist%set('beam-param', laser_beam_param)
    end select

  end subroutine read_ded_head_namelist

end module ded_head_namelist
