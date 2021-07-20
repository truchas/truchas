!!
!! THERMAL_BC_NAMELIST
!!
!! Provides a subroutine for reading the THERMAL_BC namelists for a file and
!! storing the data in a parameter list for later use.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module thermal_bc_namelist

  implicit none
  private

  public :: read_thermal_bc_namelists

contains

  subroutine read_thermal_bc_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use string_utilities, only: i_to_c, lower_case
    use parameter_list_type
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: face_set_ids(100)
    real(r8) :: htc, ambient_temp, emissivity, temp, flux, vflux(3)
    character(32) :: name, type, htc_func, ambient_temp_func, &
        emissivity_func, temp_func, flux_func, vflux_func
    namelist /thermal_bc/ name, type, face_set_ids, &
        temp, temp_func, flux, flux_func, htc, htc_func, ambient_temp, ambient_temp_func, &
        emissivity, emissivity_func, vflux, vflux_func

    call TLS_info('Reading THERMAL_BC namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all THERMAL_BC namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'thermal_bc', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'THERMAL_BC[' // i_to_c(n) // ']'

      name = NULL_C
      type = NULL_C
      face_set_ids = NULL_I
      htc = NULL_R
      htc_func = NULL_C
      ambient_temp = NULL_R
      ambient_temp_func = NULL_C
      emissivity = NULL_R
      emissivity_func = NULL_C
      temp = NULL_R
      temp_func = NULL_C
      flux = NULL_R
      flux_func = NULL_C
      vflux = NULL_R
      vflux_func = NULL_C

      if (is_IOP) read(lun,nml=thermal_bc,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(type)
      call broadcast(face_set_ids)
      call broadcast(htc)
      call broadcast(htc_func)
      call broadcast(ambient_temp)
      call broadcast(ambient_temp_func)
      call broadcast(emissivity)
      call broadcast(emissivity_func)
      call broadcast(temp)
      call broadcast(temp_func)
      call broadcast(flux)
      call broadcast(flux_func)
      call broadcast(vflux)
      call broadcast(vflux_func)

      !! A unique NAME is required; becomes the BC sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another THERMAL_BC namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      !! FACE_SET_IDS is required; cannot check that they are valid at this point.
      if (count(face_set_ids /= NULL_I) == 0) then
        call TLS_fatal(label // ': FACE_SET_IDS not specified')
      else
        call plist%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids/=NULL_I)))
      end if

      !! Check the required TYPE value.
      select case (lower_case(type))
      case ('temperature', 'flux', 'oriented-flux', 'htc', 'radiation')
      case ('interface-htc', 'gap-radiation')
      case (NULL_C)
        call TLS_fatal(label // ': TYPE not specified')
      case default
        call TLS_fatal(label // ': unknown TYPE: ' // trim(type))
      end select
      call plist%set('type', trim(type))

      select case (lower_case(type))
      case ('temperature')

        if (temp /= NULL_R .and. temp_func /= NULL_C) then
          call TLS_fatal(label // ': both TEMP and TEMP_FUNC specified')
        else if (temp /= NULL_R) then
          call plist%set('temp', temp)
        else if (temp_func /= NULL_C) then
          call plist%set('temp', temp_func)
        else
          call TLS_fatal(label // ': neither TEMP or TEMP_FUNC specified')
        end if

      case ('flux')

        if (flux /= NULL_R .and. flux_func /= NULL_C) then
          call TLS_fatal(label // ': both FLUX and FLUX_FUNC specified')
        else if (flux /= NULL_R) then
          call plist%set('flux', flux)
        else if (flux_func /= NULL_C) then
          call plist%set('flux', flux_func)
        else
          call TLS_fatal(label // ': neither FLUX or FLUX_FUNC specified')
        end if

      case ('oriented-flux')

        if (any(vflux /= NULL_R) .and. vflux_func /= NULL_C) then
          call TLS_fatal(label // ': both VFLUX and VFLUX_FUNC specified')
        else if (any(vflux /= NULL_R)) then
          if (any(vflux == NULL_R)) then
            call TLS_fatal(label // ': VFLUX only partially specified')
          else
            call plist%set('flux', vflux)
          end if
        else if (vflux_func /= NULL_C) then
          call plist%set('flux', vflux_func)
        else
          call TLS_fatal(label // ': neither VFLUX or VFLUX_FUNC specified')
        end if

      case ('htc')

        if (htc /= NULL_R .and. htc_func /= NULL_C) then
          call TLS_fatal(label // ': both HTC and HTC_FUNC specified')
        else if (htc /= NULL_R) then
          if (htc < 0) call TLS_fatal(label // ': HTC < 0.0')
          call plist%set('htc', htc)
        else if (htc_func /= NULL_C) then
          call plist%set('htc', htc_func)
        else
          call TLS_fatal(label // ': neither HTC or HTC_FUNC specified')
        end if

        if (ambient_temp /= NULL_R .and. ambient_temp_func /= NULL_C) then
          call TLS_fatal(label // ': both AMBIENT_TEMP and AMBIENT_TEMP_FUNC specified')
        else if (ambient_temp /= NULL_R) then
          call plist%set('ambient-temp', ambient_temp)
        else if (ambient_temp_func /= NULL_C) then
          call plist%set('ambient-temp', ambient_temp_func)
        else
          call TLS_fatal(label // ': neither AMBIENT_TEMP or AMBIENT_TEMP_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      case ('radiation')

        if (emissivity /= NULL_R .and. emissivity_func /= NULL_C) then
          call TLS_fatal(label // ': both EMISSIVITY and EMISSIVITY_FUNC specified')
        else if (emissivity /= NULL_R) then
          if (emissivity < 0) call TLS_fatal(label // ': EMISSIVITY < 0.0')
          if (emissivity > 1) call TLS_fatal(label // ': EMISSIVITY > 1.0')
          call plist%set('emissivity', emissivity)
        else if (emissivity_func /= NULL_C) then
          call plist%set('emissivity', emissivity_func)
        else
          call TLS_fatal(label // ': neither EMISSIVITY OR EMISSIVITY_FUNC specified')
        end if

        if (ambient_temp /= NULL_R .and. ambient_temp_func /= NULL_C) then
          call TLS_fatal(label // ': both AMBIENT_TEMP and AMBIENT_TEMP_FUNC specified')
        else if (ambient_temp /= NULL_R) then
          call plist%set('ambient-temp', ambient_temp)
        else if (ambient_temp_func /= NULL_C) then
          call plist%set('ambient-temp', ambient_temp_func)
        else
          call TLS_fatal(label // ': neither AMBIENT_TEMP OR AMBIENT_TEMP_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      case ('interface-htc')

        if (htc /= NULL_R .and. htc_func /= NULL_C) then
          call TLS_fatal(label // ': both HTC and HTC_FUNC specified')
        else if (htc /= NULL_R) then
          if (htc < 0) call TLS_fatal(label // ': htc < 0.0')
          call plist%set('htc', htc)
        else if (htc_func /= NULL_C) then
          call plist%set('htc', htc_func)
        else
          call TLS_fatal(label // ': neither HTC or HTC_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      case ('gap-radiation')

        if (emissivity /= NULL_R .and. emissivity_func /= NULL_C) then
          call TLS_fatal(label // ': both EMISSIVITY and EMISSIVITY_FUNC specified')
        else if (emissivity /= NULL_R) then
          if (emissivity < 0) call TLS_fatal(label // ': EMISSIVITY < 0.0')
          if (emissivity > 1) call TLS_fatal(label // ': EMISSIVITY > 1.0')
          call plist%set('emissivity', emissivity)
        else if (emissivity_func /= NULL_C) then
          call plist%set('emissivity', emissivity_func)
        else
          call TLS_fatal(label // ': neither EMISSIVITY OR EMISSIVITY_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      end select

      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    if (n == 0) call TLS_info('  none found')

  end subroutine read_thermal_bc_namelists

end module thermal_bc_namelist
