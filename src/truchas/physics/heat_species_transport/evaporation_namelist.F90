module evaporation_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_evaporation_namelist

  type(parameter_list), allocatable, public :: params

contains

  subroutine read_evaporation_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use input_utilities, only: seek_to_namelist, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios
    logical :: found
    character(128) :: iom
    integer, parameter :: MAX_FACE_SET_IDS = 32

    integer  :: face_set_ids(MAX_FACE_SET_IDS)
    real(r8) :: vaporization_heat, vaporization_temp, molar_mass, ambient_pressure, &
        condensation_factor
    namelist /evaporation/ vaporization_heat, vaporization_temp, molar_mass, &
        ambient_pressure, condensation_factor, face_set_ids

    !! Locate the optional EVAPORATION namelist (first occurrence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'EVAPORATION', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) return

    call TLS_info('Reading EVAPORATION namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      vaporization_heat = NULL_R
      vaporization_temp = NULL_R
      molar_mass = NULL_R
      ambient_pressure = NULL_R
      condensation_factor = NULL_R
      face_set_ids = NULL_I
      read(lun,nml=evaporation,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading EVAPORATION namelist: ' // trim(iom))

    !! Broadcast the namelist variables.
    call broadcast(vaporization_heat)
    call broadcast(vaporization_temp)
    call broadcast(molar_mass)
    call broadcast(ambient_pressure)
    call broadcast(condensation_factor)
    call broadcast(face_set_ids)

    !! Verified variables will be stuffed into this parameter list.
    allocate(params)

    if (vaporization_heat == NULL_R) then
      call TLS_fatal('no value assigned to VAPORIZATION_HEAT')
    else if (vaporization_heat < 0.0_r8) then
      call TLS_fatal('VAPORIZATION_HEAT must be >= 0')
    endif
    call params%set('vaporization-heat', vaporization_heat)

    if (vaporization_temp == NULL_R) then
      call TLS_fatal('no value assigned to VAPORIZATION_TEMP')
    else if (vaporization_temp <= 0.0_r8) then
      call TLS_fatal('VAPORIZATION_TEMP must be > 0')
    endif
    call params%set('vaporization-temp', vaporization_temp)

    if (molar_mass == NULL_R) then
      call TLS_fatal('no value assigned to MOLAR_MASS')
    else if (molar_mass <= 0.0_r8) then
      call TLS_fatal('MOLAR_MASS must be > 0')
    endif
    call params%set('molar-mass', molar_mass)

    if (ambient_pressure /= NULL_R) then
      if (ambient_pressure <= 0.0_r8) call TLS_fatal('AMBIENT_PRESSURE must be > 0')
      call params%set('ambient-pressure', ambient_pressure)
    endif

    if (condensation_factor /= NULL_R) then
      if (condensation_factor < 0.0_r8) call TLS_fatal('condensation_factor must be >= 0')
      call params%set('condensation-factor', condensation_factor)
    endif

    !! Check for a non-empty FACE_SET_IDS.
    if (count(face_set_ids /= NULL_I) == 0) call TLS_fatal('no values assigned to FACE_SET_IDS')
    call params%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids /= NULL_I)))

  end subroutine read_evaporation_namelist

end module evaporation_namelist
