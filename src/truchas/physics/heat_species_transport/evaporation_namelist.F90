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
    real(r8) :: prefactor, temp_exponent, activation_energy
    namelist /evaporation/ prefactor, temp_exponent, activation_energy, face_set_ids

    !! Locate the optional EVAPORATION namelist (first occurrence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'EVAPORATION', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) return

    call TLS_info('')
    call TLS_info('Reading EVAPORATION namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      prefactor = NULL_R
      temp_exponent = NULL_R
      activation_energy = NULL_R
      face_set_ids = NULL_I
      read(lun,nml=evaporation,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading EVAPORATION namelist: ' // trim(iom))

    !! Broadcast the namelist variables.
    call broadcast(prefactor)
    call broadcast(temp_exponent)
    call broadcast(activation_energy)
    call broadcast(face_set_ids)
    
    !! Verified variables will be stuffed into this parameter list.
    allocate(params)

    !! Postive PREFACTOR value required.
    if (prefactor == NULL_R) then
      call TLS_fatal('no value assigned to PREFACTOR')
    else if (prefactor <= 0.0_r8) then
      call TLS_fatal('PREFACTOR must be > 0')
    endif
    call params%set('prefactor', prefactor)
    
    !! TEMP_EXPONENT is optional; defaults to 0.
    if (temp_exponent /= NULL_R) call params%set('temp-exponent', temp_exponent)
    
    !! Positive ACTIVATION_ENERGY value required.
    if (activation_energy == NULL_R) then
      call TLS_fatal('no value assigned to ACTIVATION_ENERGY')
    else if (activation_energy <= 0.0_r8) then
      call TLS_fatal('ACTIVATION_ENERGY must be > 0')
    endif
    call params%set('activation-energy', activation_energy)
    
    !! Check for a non-empty FACE_SET_IDS.
    if (count(face_set_ids /= NULL_I) == 0) call TLS_fatal('no values assigned to FACE_SET_IDS')
    call params%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids /= NULL_I)))

  end subroutine read_evaporation_namelist

end module evaporation_namelist
