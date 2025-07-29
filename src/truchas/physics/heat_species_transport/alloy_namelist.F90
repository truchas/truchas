module alloy_namelist

  implicit none
  private

  public :: read_alloy_namelist

contains

  subroutine read_alloy_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parameter_list_type
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use string_utilities, only: i_to_c
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: ios
    logical :: found
    character(128) :: iom

    !! Namelist variables
    character(32) :: material
    integer :: num_comp
    real(r8) :: temp_fusion, temp_eutectic, liq_slope(16), part_coef(16), concentration(16), gamma
    namelist /alloy/ material, num_comp, temp_fusion, temp_eutectic, liq_slope, part_coef, concentration, gamma

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'alloy', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call tls_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) return

    call tls_info('Reading ALLOY namelist ...')

    material = NULL_C
    num_comp = NULL_I
    temp_fusion = NULL_R
    temp_eutectic = NULL_R
    liq_slope = NULL_R
    part_coef = NULL_R
    concentration = NULL_R
    gamma = NULL_R

    if (is_IOP) read(lun,nml=alloy,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call tls_fatal('error reading ALLOY namelist: ' // trim(iom))

    call broadcast(material)
    call broadcast(num_comp)
    call broadcast(temp_fusion)
    call broadcast(temp_eutectic)
    call broadcast(liq_slope)
    call broadcast(part_coef)
    call broadcast(concentration)
    call broadcast(gamma)

    if (material == NULL_C) then
      call tls_fatal('MATERIAL not specified')
    else
      call params%set('material', trim(material))
    end if

    if (num_comp == NULL_I) then
      call tls_fatal('NUM_COMP not specified')
    else if (num_comp <= 0) then
      call tls_fatal('NUM_COMP is non-positive')
    else
      call params%set('num-comp', num_comp)
    end if

    if (temp_fusion == NULL_R) then
      call tls_fatal('TEMP_FUSION not specified')
    else
      call params%set('temp-fusion', temp_fusion)
    end if

    if (temp_eutectic /= NULL_R) call params%set('temp-eutectic', temp_eutectic)

    if (all(liq_slope == NULL_R)) then
      call tls_fatal('LIQ_SLOPE not specified')
    else if (any(liq_slope(:num_comp) == NULL_R) .or. any(liq_slope(num_comp+1:) /= NULL_R)) then
      call tls_fatal('invalid number of LIQ_SLOPE values specified')
    else
      call params%set('liq-slope', liq_slope(:num_comp))
    end if

    if (all(part_coef == NULL_R)) then
      call tls_fatal('PART_COEF not specified')
    else if (any(part_coef(:num_comp) == NULL_R) .or. any(part_coef(num_comp+1:) /= NULL_R)) then
      call tls_fatal('invalid number of PART_COEF values specified')
    else
      call params%set('part-coef', part_coef(:num_comp))
    end if

    if (all(concentration == NULL_R)) then
      call tls_fatal('CONCENTRATION not specified')
    else if (any(concentration(:num_comp) == NULL_R) .or. any(concentration(num_comp+1:) /= NULL_R)) then
      call tls_fatal('invalid number of CONCENTRATION values specified')
    else
      call params%set('concentration', concentration(:num_comp))
    end if

    if (gamma /= NULL_R) call params%set('gamma', gamma)

  end subroutine read_alloy_namelist

end module alloy_namelist
