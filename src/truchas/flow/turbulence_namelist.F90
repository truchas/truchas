!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module turbulence_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_turbulence_namelist

  type(parameter_list), pointer, public :: params => null()

contains

  subroutine read_turbulence_namelist(lun, plist)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use input_utilities, only: seek_to_namelist, NULL_R
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout), target, optional :: plist

    integer :: ios
    logical :: found
    character(80) :: iom

    real(r8) :: length, cmu, ke_fraction
    namelist /turbulence/ length, cmu, ke_fraction

    if (associated(params)) deallocate(params)

    !! Locate the TURBULENCE namelist (optional)
    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'turbulence', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) return  ! the namelist is optional

    call TLS_info('')
    call TLS_info('Reading TURBULENCE namelist ...')

    !! Read the namelist variables
    if (is_IOP) then
      length = NULL_R
      cmu = NULL_R
      ke_fraction = NULL_R
      read(lun,nml=turbulence,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal ('error reading TURBULENCE namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(length)
    call broadcast(cmu)
    call broadcast(ke_fraction)

    if (present(plist)) then
      params => plist
    else
      allocate(params)
    end if

    call params%set('type', 'alg') ! only support a single type

    if (length == NULL_R) then
      call TLS_fatal('LENGTH must be assigned a value.')
    else if (length <= 0.0_r8) then
      call TLS_fatal('LENGTH must be > 0.0')
    end if
    call params%set('length', length)

    if (cmu /= NULL_R) then
      if (cmu <= 0.0_r8) call TLS_fatal('CMU must be > 0.0')
      call params%set('cmu', cmu)
    end if

    if (ke_fraction /= NULL_R) then
      if (ke_fraction <= 0.0_r8 .or. ke_fraction >= 1.0_r8) &
          call TLS_fatal('KE_FRACTION must be in (0,1)')
      call params%set('ke fraction', ke_fraction)
    end if

    if (present(plist)) params => null()

  end subroutine read_turbulence_namelist

end module turbulence_namelist
