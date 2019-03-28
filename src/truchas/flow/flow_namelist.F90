!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flow_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_flow_namelist

  type(parameter_list), public :: params

contains

  subroutine read_flow_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I
    use string_utilities, only: i_to_c
    use parameter_module, only: nmat
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios, n
    logical :: found
    character(80) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer  :: subcycles, material_priority(16), fischer_dim
    logical  :: inviscid, track_interfaces, nested_dissection
    real(r8) :: viscous_implicitness, viscous_number, courant_number
    real(r8) :: fluid_cutvof, min_face_fraction, cutoff
    namelist /flow/ inviscid, &
        viscous_implicitness, viscous_number, courant_number, &
        fluid_cutvof, min_face_fraction, &
        track_interfaces, nested_dissection, subcycles, &
        cutoff, material_priority, fischer_dim

    call TLS_info('')
    call TLS_info('Reading FLOW namelist ...')

    !! Locate the FLOW namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'flow', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('FLOW namelist not found')

    !! Default values
    track_interfaces = (nmat > 1)
    nested_dissection = .true.
    subcycles = NULL_I
    cutoff = NULL_R
    material_priority = NULL_I

    inviscid = .false.
    viscous_implicitness = NULL_R
    viscous_number = NULL_R
    courant_number = NULL_R

    fluid_cutvof = NULL_R
    min_face_fraction = NULL_R
    fischer_dim = NULL_I

    !! Read the FLOW namelist
    if (is_IOP) read(lun,nml=flow,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading FLOW namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(track_interfaces)
    call broadcast(nested_dissection)
    call broadcast(subcycles)
    call broadcast(cutoff)
    call broadcast(material_priority)

    call broadcast(inviscid)
    call broadcast(viscous_implicitness)
    call broadcast(viscous_number)
    call broadcast(courant_number)

    call broadcast(fluid_cutvof)
    call broadcast(min_face_fraction)
    call broadcast(fischer_dim)

    !! Check values and stuff into a parameter list for later use.

    plist => params%sublist('volume-tracker')
    call plist%set('track_interfaces', track_interfaces)
    call plist%set('nested_dissection', nested_dissection)

    n = count(material_priority /= NULL_I)
    if (n > 0) then
      call plist%set('material_priority', material_priority(:n))
    end if

    if (subcycles /= NULL_I) then
      if (subcycles < 1) call TLS_fatal('SUBCYCLES must be > 0')
      call plist%set('subcycles', subcycles)
    end if

    if (cutoff /= NULL_R) then
      if (cutoff <= 0.0_r8) call TLS_fatal('CUTOFF must be > 0.0')
      call plist%set('cutoff', cutoff)
    end if

    plist => params%sublist('options')
    call plist%set('inviscid', inviscid)

    if (viscous_implicitness /= NULL_R) then
      if (viscous_implicitness < 0 .or. viscous_implicitness > 1) &
          call TLS_fatal('VISCOUS_IMPLICITNESS must be in [0,1]')
      call plist%set('viscous-implicitness', viscous_implicitness)
    end if

    if (viscous_number /= NULL_R) then
      if (viscous_number < 0.0_r8) call TLS_fatal('VISCOUS_NUMBER must be >= 0.0')
      call plist%set('viscous-number', viscous_number)
    end if

    if (courant_number /= NULL_R) then
      if (courant_number <= 0 .or. courant_number > 1) &
          call TLS_fatal('COURANT_NUMBER must be in (0,1]')
      call plist%set('courant number', courant_number)
    end if

    if (fischer_dim /= NULL_I) then
      if (fischer_dim < 0) call TLS_fatal('FISCHER_DIM must be >= 0')
      call plist%set('fischer-dim', fischer_dim)
    end if

    plist => params%sublist('cutoffs')
    if (fluid_cutvof /= NULL_R) then
      if (fluid_cutvof <= 0.0_r8) call TLS_fatal('FLUID_CUTVOF must be > 0.0')
      call plist%set('cutvof', fluid_cutvof)
    end if

    if (min_face_fraction /= NULL_R) then
      if (min_face_fraction <= 0.0_r8) call TLS_fatal('MIN_FACE_FRACTION must be > 0.0')
      call plist%set('min face fraction', min_face_fraction)
    end if

  end subroutine read_flow_namelist

end module flow_namelist
