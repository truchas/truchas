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
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I, NULL_C
    use string_utilities, only: i_to_c
    use material_model_driver, only: matl_model
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios, n, j
    logical :: found
    character(80) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer  :: vol_track_subcycles, fischer_dim
    logical  :: inviscid, track_interfaces, nested_dissection, void_collapse
    logical  :: wisp_redistribution
    real(r8) :: viscous_implicitness, viscous_number, courant_number
    real(r8) :: fluid_frac_threshold, min_face_fraction, vol_frac_cutoff
    real(r8) :: void_collapse_relaxation, wisp_cutoff, wisp_absorption_fraction
    real(r8) :: wisp_neighbor_cutoff
    character(32) :: material_priority(16)
    namelist /flow/ inviscid, &
        viscous_implicitness, viscous_number, courant_number, &
        fluid_frac_threshold, min_face_fraction, &
        track_interfaces, nested_dissection, vol_track_subcycles, &
        vol_frac_cutoff, material_priority, fischer_dim, void_collapse, &
        void_collapse_relaxation, wisp_redistribution, wisp_cutoff, &
        wisp_absorption_fraction, wisp_neighbor_cutoff

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
    track_interfaces = (matl_model%nphase > 1)
    nested_dissection = .true.
    vol_track_subcycles = NULL_I
    vol_frac_cutoff = NULL_R
    material_priority = NULL_C
    wisp_redistribution = .false.
    wisp_cutoff = NULL_R
    wisp_absorption_fraction = NULL_R
    wisp_neighbor_cutoff = NULL_R

    inviscid = .false.
    viscous_implicitness = NULL_R
    viscous_number = NULL_R
    courant_number = NULL_R

    fluid_frac_threshold = NULL_R
    min_face_fraction = NULL_R
    fischer_dim = NULL_I

    void_collapse = .false.
    void_collapse_relaxation = NULL_R

    !! Read the FLOW namelist
    if (is_IOP) read(lun,nml=flow,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading FLOW namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(track_interfaces)
    call broadcast(nested_dissection)
    call broadcast(vol_track_subcycles)
    call broadcast(vol_frac_cutoff)
    call broadcast(material_priority)

    call broadcast(inviscid)
    call broadcast(viscous_implicitness)
    call broadcast(viscous_number)
    call broadcast(courant_number)

    call broadcast(fluid_frac_threshold)
    call broadcast(min_face_fraction)
    call broadcast(fischer_dim)

    call broadcast(void_collapse)
    call broadcast(void_collapse_relaxation)

    call broadcast(wisp_redistribution)
    call broadcast(wisp_cutoff)
    call broadcast(wisp_absorption_fraction)
    call Broadcast(wisp_neighbor_cutoff)
    !! Check values and stuff into a parameter list for later use.

    plist => params%sublist('volume-tracker')
    call plist%set('track_interfaces', track_interfaces)
    call plist%set('nested_dissection', nested_dissection)
    call plist%set('wisp_redistribution', wisp_redistribution)

    n = count(material_priority /= NULL_C)
    if (n > 0) then
      do j = 1, n
        if (material_priority(j) == 'SOLID') cycle
        if (matl_model%phase_index(material_priority(j)) == 0) &
            call TLS_fatal('unknown material for MATERIAL_PRIORITY: ' // trim(material_priority(j)))
      end do
      call plist%set('material_priority', material_priority(:n))
    end if

    if (vol_track_subcycles /= NULL_I) then
      if (vol_track_subcycles < 1) call TLS_fatal('SUBCYCLES must be > 0')
      call plist%set('subcycles', vol_track_subcycles)
    end if

    if (vol_frac_cutoff /= NULL_R) then
      if (vol_frac_cutoff <= 0.0_r8) call TLS_fatal('VOL_FRAC_CUTOFF must be > 0.0')
      call plist%set('cutoff', vol_frac_cutoff)
    end if

    if (wisp_cutoff /= NULL_R) then
      if (wisp_cutoff < 0.0_r8 .or. wisp_cutoff > 1.0_r8) &
          call TLS_fatal("WISP_CUTOFF must be in [0, 1]")
      call plist%set('wisp_cutoff', wisp_cutoff)
    end if

    if (wisp_absorption_fraction /= NULL_R) then
      if (wisp_absorption_fraction <= 0.0_r8 .or. wisp_absorption_fraction > 1.0_r8) &
          call TLS_fatal("WISP_ABSORBTION_FRACTION must be in (0, 1]")
      call plist%set('wisp_absorbtion_fraction', wisp_absorption_fraction)
    end if

    if (wisp_neighbor_cutoff /= NULL_R) then
      if (wisp_neighbor_cutoff <= 0.0_r8 .or. wisp_neighbor_cutoff > 1.0_r8) &
          call TLS_fatal("WISP_NEIGHBOR_CUTOFF must be in (0, 1]")
      call plist%set('wisp_neighbor_cutoff', wisp_neighbor_cutoff)
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

    call plist%set('void-collapse', void_collapse)
    if (void_collapse_relaxation /= NULL_R) then
      if (void_collapse_relaxation < 0 .or. void_collapse_relaxation > 1) &
          call TLS_fatal('VOID COLLAPSE RELAXATION must be in [0,1]')
      call plist%set('void-collapse-relaxation', void_collapse_relaxation)
    end if

    plist => params%sublist('cutoffs')
    if (fluid_frac_threshold /= NULL_R) then
      if (fluid_frac_threshold <= 0.0_r8) call TLS_fatal('FLUID_FRAC_THRESHOLD must be > 0.0')
      call plist%set('cutvof', fluid_frac_threshold)
    end if

    if (min_face_fraction /= NULL_R) then
      if (min_face_fraction <= 0.0_r8) call TLS_fatal('MIN_FACE_FRACTION must be > 0.0')
      call plist%set('min face fraction', min_face_fraction)
    end if

  end subroutine read_flow_namelist

end module flow_namelist
