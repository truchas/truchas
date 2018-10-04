!!
!! FLOW_DRIVER
!!
!! This code is meant to be deleted and subsumed by high-level cycle driver.
!! At the moment, this driver serves as a mediator between the flow
!! code with the desired api and the rest of truchas
!!
!! Peter Brady <ptb@lanl.gov>
!! 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMMING INTERFACE
!!
!!  This module provides procedures for driving the microstructure modeling
!!  kernel.  It serves as an adapter between top-level Truchas drivers (input,
!!  initialization, cycle, and output) and the microstructure modeling object.
!!  The existing top-level drivers manage no data; they merely provide for
!!  orchestrating the sequence of procedural steps.  Thus a major role of
!!  these subroutines is to assemble the data from various Truchas components
!!  that must be passed the the microstructure modeling object.
!!
!!  CALL READ_FLOW_NAMELIST (LUN) reads several different flow
!!    related namelists:
!!      flow_driver -> must be specified to use new flow implimentation
!!      flow_cutoffs -> various magic numbers to smooth over numerical difficulties
!!      flow_options -> options for flow that don't really fit into
!!      the the predictor/corrector split.
!!      flow_turbulence_model
!!      flow_predictor
!!      flow_corrector
!!      flow_bc
!!    After stuffing all the namelist inputs into a parameter list, read_flow_params
!!    is then called with that parameter list.  Thus if you have a parameter list
!!    you may directly call read_flow_params
!!
!!
!!  CALL FLOW_DRIVER_INIT (T) initializes the driver.  T is the initial time.
!!    This should be called only if heat transfer physics is enabled, and after
!!    its initialization.  If microstructure analysis has not been enabled this
!!    subroutine does nothing, and so may always be called.
!!
!!  CALL VTRACK_UPDATE (T) updates or advances the microstructure model to the
!!    next time level.  The subroutine handles collecting all the necessary
!!    mesh-based state arrays required by the model; only the current time T
!!    needs to be passed.
!!
!! NOTES
!!
!! The allocatation status of the private module variable THIS serves
!! to indicate whether the new flow is enabled.  It is allocated by
!! READ_FLOW_NAMELIST if the flow namelist is present in the input
!! file.
!!
#include "f90_assert.fpp"
#define ASDF 0

module flow_driver

  use kinds, only: r8
  use unstr_mesh_type
  use flow_mesh_type
  use flow_type
  use flow_props_type
  use parameter_list_type
  use phase_property_table
  use truchas_logging_services
  use truchas_timers
  use index_partitioning
  use scalar_func_class
  use scalar_func_tools
  use scalar_func_containers
  use flow_operators, only: flow_operators_init
  use flow_bc_type, only: read_flow_bc_namelist
  use turbulence_models, only: read_flow_turbulence_model_namelist
  use flow_prediction_type, only: read_flow_predictor_namelist
  use flow_projection_type, only: read_flow_corrector_namelist
  implicit none
  private

  public :: read_flow_namelist, flow_timestep
  public :: flow_driver_init, flow_step, flow_final, flow_enabled, flow_accept
  public :: flow_vel_fn_view, flow_vel_cc_view, flow_P_cc_view

  type :: flow_driver_data
    type(flow_mesh), pointer :: mesh => null() ! OWNING
    type(flow) :: flow
    type(flow_props) :: props
    ! The flow driver shouldn't logically need this but temperature is currently
    ! stored in Zone%Temp so we need to keep a version on the new mesh as well
    real(r8), allocatable :: temperature_cc(:)
    logical :: active
  end type flow_driver_data

  type(flow_driver_data), allocatable, target :: this

contains

  function flow_vel_fn_view() result(p)
    real(r8), pointer :: p(:)
    ASSERT(flow_enabled())
    p => this%flow%vel_fn_view()
  end function flow_vel_fn_view

  function flow_vel_cc_view() result(p)
    real(r8), pointer :: p(:,:)
    ASSERT(flow_enabled())
    p => this%flow%vel_cc_view()
  end function flow_vel_cc_view

  function flow_P_cc_view() result(p)
    real(r8), pointer :: p(:)
    ASSERT(flow_enabled())
    p => this%flow%P_cc_view()
  end function flow_P_cc_view

  subroutine flow_final
    if (allocated(this)) deallocate(this)
  end subroutine flow_final

  subroutine flow_timestep(dtc, dtv)
    real(r8), intent(out) :: dtc, dtv
    ASSERT(flow_enabled())
    call this%flow%timestep(this%props, dtc, dtv)
  end subroutine flow_timestep

  logical function flow_enabled()
    flow_enabled = allocated(this)
  end function flow_enabled

  subroutine read_flow_params(p)
    type(parameter_list), pointer, intent(in) :: p
    call this%props%read_params(p)
    call this%flow%read_params(p)
  end subroutine read_flow_params

  subroutine read_flow_options_namelist(lun, p)
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use flow_input_utils
    use physics_module, only: body_force_density

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(inout) :: p
    type(parameter_list), pointer :: pp
    integer :: ios
    logical :: found
    character(128) :: iom
    ! flow_options namelist
    logical :: inviscid, stokes
    real(r8) :: viscous_implicitness, solidify_implicitness
    real(r8) :: viscous_number, courant_number

    namelist /flow_options/ inviscid, stokes, &
        viscous_implicitness, solidify_implicitness, viscous_number, &
        courant_number

    pp => p%sublist("options")
    inviscid = .false.
    stokes = .false.
    viscous_implicitness = null_r
    solidify_implicitness = null_r
    viscous_number = null_r
    courant_number = null_r

    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'FLOW_OPTIONS', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('')
      call TLS_info('Reading FLOW_OPTIONS namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        read(lun,nml=flow_options,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOW_OPTIONS namelist: ' // trim(iom))

      if (inviscid .and. stokes) &
          call TLS_fatal("inviscid Stokes flow not supported")

      call broadcast(inviscid)
      call broadcast(stokes)
      call broadcast(viscous_implicitness)
      call broadcast(solidify_implicitness)
      call broadcast(viscous_number)
      call broadcast(courant_number)

      ! active is an option for the driver so don't put it in a plist
      call plist_set_if(pp, 'viscous implicitness', viscous_implicitness)
      call plist_set_if(pp, 'solidify implicitness', solidify_implicitness)
      call plist_set_if(pp, 'courant number', courant_number)
      call plist_set_if(pp, 'viscous number', viscous_number)
      call plist_set_if(pp, 'body force', body_force_density)
    end if

    call pp%set('inviscid', inviscid)
    call pp%set('stokes', stokes)

  end subroutine read_flow_options_namelist

  subroutine read_flow_namelist(lun)
    use string_utilities, only: i_to_c
    use flow_input_utils

    integer, intent(in) :: lun
    integer :: ios
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: p

    if (.not.allocated(this)) allocate(this)

    allocate(p)

    call read_flow_options_namelist(lun, p)

    ! This needs a better name... right now flow_props keeps track of
    ! all the cutoff values since they are used to compute rho/mu
    call read_flow_props_namelist(lun, p)
    call read_flow_bc_namelist(lun, p)
    call read_flow_turbulence_model_namelist(lun, p)
    call read_flow_predictor_namelist(lun, p)
    call read_flow_corrector_namelist(lun, p)

    call read_flow_params(p)
  end subroutine read_flow_namelist


  subroutine flow_driver_init(mesh)
    use zone_module, only: zone
    use material_interop, only: void_material_index, material_to_phase
    use physics_module, only: prescribed_flow
    use scalar_func_factories, only: alloc_const_scalar_func
    use property_data_module, only: isImmobile
    use parameter_module, only: nmat

    type(unstr_mesh), pointer, intent(in) :: mesh
    !real(r8), intent(in) :: vof(:,:), flux_vol(:,:)
    !-
    integer :: i, property_id
    integer, allocatable :: fluids(:)
    logical, allocatable :: is_real_fluid(:)
    real(r8), allocatable :: density(:), velocity_cc(:,:)
    real(r8) :: state(1)
    class(scalar_func_box), allocatable :: density_delta(:), viscosity(:)
    class(scalar_func), allocatable :: f

    if (.not.allocated(this)) return

    allocate(this%mesh)
    call this%mesh%init(mesh)

    allocate(this%temperature_cc(mesh%ncell))

    ! Some duplication here from vtrack_driver.  This should all be subsumed and handled
    ! properly by a sufficiently intelligent physics driver at some point
    if (all(isImmobile(:nmat))) then ! no fluids
      deallocate(this)
      return
    end if

    ! identify non-void fluids
    is_real_fluid = .not.isImmobile(:nmat)
    if (void_material_index > 0) is_real_fluid(void_material_index) = .false.
    fluids = pack(material_to_phase, is_real_fluid)

    ! get density
    state(1) = 0.0_r8
    property_id = ppt_property_id("density")
    allocate(density(size(fluids)))
    do i = 1, size(fluids)
      call ppt_get_phase_property(fluids(i), property_id, f)
      ASSERT(allocated(f))
      ASSERT(is_const(f))
      density(i) = f%eval(state)
    end do

    ! get viscosity
    ! if not given, assume 0-valued viscosity
    allocate(viscosity(size(fluids)))
    if (ppt_has_property("viscosity")) then
      property_id = ppt_property_id("viscosity")
      do i = 1, size(fluids)
        if (ppt_has_phase_property(fluids(i), property_id)) then
          call ppt_get_phase_property(fluids(i), property_id, f)
          call move_alloc(f, viscosity(i)%f)
        else
          call alloc_const_scalar_func(viscosity(i)%f, 0.0_r8)
        end if
      end do
    else
      do i = 1, size(fluids)
        call alloc_const_scalar_func(viscosity(i)%f, 0.0_r8)
      end do
    end if

    ! get density deviations
    ! if not given, assume 0-valued density deviation
    allocate(density_delta(size(fluids)))
    if (ppt_has_property("density deviation")) then
      property_id = ppt_property_id("density deviation")
      do i = 1, size(fluids)
        if (ppt_has_phase_property(fluids(i), property_id)) then
          call ppt_get_phase_property(fluids(i), property_id, f)
          call move_alloc(f, density_delta(i)%f)
        else
          call alloc_const_scalar_func(density_delta(i)%f, 0.0_r8)
        end if
      end do
    else
      do i = 1, size(fluids)
        call alloc_const_scalar_func(density_delta(i)%f, 0.0_r8)
      end do
    end if

    call flow_operators_init(this%mesh)
    call this%props%init(this%mesh, density, density_delta, viscosity)

    ! the initial velocity is provided from the velocity_init routine
    allocate(velocity_cc(3, this%mesh%mesh%ncell_onP))
    do i = 1,this%mesh%mesh%ncell_onP
      velocity_cc(:,i) = zone(i)%vc
    end do
    call this%flow%init(this%mesh, prescribed_flow, velocity_cc)

  end subroutine flow_driver_init

  subroutine flow_step(t, dt, vof, flux_vol, initial)
    use zone_module, only: Zone
    use physics_module, only: prescribed_flow
    use advection_velocity_namelist, only: adv_vel
    use truchas_timers
    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vof(:,:), flux_vol(:,:)
    logical, optional, intent(in) :: initial

    call start_timer('Flow')
    if (prescribed_flow) then
      block
        integer :: j
        real(r8) :: args(0:3)
        args(0) = t
        do j = 1, this%mesh%mesh%ncell
          args(1:3) = this%mesh%cell_centroid(:,j)
          this%flow%vel_cc(:,j) = adv_vel%eval(args)
        end do
        do j = 1, this%mesh%mesh%nface
          args(1:3) = this%mesh%face_centroid(:,j)
          this%flow%vel_fn(j) = dot_product(adv_vel%eval(args), &
              this%mesh%mesh%normal(:,j))/this%mesh%mesh%area(j)
        end do
      end block
      ! This will be needed by timestep
      this%temperature_cc(1:this%mesh%mesh%ncell_onP) = Zone%Temp
      call gather_boundary(this%mesh%mesh%cell_ip, this%temperature_cc)
      call this%props%update_cc(vof, this%temperature_cc, initial=initial)
      call stop_timer('Flow')
      return
    end if
#if ASDF
    associate (m => this%mesh%mesh)
      write(*, "('TOP LEVEL flux_vol[',i3,']: ',6es15.5): ") 771, flux_vol(1,m%xcface(771):m%xcface(772)-1)
    end associate
#endif
    this%temperature_cc(1:this%mesh%mesh%ncell_onP) = Zone%Temp
    call gather_boundary(this%mesh%mesh%cell_ip, this%temperature_cc)

    call this%props%update_cc(vof, this%temperature_cc, initial=initial)
    call this%flow%step(t, dt, this%props, flux_vol, initial=initial)
    call stop_timer('Flow')

  end subroutine flow_step

  subroutine flow_accept()
    call this%props%accept()
    call this%flow%accept()
  end subroutine flow_accept

end module flow_driver
