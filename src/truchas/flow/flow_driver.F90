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
  implicit none
  private

  public :: read_flow_namelists, flow_timestep
  public :: flow_driver_init, flow_step, flow_final, flow_enabled, flow_accept
  public :: flow_vel_fn_view, flow_vel_cc_view, flow_P_cc_view
  public :: read_fluxing_velocity, get_legacy_flux_vel

  type :: flow_driver_data
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- not owned
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
    call this%flow%timestep(dtc, dtv)
  end subroutine flow_timestep

  logical function flow_enabled()
    flow_enabled = allocated(this)
  end function flow_enabled

  subroutine read_flow_namelists(lun)

    use flow_namelist
    use flow_predictor_namelist
    use flow_corrector_namelist
    use flow_bc_namelist, only: read_flow_bc_namelists
    use turbulence_namelist, only: read_turbulence_namelist
    use physics_module, only: body_force_density, prescribed_flow
    use advection_velocity_namelist

    integer, intent(in) :: lun
    integer :: ios
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: plist

    if (.not.allocated(this)) allocate(this)

    call read_flow_namelist(lun)

    ! Insert body_force_density from the PHYSICS namelist
    plist => params%sublist('options')
    call plist%set('body force', body_force_density)

    ! Add additional sublists from other namelists
    plist => params%sublist('bc')
    call read_flow_bc_namelists(lun, plist)
    call read_flow_predictor_namelist(lun, params)
    call read_flow_corrector_namelist(lun, params)

    plist => params%sublist('turbulence model')
    call read_turbulence_namelist(lun, plist)

    if (prescribed_flow) call read_advection_velocity_namelist(lun)

  end subroutine read_flow_namelists


  subroutine flow_driver_init(vel_fn)
    use mesh_manager, only: unstr_mesh_ptr
    use flow_namelist, only: params
    use zone_module, only: zone
    use material_interop, only: void_material_index, material_to_phase
    use physics_module, only: prescribed_flow
    use scalar_func_factories, only: alloc_const_scalar_func
    use property_data_module, only: isImmobile
    use parameter_module, only: nmat
    use vtrack_driver, only: vtrack_driver_init

    real(r8), intent(in), optional :: vel_fn(:)

    !real(r8), intent(in) :: vof(:,:), flux_vol(:,:)
    !-
    integer :: i, property_id
    integer, allocatable :: fluids(:)
    logical, allocatable :: is_real_fluid(:)
    real(r8), allocatable :: density(:), velocity_cc(:,:)
    real(r8) :: state(1)
    class(scalar_func_box), allocatable :: density_delta(:), viscosity(:)
    class(scalar_func), allocatable :: f
    type(parameter_list), pointer :: plist

    INSIST(allocated(this))

    plist => params%sublist('volume-tracker')
    call vtrack_driver_init(plist)

    this%mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(this%mesh))
    call this%mesh%init_cell_centroid
    call this%mesh%init_face_centroid
    call this%mesh%init_face_normal_dist

    allocate(this%temperature_cc(this%mesh%ncell))

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
    call this%props%init(this%mesh, density, density_delta, viscosity, params)

    ! the initial velocity is provided from the velocity_init routine
    allocate(velocity_cc(3, this%mesh%ncell_onP))
    do i = 1,this%mesh%ncell_onP
      velocity_cc(:,i) = zone(i)%vc
    end do
    !!FIXME? Optional argument P_CC is missing -- is it needed?
    call this%flow%init(this%mesh, this%props, prescribed_flow, velocity_cc, zone%p, vel_fn, params=params)

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
        do j = 1, this%mesh%ncell
          args(1:3) = this%mesh%cell_centroid(:,j)
          this%flow%vel_cc(:,j) = adv_vel%eval(args)
        end do
        do j = 1, this%mesh%nface
          args(1:3) = this%mesh%face_centroid(:,j)
          this%flow%vel_fn(j) = dot_product(adv_vel%eval(args), &
              this%mesh%normal(:,j))/this%mesh%area(j)
        end do
      end block
      ! This will be needed by timestep
      this%temperature_cc(1:this%mesh%ncell_onP) = Zone%Temp
      call gather_boundary(this%mesh%cell_ip, this%temperature_cc)
      call this%props%update_cc(vof, this%temperature_cc, initial=initial)
      call stop_timer('Flow')
      return
    end if
#if ASDF
    associate (m => this%mesh)
      write(*, "('TOP LEVEL flux_vol[',i3,']: ',6es15.5): ") 771, flux_vol(1,m%xcface(771):m%xcface(772)-1)
    end associate
#endif
    this%temperature_cc(1:this%mesh%ncell_onP) = Zone%Temp
    call gather_boundary(this%mesh%cell_ip, this%temperature_cc)

    call this%props%update_cc(vof, this%temperature_cc, initial=initial)
    call this%flow%step(t, dt, flux_vol, initial=initial)
    call stop_timer('Flow')

  end subroutine flow_step

  subroutine flow_accept()
    call this%props%accept()
    call this%flow%accept()
  end subroutine flow_accept

  !! This procedure is the counterpart to fluid_data_module::read_flow_data,
  !! which reads the fluxing velocity record from the restart file. The
  !! challenge here is that the HDF5 output and restart file use the legacy
  !! degenerate hex format, and for the fluxing velocities, holds the data as
  !! the outward normal side velocities for each cell. So we need to both map
  !! legacy mesh sides to new mesh sides, and then to faces. The data is
  !! essentially duplicated for interior faces (if it is consistent). We use
  !! data from the side whose outward orientation is the same as that for the
  !! face.

  subroutine read_fluxing_velocity(unit, version, vel_fn)

    use legacy_mesh_api, only: ncells, nfc, pcell => unpermute_mesh_vector
    use common_impl, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP
    use restart_utilities, only: read_dist_array
    use index_partitioning, only: gather_boundary

    integer, intent(in) :: unit, version
    real(r8), allocatable, intent(out) :: vel_fn(:)

    integer :: j, k, n
    real(r8) :: old_flux_vel(nfc,this%mesh%ncell)

    INSIST(this%mesh%ncell_onP == ncells) ! may fail if gap elements are present

    call read_dist_array(unit, old_flux_vel(:,:ncells), pcell, 'READ_FLUXING_VELOCITY: error reading Fluxing_Velocity records')
    call gather_boundary(this%mesh%cell_ip, old_flux_vel)

    allocate(vel_fn(this%mesh%nface_onP))
    vel_fn = 0

!    do j = 1, this%mesh%ncell
!      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
!        select case (this%mesh%xcnode(j+1)-this%mesh%xcnode(j)) ! number of nodes for cell
!        case (4)  ! tet cell
!          do k = 1, size(cface)
!            if (btest(this%mesh%cfpar(j),pos=k)) cycle  ! opposite orientation
!            vel_fn(cface(k)) = old_flux_vel(NEW_TET_SIDE_MAP(k),j)
!          end do
!        case (5)  ! pyramid cell
!          do k = 1, size(cface)
!            if (btest(this%mesh%cfpar(j),pos=k)) cycle  ! opposite orientation
!            vel_fn(cface(k)) = old_flux_vel(NEW_PYR_SIDE_MAP(k),j)
!          end do
!        case (6)  ! prism cell
!          do k = 1, size(cface)
!            if (btest(this%mesh%cfpar(j),pos=k)) cycle  ! opposite orientation
!            vel_fn(cface(k)) = old_flux_vel(NEW_PRI_SIDE_MAP(k),j)
!          end do
!        case (8)  ! hex cell
!          do k = 1, size(cface)
!            if (btest(this%mesh%cfpar(j),pos=k)) cycle  ! opposite orientation
!            vel_fn(cface(k)) = old_flux_vel(NEW_HEX_SIDE_MAP(k),j)
!          end do
!        case default
!          INSIST(.false.)
!        end select
!      end associate
!    end do
!    call gather_boundary(this%mesh%face_ip, vel_fn)

    do j = 1, this%mesh%nface_onP
      n = this%mesh%fcell(1,j)
      associate (cface => this%mesh%cface(this%mesh%xcface(n):this%mesh%xcface(n+1)-1))
        do k = 1, size(cface)
          if (cface(k) == j) exit
        end do
        ASSERT(k <= size(cface))
        select case (this%mesh%xcnode(n+1)-this%mesh%xcnode(n)) ! number of nodes for cell
        case (4)  ! tet cell
          k = NEW_TET_SIDE_MAP(k)
        case (5)  ! pyramid cell
          k = NEW_PYR_SIDE_MAP(k)
        case (6)  ! prism cell
          k = NEW_PRI_SIDE_MAP(k)
        case (8)  ! hex cell
          k = NEW_HEX_SIDE_MAP(k)
        case default
          INSIST(.false.)
        end select
        vel_fn(j) = old_flux_vel(k,j)
      end associate
    end do

  end subroutine read_fluxing_velocity

  !! This is the companion of read_fluxing_velocity. The HDF5 output still uses
  !! the legacy degenerate hex mesh format and cell-based outward fluxing
  !! velocities used by the legacy flow solver. This routine returns the new
  !! flow solvers face normal velocities in that format to be used in output.

  subroutine get_legacy_flux_vel(fluxing_velocity)

    use legacy_mesh_api, only: ncells, nfc
    use common_impl, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP

    real(r8), intent(out) :: fluxing_velocity(:,:)

    real(r8), pointer :: vel_fn(:)
    real(r8) :: tmp
    integer :: j, k

    ASSERT(size(fluxing_velocity,dim=1) == nfc)
    ASSERT(size(fluxing_velocity,dim=2) == ncells)
    INSIST(ncells == this%mesh%ncell_onP) ! may fail if gap elements are present

    vel_fn => flow_vel_fn_view()
    ASSERT(size(vel_fn) == this%mesh%nface)  ! we assume off-process are up-to-date
    fluxing_velocity = 0.0_r8
    do j = 1, ncells
      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        select case (this%mesh%xcnode(j+1)-this%mesh%xcnode(j)) ! number of nodes for cell
        case (4)  ! tet cell
          do k = 1, size(cface)
            tmp = vel_fn(cface(k))
            if (btest(this%mesh%cfpar(j),pos=k)) tmp = -tmp
            fluxing_velocity(NEW_TET_SIDE_MAP(k),j) = tmp
          end do
        case (5)  ! pyramid cell
          do k = 1, size(cface)
            tmp = vel_fn(cface(k))
            if (btest(this%mesh%cfpar(j),pos=k)) tmp = -tmp
            fluxing_velocity(NEW_PYR_SIDE_MAP(k),j) = tmp
          end do
        case (6)  ! prism cell
          do k = 1, size(cface)
            tmp = vel_fn(cface(k))
            if (btest(this%mesh%cfpar(j),pos=k)) tmp = -tmp
            fluxing_velocity(NEW_PRI_SIDE_MAP(k),j) = tmp
          end do
        case (8)  ! hex cell
          do k = 1, size(cface)
            tmp = vel_fn(cface(k))
            if (btest(this%mesh%cfpar(j),pos=k)) tmp = -tmp
            fluxing_velocity(NEW_HEX_SIDE_MAP(k),j) = tmp
          end do
        case default
          INSIST(.false.)
        end select
      end associate
    end do

  end subroutine get_legacy_flux_vel

end module flow_driver
