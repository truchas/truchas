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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
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

module flow_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use flow_type
  use flow_props_type
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
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
  public :: flow_driver_set_initial_state
  public :: flow_driver_dump_state
  public :: flow_set_pre_solidification_density
  public :: inflow_plist

  type :: flow_driver_data
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- not owned
    type(flow) :: flow
    type(flow_props) :: props
    ! The flow driver shouldn't logically need this but temperature is currently
    ! stored in Zone%Temp so we need to keep a version on the new mesh as well
    real(r8), allocatable :: temperature_cc(:)

    ! a copy of the diffusion solver's variable for surface tension to hold as a reference
    real(r8), allocatable :: temperature_fc(:)

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

  subroutine flow_timestep(dt, dt_tag)
    real(r8), intent(out) :: dt
    character(:), allocatable, intent(inout) :: dt_tag
    ASSERT(flow_enabled())
    call this%flow%timestep(dt, dt_tag)
  end subroutine flow_timestep

  logical function flow_enabled()
    flow_enabled = allocated(this)
  end function flow_enabled

  function inflow_plist()
    type(parameter_list), pointer :: inflow_plist
    inflow_plist => this%flow%bc%inflow_plist
  end function inflow_plist

  subroutine read_flow_namelists(lun)

    use flow_namelist
    use flow_solver_namelists
    use flow_bc_namelist, only: read_flow_bc_namelists
    use turbulence_namelist, only: read_turbulence_namelist
    use physics_module, only: body_force_density, prescribed_flow
    use advection_velocity_namelist

    integer, intent(in) :: lun

    type(parameter_list), pointer :: plist

    if (.not.allocated(this)) allocate(this)

    call read_flow_namelist(lun)

    ! Insert body_force_density from the PHYSICS namelist
    plist => params%sublist('options')
    call plist%set('body force', body_force_density)

    ! Add additional sublists from other namelists
    plist => params%sublist('bc')
    call read_flow_bc_namelists(lun, plist)
    plist => params%sublist("predictor")
    call read_flow_viscous_solver_namelist(lun, plist)
    plist => params%sublist('corrector')
    call read_flow_pressure_solver_namelist(lun, plist)

    plist => params%sublist('turbulence model')
    call read_turbulence_namelist(lun, plist)

    if (prescribed_flow) call read_advection_velocity_namelist(lun)

  end subroutine read_flow_namelists


  subroutine flow_driver_init

    use mesh_manager, only: unstr_mesh_ptr
    use flow_namelist, only: params
    use material_model_driver, only: matl_model
    use physics_module, only: prescribed_flow
    use scalar_func_factories, only: alloc_const_scalar_func
    use vtrack_driver, only: vtrack_driver_init, vtrack_set_inflow_bc
    use flow_bc_type
    use truchas_logging_services

    integer :: n, i, stat
    integer, allocatable :: fluids(:)
    real(r8), allocatable :: density(:)
    class(scalar_func_box), allocatable :: density_delta(:), viscosity(:)
    class(scalar_func), allocatable :: f
    type(parameter_list), pointer :: plist
    type(flow_bc), pointer :: flowbc
    character(:), allocatable :: errmsg

    INSIST(allocated(this))

    plist => params%sublist('volume-tracker')
    call vtrack_driver_init(plist)

    call TLS_info('')
    call TLS_info('Initializing fluid flow solver ...')

    this%mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(this%mesh))
    call this%mesh%init_cell_centroid
    call this%mesh%init_face_centroid
    call this%mesh%init_face_normal_dist

    allocate(this%temperature_cc(this%mesh%ncell), this%temperature_fc(this%mesh%nface))

    ! Some duplication here from vtrack_driver.  This should all be subsumed and handled
    ! properly by a sufficiently intelligent physics driver at some point
    do i = matl_model%nphase_real, 1, -1
      if (matl_model%is_fluid(i)) exit
    end do
    if (i == 0) then ! no fluids
      deallocate(this)
      return
    end if

    ! identify non-void fluids
    n = matl_model%nphase_real
    fluids = pack([(i,i=1,n)], matl_model%is_fluid(:n))

    ! get density
    allocate(density(size(fluids)))
    do i = 1, size(fluids)
      call matl_model%get_phase_prop(fluids(i), 'density', f)
      ASSERT(allocated(f))
      ASSERT(is_const(f))
      density(i) = f%eval([real(r8)::])
    end do

    ! get viscosity
    ! if not given, assume 0-valued viscosity
    allocate(viscosity(size(fluids)))
    do i = 1, size(fluids)
      call matl_model%get_phase_prop(fluids(i), 'viscosity', viscosity(i)%f)
      if (.not.allocated(viscosity(i)%f)) call alloc_const_scalar_func(viscosity(i)%f, 0.0_r8)
    end do

    ! get density deltas
    ! if not given, assume 0-valued density delta
    allocate(density_delta(size(fluids)))
    do i = 1, size(fluids)
      call matl_model%get_phase_prop(fluids(i), 'density-delta', density_delta(i)%f)
      if (.not.allocated(density_delta(i)%f)) call alloc_const_scalar_func(density_delta(i)%f, 0.0_r8)
    end do

    call flow_operators_init(this%mesh)
    i = size(fluids) + merge(1, 0, matl_model%have_void) ! get number of fluids
    call this%props%init(this%mesh, i, density, density_delta, viscosity, params)

    allocate(flowbc)
    if (.not.prescribed_flow) then
      call flowbc%init(this%mesh, params, this%props%vof, this%temperature_fc, stat, errmsg)
      if (stat /= 0) call TLS_fatal('FLOW initialization error: ' // errmsg)
      call vtrack_set_inflow_bc(flowbc%inflow_plist, stat, errmsg)
    end if

    call this%flow%init(this%mesh, this%props, flowbc, prescribed_flow, params=params)

  end subroutine flow_driver_init

  subroutine flow_driver_set_initial_state(t, dt, temperature_fc, vel_fn)

    use zone_module, only: zone
    use vtrack_driver, only: vtrack_vof_view
    use physics_module, only: prescribed_flow
    use advection_velocity_namelist, only: adv_vel

    real(r8), intent(in) :: t, dt
    real(r8), intent(in), pointer :: temperature_fc(:)
    real(r8), intent(in), optional :: vel_fn(:)

    integer :: j
    real(r8) :: vcell(3,this%mesh%ncell_onP)
    real(r8), pointer :: vof(:,:)

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
      vof => vtrack_vof_view()
      call this%mesh%cell_imap%gather_offp(this%temperature_cc)
      call this%props%set_initial_state(vof, this%temperature_cc)
      call stop_timer('Flow')
      return
    end if

    do j = 1, this%mesh%ncell_onP
      vcell(:,j) = zone(j)%vc
    end do

    this%temperature_cc(1:this%mesh%ncell_onP) = Zone%Temp
    call this%mesh%cell_imap%gather_offp(this%temperature_cc)

    if (associated(temperature_fc)) then
      this%temperature_fc(:this%mesh%nface_onP) = temperature_fc(:this%mesh%nface_onP)
      call this%mesh%face_imap%gather_offp(this%temperature_fc)
    end if

    vof => vtrack_vof_view()

    if (present(vel_fn)) then ! RESTART
      call this%flow%set_initial_state(t, zone%p, vcell, vel_fn, vof, this%temperature_cc)
    else
      call this%flow%set_initial_state(t, dt, vcell, vof, this%temperature_cc)
    end if

    call stop_timer('Flow')

  end subroutine flow_driver_set_initial_state


  subroutine flow_step(t, dt, vof, flux_vol, temperature_fc)

    use zone_module, only: Zone
    use physics_module, only: prescribed_flow
    use advection_velocity_namelist, only: adv_vel

    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vof(:,:), flux_vol(:,:)
    real(r8), intent(in), pointer :: temperature_fc(:)

    call start_timer('Flow')

    this%temperature_cc(1:this%mesh%ncell_onP) = Zone%Temp
    call this%mesh%cell_imap%gather_offp(this%temperature_cc)

    if (prescribed_flow) then
      call this%props%update_cc(vof, this%temperature_cc)
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
    else
      if (associated(temperature_fc)) then
        this%temperature_fc(:this%mesh%nface_onP) = temperature_fc(:this%mesh%nface_onP)
        call this%mesh%face_imap%gather_offp(this%temperature_fc)
      end if
      call this%flow%step(t, dt, vof, flux_vol, this%temperature_cc)
    end if

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

    use degen_hex_topology, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP
    use restart_utilities, only: read_dist_array
    use mesh_manager, only: unstr_mesh_ptr

    integer, intent(in) :: unit, version
    real(r8), allocatable, intent(out) :: vel_fn(:)

    integer :: j, k, n
    real(r8), allocatable :: old_flux_vel(:,:)
    type(unstr_mesh), pointer :: mesh

    mesh => unstr_mesh_ptr('MAIN') !NB: cannot use this%mesh; it is not initialized when this subroutine is called!
    INSIST(associated(mesh))

    allocate(old_flux_vel(6,mesh%ncell))
    call read_dist_array(unit, old_flux_vel(:,:mesh%ncell_onP), mesh%xcell(:mesh%ncell_onP), 'READ_FLUXING_VELOCITY: error reading Fluxing_Velocity records')
    call mesh%cell_imap%gather_offp(old_flux_vel)

    allocate(vel_fn(mesh%nface_onP))
    vel_fn = 0

    do j = 1, mesh%nface_onP
      n = mesh%fcell(1,j)
      associate (cface => mesh%cface(mesh%xcface(n):mesh%xcface(n+1)-1))
        do k = 1, size(cface)
          if (cface(k) == j) exit
        end do
        ASSERT(k <= size(cface))
        select case (mesh%xcnode(n+1)-mesh%xcnode(n)) ! number of nodes for cell
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
        vel_fn(j) = old_flux_vel(k,n)
      end associate
    end do

  end subroutine read_fluxing_velocity

  !! This is the companion of read_fluxing_velocity. The HDF5 output still uses
  !! the legacy degenerate hex mesh format and cell-based outward fluxing
  !! velocities used by the legacy flow solver. This routine returns the new
  !! flow solvers face normal velocities in that format to be used in output.

  subroutine get_legacy_flux_vel(fluxing_velocity)

    use degen_hex_topology, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP

    real(r8), intent(out) :: fluxing_velocity(:,:)

    real(r8), pointer :: vel_fn(:)
    real(r8) :: tmp
    integer :: j, k

    ASSERT(size(fluxing_velocity,dim=1) == 6)
    ASSERT(size(fluxing_velocity,dim=2) == this%mesh%ncell_onP)

    vel_fn => flow_vel_fn_view()
    ASSERT(size(vel_fn) == this%mesh%nface)  ! we assume off-process are up-to-date
    fluxing_velocity = 0.0_r8
    do j = 1, this%mesh%ncell_onP
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

  subroutine flow_driver_dump_state
    call this%flow%dump_state
  end subroutine flow_driver_dump_state

  subroutine flow_set_pre_solidification_density(vof)
    real(r8), intent(in) :: vof(:,:)
    call this%props%set_pre_solidification_density(vof)
  end subroutine flow_set_pre_solidification_density

end module flow_driver
