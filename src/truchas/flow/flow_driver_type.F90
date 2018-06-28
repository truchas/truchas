#include "f90_assert.fpp"

module flow_driver_type

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
  implicit none
  private

  public :: read_flow_namelist, read_flow_params
  public :: flow_driver_init, flow_step, flow_final, flow_enabled

  type :: flow_driver
    type(flow_mesh), pointer :: mesh => null()
    type(flow) :: flow
    type(flow_props) :: props
  end type flow_driver
  type(flow_driver), allocatable :: this

contains

  subroutine flow_final
    if (allocated(this)) deallocate(this)
  end subroutine flow_final


  logical function flow_enabled()
    flow_enabled = allocated(this)
  end function flow_enabled

  subroutine read_flow_params(p)
    type(parameter_list), pointer, intent(in) :: p
    type(parameter_list), pointer :: pl_proj, pl_props

    allocate(this)
    pl_props => p%sublist("properties")
    pl_proj => p%sublist("projection")
    call this%props%read_params(pl_props)
    call this%flow%read_params(pl_proj)
  end subroutine read_flow_params

  subroutine read_flow_namelist(lun)
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities

    integer, intent(in) :: lun
    integer :: ios
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: p, pl_props, pl_proj, pl_fischer, pl_solver
    real(r8) :: fluid_cutvof, min_face_fraction
    integer :: fischer_history
    !- hypre related items
    real(r8) :: rel_tol, abs_tol, conv_rate_tol, amg_strong_threshold
    integer :: max_ds_iter, max_amg_iter, gmres_krylov_dim, cg_use_two_norm, logging_level
    integer :: print_level, amg_max_levels, amg_coarsen_type, amg_coarsen_sweeps
    integer :: amg_smoothing_method, amg_smoothing_sweeps, amg_interp_method

    character(:), allocatable :: krylov_method, amg_coarsen_method

    namelist /FLOW_CUTOFFS/ fluid_cutvof, min_face_fraction
    namelist /FLOW_PROJECTION/ fischer_history, rel_tol, &
        abs_tol, conv_rate_tol, max_ds_iter, max_amg_iter, &
        krylov_method, gmres_krylov_dim, cg_use_two_norm, &
        logging_level, print_level, amg_strong_threshold, &
        amg_max_levels, amg_coarsen_method, amg_coarsen_type, &
        amg_smoothing_sweeps, amg_smoothing_method, amg_interp_method

    allocate(p)
    pl_props => p%sublist("properties")
    pl_proj => p%sublist("projection")
    pl_fischer => pl_proj%sublist("fischer")
    pl_solver => pl_proj%sublist("solver")

    ! handle cutoffs
    fluid_cutvof = NULL_R
    min_face_fraction = NULL_R

    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'FLOW_CUTOFFS', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('')
      call TLS_info('Reading FLOWPROPS namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        read(lun,nml=flow_cutoffs,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOWPROPS namelist: ' // trim(iom))

      call broadcast(fluid_cutvof)
      call broadcast(min_face_fraction)
      if (fluid_cutvof /= NULL_R) call pl_props%set('cutvof', fluid_cutvof)
      if (min_face_fraction /= NULL_R) call pl_props%set('min_face_fraction', min_face_fraction)
    end if

    ! projection
    fischer_history = NULL_I
    max_ds_iter = NULL_I
    max_amg_iter = NULL_I
    gmres_krylov_dim = NULL_I
    cg_use_two_norm = NULL_I
    logging_level = NULL_I
    print_level = NULL_I
    amg_max_levels = NULL_I
    amg_coarsen_type = NULL_I
    amg_coarsen_sweeps = NULL_I
    amg_smoothing_method = NULL_I
    amg_interp_method = NULL_I
    rel_tol = NULL_R
    abs_tol = NULL_R
    conv_rate_tol = NULL_R
    amg_strong_threshold = NULL_R

    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'FLOW_PROJECTION', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('')
      call TLS_info('Reading FLOW_PROJECTION namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        read(lun,nml=flow_projection,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOW_PROJECTION namelist: ' // trim(iom))

      call broadcast(fischer_history)
      call broadcast(max_ds_iter)
      call broadcast(max_amg_iter)
      call broadcast(gmres_krylov_dim)
      call broadcast(cg_use_two_norm)
      call broadcast(logging_level)
      call broadcast(print_level)
      call broadcast(amg_max_levels)
      call broadcast(amg_coarsen_type)
      call broadcast(amg_coarsen_sweeps)
      call broadcast(amg_smoothing_method)
      call broadcast(amg_interp_method)
      call broadcast(rel_tol)
      call broadcast(abs_tol)
      call broadcast(conv_rate_tol)
      call broadcast(amg_strong_threshold)

      if (fischer_history /= NULL_I) call pl_fischer%set('history', fischer_history)
      if (rel_tol /= NULL_R) call pl_solver%set('rel-tol', rel_tol)
      if (abs_tol /= NULL_R) call pl_solver%set('abs-tol', abs_tol)
      if (conv_rate_tol /= NULL_R) call pl_solver%set('conv-rate-tol', conv_rate_tol)
      if (max_ds_iter /= NULL_I) call pl_solver%set('max-ds-iter', max_ds_iter)
      if (max_amg_iter /= NULL_I) call pl_solver%set('max-amg-iter', max_amg_iter)
      if (allocated(krylov_method)) call pl_solver%set('krylov-method', krylov_method)
      if (gmres_krylov_dim /= NULL_I) call pl_solver%set('gmres-krylov-dim', gmres_krylov_dim)
      if (cg_use_two_norm /= NULL_I) call pl_solver%set('cg-use-two-norm', cg_use_two_norm /= 0)
      if (logging_level /= NULL_I) call pl_solver%set('logging-level', logging_level)
      if (print_level /= NULL_I) call pl_solver%set('print-level', print_level)
      if (amg_strong_threshold /= NULL_R) call pl_solver%set('amg-strong-threshold', amg_strong_threshold)
      if (amg_max_levels /= NULL_I) call pl_solver%set('amg-max-levels', amg_max_levels)
      if (allocated(amg_coarsen_method)) call pl_solver%set('amg-coarsen-method', amg_coarsen_method)
      if (amg_coarsen_type /= NULL_I) call pl_solver%set('amg-coarsen-type', amg_coarsen_type)
      if (amg_smoothing_sweeps /= NULL_I) call pl_solver%set('amg-smoothing-sweeps', amg_smoothing_sweeps)
      if (amg_smoothing_method /= NULL_I) call pl_solver%set('amg-smoothing-method', amg_smoothing_method)
      if (amg_interp_method /= NULL_I) call pl_solver%set('amg-interp-method', amg_interp_method)

    end if

    call read_flow_params(p)
  end subroutine read_flow_namelist


  subroutine flow_driver_init(mesh)
    type(unstr_mesh), pointer, intent(in) :: mesh
    !-
    integer :: void, mu_id, rho_id, rho_delta_id, i
    integer, pointer :: phases(:)
    integer, allocatable :: fluids(:)
    real(r8), allocatable :: density(:)
    real(r8) :: state(1)
    class(scalar_func_box), allocatable :: density_delta(:), viscosity(:)
    class(scalar_func), allocatable :: f

    if (.not.allocated(this)) return

    call this%mesh%init(mesh)

    ! lots of duplication here from vtrack_driver.  This should all be subsumed and handled
    ! properly by a sufficiently intelligent physics driver at some point
    if (ppt_has_property("viscosity")) then
      mu_id = ppt_property_id("viscosity")
    else ! no fluids
      deallocate(this)
      return
    end if
    rho_id = ppt_property_id("density")
    rho_delta_id = ppt_property_id("density deviation")

    void = 0
    state(1) = 0.0_r8
    call ppt_get_phase_ids(phases)
    do i = 1, size(phases)
      ! all fluids have a viscosity and density
      if (ppt_has_phase_property(phases(i), mu_id) .and. &
          ppt_has_phase_property(phases(i), rho_id)) then
        ! void has zero density
        call ppt_get_phase_property(phases(i), rho_id, f)
        ASSERT(allocated(f))
        ASSERT(is_const(f))
        if (f%eval(state) == 0.0_r8) then
          void = void + 1
          cycle
        end if
        fluids = [fluids, phases(i)]
      end if
    end do

    allocate(density(size(fluids)))
    allocate(density_delta(size(fluids)))
    allocate(viscosity(size(fluids)))

    do i = 1, size(fluids)
      call ppt_get_phase_property(fluids(i), rho_id, f)
      density(i) = f%eval(state)
      call ppt_get_phase_property(fluids(i), mu_id, f)
      call move_alloc(f, viscosity(i)%f)
      call ppt_get_phase_property(fluids(i), rho_delta_id, f)
      call move_alloc(f, density_delta(i)%f)
    end do

    call this%props%init(this%mesh, density, density_delta, viscosity, void > 0)
  end subroutine flow_driver_init

  subroutine flow_step(t, dt, vof, temperature_cc)
    real(r8), intent(in) :: t, dt, vof(:,:), temperature_cc(:)
    !-

    call this%props%update_cc(vof, temperature_cc, initial=.false.)
    call this%flow%step(...)

  end subroutine flow_step

  subroutine flow_accent()

  end subroutine flow_accent

end module flow_driver_type
