!!
!! THERMAL_COMPONENT_FACTORY
!!
!! This module defines the factory procedure that constructs a THERMAL_COMPONENT
!! for the enthalpy transport part of the HT, HTSD, and FHT models. It gathers
!! the thermal material properties, source terms, and standard thermal boundary
!! and interface conditions from the mesh, material, and input-parameter
!! factories; model-level coupling, such as view factor radiation, is handled
!! separately.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module thermal_component_factory

  use thermal_component_type
  use prop_cell_func_type, only: prop_cell_func
  use unstr_mesh_type
  use matl_mesh_func_type
  use thermal_bc_factory_type
  use thermal_source_factory_type
  implicit none
  private

  public :: define_thermal_component

contains

  subroutine define_thermal_component(mesh, mmf, bc_fac, src_fac, comp, stat, errmsg)

    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(thermal_bc_factory), intent(inout) :: bc_fac
    type(thermal_source_factory), intent(inout) :: src_fac
    type(thermal_component), intent(out) :: comp
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call define_thermal_parameters(mesh, mmf, comp, src_fac, stat, errmsg)
    if (stat /= 0) return

    call define_thermal_bc(mesh, comp, bc_fac, stat, errmsg)
    if (stat /= 0) return

  end subroutine

  subroutine define_thermal_parameters(mesh, mmf, comp, src_fac, stat, errmsg)

    use parallel_communication, only: global_any
    use material_model_driver, only: matl_model
    use material_utilities

    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(thermal_component), intent(inout) :: comp
    type(thermal_source_factory), intent(inout) :: src_fac
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    !! Enthalpy density.
    call required_property_check(matl_model, 'enthalpy', stat, errmsg)
    if (stat /= 0) return
    allocate(prop_cell_func :: comp%H_of_T)
    select type (H_of_T => comp%H_of_T)
    type is (prop_cell_func)
      call H_of_T%init(mmf, 'enthalpy', stat, errmsg)
    end select
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = 'unexpected error defining H_of_T: ' // errmsg
      return
    end if

    !! Thermal conductivity.
    call required_property_check(matl_model, 'conductivity', stat, errmsg)
    if (stat /= 0) return
    allocate(prop_cell_func :: comp%conductivity)
    select type (conductivity => comp%conductivity)
    type is (prop_cell_func)
      call conductivity%init(mmf, 'conductivity', stat, errmsg)
    end select
    !call pmf_set_harmonic_average (comp%conductivity)
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = 'unexpected error defining conductivity: ' // errmsg
      return
    end if

    !! User-defined heat source.
    call src_fac%alloc_source_funcs(comp%src, stat, errmsg)
    if (stat /= 0) return

    stat = 0

  end subroutine define_thermal_parameters

  subroutine define_thermal_bc(mesh, comp, bc_fac, stat, errmsg)

    use evaporation_namelist, only: evap_params => params
    use evap_heat_flux_type
    use parallel_communication, only: global_all, global_any
    use truchas_logging_services

    type(unstr_mesh), intent(in), target :: mesh
    type(thermal_component), intent(inout) :: comp
    type(thermal_bc_factory), intent(inout) :: bc_fac
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical, allocatable :: mask(:), fmask(:)
    type(evap_heat_flux), allocatable :: evap_flux

    allocate(mask(mesh%nface))

    mask = .false.

    !! Define the internal HTC interface conditions.
    call bc_fac%alloc_htc_ic(comp%ic_htc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%ic_htc)) then
      mask(comp%ic_htc%index(1,:)) = .true.
      mask(comp%ic_htc%index(2,:)) = .true.
    end if

    !! Define the gap radiation interface conditions;
    !! may be superimposed with HTC conditions.
    call bc_fac%alloc_rad_ic(comp%ic_rad, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%ic_rad)) then
      mask(comp%ic_rad%index(1,:)) = .true.
      mask(comp%ic_rad%index(2,:)) = .true.
    end if
    if (allocated(comp%ic_htc) .or. allocated(comp%ic_rad)) then
      call mesh%face_imap%gather_offp(mask)
    end if

    !! Flux-type boundary conditions.  These may be superimposed.
    allocate(fmask(mesh%nface))
    fmask = .false.

    !! Define the simple flux boundary conditions.
    call bc_fac%alloc_flux_bc(comp%bc_flux, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%bc_flux)) then
      if (global_any(mask(comp%bc_flux%index))) &
          call TLS_info('    NOTE: flux condition is superimposed with interface conditions')
      fmask(comp%bc_flux%index) = .true. ! mark the simple flux faces
    end if

    !! Define the oriented flux boundary conditions.
    call bc_fac%alloc_vflux_bc(comp%bc_vflux, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%bc_vflux)) then
      if (global_any(mask(comp%bc_vflux%index))) &
          call TLS_info('    NOTE: oriented-flux condition is superimposed with interface conditions')
      fmask(comp%bc_vflux%index) = .true. ! mark the oriented flux faces
    end if

    !! Define the external HTC boundary conditions.
    call bc_fac%alloc_htc_bc(comp%bc_htc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%bc_htc)) then
      if (global_any(mask(comp%bc_htc%index))) &
          call TLS_info('    NOTE: htc condition is superimposed with interface conditions')
      fmask(comp%bc_htc%index) = .true. ! mark the HTC faces
    end if

    !! Define the (simple) radiation boundary conditions.
    call bc_fac%alloc_rad_bc(comp%bc_rad, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%bc_rad)) then
      if (global_any(mask(comp%bc_rad%index))) &
          call TLS_info('    NOTE: radiation condition is superimposed with interface conditions')
      fmask(comp%bc_rad%index) = .true. ! mark the radiation faces
    end if

    !! Define the evaporation heat flux boundary condition
    if (allocated(evap_params)) then
      allocate(evap_flux)
      call evap_flux%init(mesh, evap_params, stat, errmsg)
      if (stat /= 0) return
      if (.not.global_all(fmask(evap_flux%index))) then
        stat = -1
        errmsg = 'evaporation heat flux applied to non-flux boundary'
        return
      end if
      call move_alloc(evap_flux, comp%evap_flux)
    end if

    !! Merge flux mask into main mask.
    mask = mask .or. fmask
    deallocate(fmask)

    !! Define the Dirichlet boundary conditions.
    call bc_fac%alloc_dir_bc(comp%bc_dir, stat, errmsg)
    if (stat /= 0) return
    if (allocated(comp%bc_dir)) then
      if (global_any(mask(comp%bc_dir%index))) then
        stat = -1
        errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
        return
      end if
      mask(comp%bc_dir%index) = .true. ! mark the dirichlet faces
    end if

    stat = 0

  end subroutine define_thermal_bc

end module thermal_component_factory
