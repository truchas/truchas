!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_model_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use alloy_model_type
  use unstr_mesh_type
  use mfd_disc_type
  use matl_mesh_func_type
  use thermal_bc_factory_class
  use thermal_source_factory_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  public :: create_alloy_model

contains

  function create_alloy_model(tinit, disc, mmf, bc_fac, src_fac, params, stat, errmsg) result(model)

    use enclosure_radiation_namelist, only: er_params => params

    real(r8), intent(in) :: tinit
    type(mfd_disc), intent(in), target :: disc
    type(matl_mesh_func), intent(in), target :: mmf
    class(thermal_bc_factory), intent(inout) :: bc_fac
    type(thermal_source_factory), intent(inout) :: src_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable , intent(out):: errmsg
    type(alloy_model), pointer :: model

    allocate(model)

    call model%init(disc)

    !! Defines the equation parameter components of MODEL.
    call define_system_parameters(disc%mesh, mmf, model, stat, errmsg)
    if (stat /= 0) return

    !! Defines the boundary condition components of MODEL.
    call define_system_bc(disc%mesh, model, stat, errmsg)
    if (stat /= 0) return

    !! Define the alloy component of MODEL
    block
      use material_model_driver, only: matl_model
      use material_class
      integer :: n
      class(material), pointer :: matl
      character(:), allocatable :: name
      call params%get('material', name, stat, errmsg)
      if (stat /= 0) return
      n = matl_model%matl_index(name)
      if (n == 0) then
        stat = 1
        errmsg = 'unknown material: ' // name
        return
      end if
      call matl_model%get_matl_ref(n, matl)
      call model%alloy%init(matl, params, stat, errmsg)
      if (stat /= 0) return
    end block

  contains

    subroutine define_system_parameters (mesh, mmf, model, stat, errmsg)

      use matl_mesh_func_type
      use parallel_communication, only: global_any
      use material_model_driver, only: matl_model
      use material_utilities

      type(unstr_mesh), intent(in), target :: mesh
      type(matl_mesh_func), intent(in), target :: mmf
      type(alloy_model), intent(inout) :: model
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      !integer, allocatable :: matid(:)
      character(128) :: errmsg2

      !! Retrieve a list of all the material IDs that may be encountered.
      !call mmf%get_all_matl(matid, drop_void=.true.)

      !! Enthalpy density.
      call required_property_check(matl_model, 'enthalpy', stat, errmsg)
      if (stat /= 0) return
      call model%H_of_T%init(mmf, 'enthalpy', stat, errmsg2)
      if (global_any(stat /= 0)) then
        stat = -1
        errmsg = 'unexpected error defining H_of_T: ' // trim(errmsg2)
        return
      end if

      !! Thermal conductivity.
      call required_property_check(matl_model, 'conductivity', stat, errmsg)
      if (stat /= 0) return
      call model%conductivity%init(mmf, 'conductivity', stat, errmsg2)
      !call pmf_set_harmonic_average (model%conductivity)
      if (global_any(stat /= 0)) then
        stat = -1
        errmsg = 'unexpected error defining conductivity: ' // trim(errmsg2)
        return
      end if

      !! External heat source.
      call src_fac%alloc_source_funcs(model%src, stat, errmsg)
      if (stat /= 0) return

      stat = 0
      errmsg = ''

    end subroutine define_system_parameters

    subroutine define_system_bc(mesh, model, stat, errmsg)

      use evaporation_namelist, only: evap_params => params
      use evap_heat_flux_type
      use bitfield_type
      use parallel_communication, only: global_all, global_any, global_count
      use string_utilities, only: i_to_c

      type(unstr_mesh), intent(in), target :: mesh
      type(alloy_model), intent(inout) :: model
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      integer :: j, n
      logical, allocatable :: mask(:), rmask(:), fmask(:)
      integer, allocatable :: setids(:)
      character(160) :: string1, string2
      type(bitfield) :: bitmask
      type(evap_heat_flux), allocatable :: evap_flux

      allocate(mask(mesh%nface))

      mask = .false.

      !! Define the internal HTC interface conditions.
      call bc_fac%alloc_htc_ic(model%ic_htc, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%ic_htc)) then
        mask(model%ic_htc%index(1,:)) = .true.
        mask(model%ic_htc%index(2,:)) = .true.
        call mesh%face_imap%gather_offp(mask)
      end if

      !! Define the gap radiation interface conditions;
      !! may be superimposed with HTC conditions.
      call bc_fac%alloc_rad_ic(model%ic_rad, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%ic_rad)) then
        mask(model%ic_rad%index(1,:)) = .true.
        mask(model%ic_rad%index(2,:)) = .true.
        call mesh%face_imap%gather_offp(mask)
      end if

      !! Flux-type boundary conditions.  These may be superimposed.
      allocate(fmask(mesh%nface))
      fmask = .false.

      !! Define the simple flux boundary conditions.
      call bc_fac%alloc_flux_bc(model%bc_flux, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%bc_flux)) then
        if (global_any(mask(model%bc_flux%index))) then
          stat = -1
          errmsg = 'temperature flux boundary condition overlaps with interface conditions'
          return
        end if
        fmask(model%bc_flux%index) = .true. ! mark the simple flux faces
      end if

      !! Define the oriented flux boundary conditions.
      call bc_fac%alloc_vflux_bc(model%bc_vflux, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%bc_vflux)) then
        if (global_any(mask(model%bc_vflux%index))) &
            call TLS_info('    NOTE: oriented-flux condition is superimposed with interface conditions')
        do j = 1, size(model%bc_vflux%index)  ! index not necessarily 1-1
          fmask(model%bc_vflux%index(j)) = .true. ! mark the oriented flux faces
        end do
      end if

      !! Define the external HTC boundary conditions.
      call bc_fac%alloc_htc_bc(model%bc_htc, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%bc_htc)) then
        if (global_any(mask(model%bc_htc%index))) then
          stat = -1
          errmsg = 'temperature HTC boundary condition overlaps with interface conditions'
          return
        end if
        fmask(model%bc_htc%index) = .true. ! mark the HTC faces
      end if

      allocate(rmask(mesh%nface))
      rmask = .false.

      !! Define the (simple) radiation boundary conditions.
      call bc_fac%alloc_rad_bc(model%bc_rad, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%bc_rad)) then
        if (global_any(rmask(model%bc_rad%index))) then
          stat = -1
          errmsg = 'temperature radiation boundary condition overlaps with enclosure radiation'
          return
        end if
      end if
      if (allocated(model%bc_rad)) then
        if (global_any(mask(model%bc_rad%index))) then
          stat = -1
          errmsg = 'temperature radiation boundary condition overlaps with interface conditions'
          return
        end if
        fmask(model%bc_rad%index) = .true. ! mark the radiation faces
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
        call move_alloc(evap_flux, model%evap_flux)
      end if

      !! Merge flux mask into main mask.
      mask = mask .or. fmask
      deallocate(fmask)

      !! Define the Dirichlet boundary conditions.
      call bc_fac%alloc_dir_bc(model%bc_dir, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%bc_dir)) then
        if (global_any(mask(model%bc_dir%index))) then
          stat = -1
          errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
          return
        end if
        mask(model%bc_dir%index) = .true. ! mark the dirichlet faces
      end if

      mask = mask .or. rmask

      !! Finally verify that a condition has been applied to every boundary face.
      mask = mask .neqv. btest(mesh%face_set_mask,0)
      if (global_any(mask)) then
        call mesh%get_face_set_ids(pack([(j,j=1,mesh%nface)], mask), setids)
        if (size(setids) == 0) then
          string1 = '(none)'
        else
          write(string1,'(i0,*(:,", ",i0))') setids
        end if
        call mesh%get_link_set_ids(mask, setids)
        if (size(setids) == 0) then
          string2 = '(none)'
        else
          write(string2,'(i0,*(:,", ",i0))') setids
        end if
        errmsg = 'incomplete temperature boundary/interface specification;' // &
            ' remaining boundary faces belong to face sets ' // trim(string1) // &
            '; and interface link sets ' // trim(string2)
        bitmask = ibset(ZERO_BITFIELD, 0)
        mask = mask .and. (mesh%face_set_mask == bitmask)
        mask(mesh%lface(1,:)) = .false.
        mask(mesh%lface(2,:)) = .false.
        n = global_count(mask(:mesh%nface_onP))
        if (n > 0) errmsg = errmsg // '; ' // i_to_c(n) // ' faces belong to neither'
        stat = -1
        return
      end if

      stat = 0
      errmsg = ''

    end subroutine define_system_bc

  end function create_alloy_model

end module alloy_model_factory
