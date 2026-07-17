!!
!! SPECIES_COMPONENT_FACTORY
!!
!! This module defines the factory procedure that constructs a SPECIES_COMPONENT
!! for the species transport part of the SD and HTSD models. It gathers the
!! species material properties, source term, and standard species boundary and
!! interface conditions from the mesh, material, and input-parameter factories.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module species_component_factory

  use species_component_type
  use prop_cell_func_type, only: prop_cell_func
  use unstr_mesh_type
  use matl_mesh_func_type
  use species_bc_factory_type
  use species_source_factory_type
  implicit none
  private

  public :: define_species_component

contains

  subroutine define_species_component(mesh, mmf, bc_fac, src_fac, species_index, species, &
      stat, errmsg, include_soret)

    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(species_bc_factory), intent(inout) :: bc_fac
    type(species_source_factory), intent(inout) :: src_fac
    integer, intent(in) :: species_index
    type(species_component), intent(out) :: species
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: include_soret

    call define_species_parameters(mmf, src_fac, species_index, species, stat, errmsg, include_soret)
    if (stat /= 0) return

    call define_species_bc(mesh, bc_fac, species_index, species, stat, errmsg)
    if (stat /= 0) return

  end subroutine define_species_component

  subroutine define_species_parameters(mmf, src_fac, species_index, species, &
      stat, errmsg, include_soret)

    use parallel_communication, only: global_any
    use material_model_driver, only: matl_model
    use material_utilities

    type(matl_mesh_func), intent(in), target :: mmf
    type(species_source_factory), intent(inout) :: src_fac
    integer, intent(in) :: species_index
    type(species_component), intent(inout) :: species
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: include_soret

    logical :: define_soret
    character(len=31) :: property

    ASSERT(species_index > 0)
    define_soret = .false.
    if (present(include_soret)) define_soret = include_soret

    !! Define the equation parameter components.
    write(property,'(a,i0)') 'diffusivity', species_index
    call required_property_check(matl_model, property, stat, errmsg)
    if (stat /= 0) return
    allocate(prop_cell_func :: species%diffusivity)
    select type (diffusivity => species%diffusivity)
    type is (prop_cell_func)
      call diffusivity%init(mmf, property, stat, errmsg)
    end select
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = 'unexpected error defining diffusivity: ' // errmsg
      return
    end if

    !! User-defined source for this species component.
    call src_fac%alloc_source_func(species_index, species%src, stat, errmsg)
    if (stat /= 0) return

    !! Define the optional Soret effect coefficients.
    if (define_soret) then
      write(property,'(a,i0)') 'soret-coef', species_index
      call optional_property_check(matl_model, property, stat, errmsg)
      if (stat < 0) return
      if (stat == 0) then
        allocate(prop_cell_func :: species%soret)
        select type (soret => species%soret)
        type is (prop_cell_func)
          call soret%init(mmf, property, stat, errmsg)
        end select
        if (global_any(stat /= 0)) then
          stat = -1
          errmsg = 'unexpected error defining soret coefficient: ' // errmsg
          return
        end if
      end if
    end if

    stat = 0

  end subroutine define_species_parameters

  subroutine define_species_bc(mesh, bc_fac, species_index, species, stat, errmsg)

    use bitfield_type, only: btest
    use parallel_communication, only: global_any

    type(unstr_mesh), intent(in), target :: mesh
    type(species_bc_factory), intent(inout) :: bc_fac
    integer, intent(in) :: species_index
    type(species_component), intent(inout) :: species
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j
    logical :: mask(mesh%nface)
    character(len=31) :: variable
    integer, allocatable :: setids(:)

    ASSERT(species_index > 0)

    !! Define the boundary condition components.
    write(variable,'(a,i0)') 'species-', species_index
    mask = .false.  ! used to tag faces where a BC has been applied

    !! Define the internal MTC interface conditions.
    call bc_fac%alloc_mtc_ic(species_index, species%ic_mtc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(species%ic_mtc)) then
      mask(species%ic_mtc%index(1,:)) = .true.
      mask(species%ic_mtc%index(2,:)) = .true.
      call mesh%face_imap%gather_offp(mask)
    end if

    !! Define the Dirichlet BC object for this concentration component.
    call bc_fac%alloc_dir_bc(species_index, species%bc_dir, stat, errmsg)
    if (stat /= 0) return
    if (allocated(species%bc_dir)) then
      if (global_any(mask(species%bc_dir%index))) then
        stat = -1
        errmsg = 'dirichlet BC overlaps with other BC'
        return
      end if
      mask(species%bc_dir%index) = .true.  ! tag the Dirichlet BC faces
    end if

    !! Define the flux BC object for this species component.
    call bc_fac%alloc_flux_bc(species_index, species%bc_flux, stat, errmsg)
    if (stat /= 0) return
    if (allocated(species%bc_flux)) then
      if (global_any(mask(species%bc_flux%index))) then
        stat = -1
        errmsg = trim(variable) // ': flux BC overlaps with other BC!'
        return
      end if
      mask(species%bc_flux%index) = .true.  ! tag the flux faces
    end if

    !! Define the MTC BC object for this species component.
    call bc_fac%alloc_mtc_bc(species_index, species%bc_mtc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(species%bc_mtc)) then
      if (global_any(mask(species%bc_mtc%index))) then
        stat = -1
        errmsg = trim(variable) // ': MTC BC overlaps with other BC!'
        return
      end if
      mask(species%bc_mtc%index) = .true.  ! tag the MTC faces
    end if

    !! Finally verify that a BC has been applied to every boundary face.
    if (global_any(mask.neqv.btest(mesh%face_set_mask,0))) then
      mask = mask .neqv. btest(mesh%face_set_mask,0)
      call mesh%get_face_set_IDs(pack([(j,j=1,mesh%nface)], mask), setids)
      stat = -1
      !! The assumption here is that no bad faces are internal; any errors of
      !! that kind should have been caught when creating the bd_data objects.
      if (allocated(errmsg)) deallocate(errmsg)
      allocate(character(len=2048) :: errmsg)
      write(errmsg,'(a,99(:,1x,i0))') &
        trim(variable) // ': incomplete BC specification; ' // &
        'remaining boundary faces belong to face sets', setids
      errmsg = trim(errmsg)
      deallocate(setids)
      return
    end if

    stat = 0

  end subroutine define_species_bc

end module species_component_factory
