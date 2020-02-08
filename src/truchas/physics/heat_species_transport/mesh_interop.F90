!!
!! MESH_INTEROP
!!
!! This module provides some data and procedures that facilitate the exchange
!! of mesh-based data between the Truchas mesh structure and the new mesh data
!! structure used by the diffusion solver.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The following routines handle the mapping between the MATL module data
!! structure used in the older parts of Truchas and a material mesh function
!! used by the diffusion solver.  Both structures describe the distribution
!! of materials across the mesh, albeit in different ways.  Note that while
!! the procedures take a general material mesh function as an argument, they
!! necessarily use and modify the MATL module variable from MATL_MODULE.
!!
!!  CALL MMF_INIT (MESH, MMF, STAT, ERRMSG) creates a material mesh function
!!    MMF over MESH that is compatible with the MATL structure, and initializes
!!    its volume fraction data using the current VoF data from MATL.  Note
!!    that this implicitly calls UPDATE_MMF_FROM_MATL so a separate explicit
!!    call to it is not needed.
!!
!!  CALL UPDATE_MMF_FROM_MATL (MMF) copies the Truchas material VoF data
!!    stored in the MATL structure into the material system volume fraction
!!    data stored in the MMF structure.  Note that though a material may
!!    transform from one phase to another within the material system, the
!!    material itself is conserved.  Thus in the absence of flow, which could
!!    transport material between cells, the material system volume fractions
!!    remain fixed, and the work done by this routine only needs to be done
!!    once.
!!
!!  CALL UPDATE_MATL_FROM_MMF (MMF, STATE) updates the Truchas MATL module
!!    variable with new VoF data derived from phase mixtures as the current
!!    solution state and material system volume fractions stored in the MMF
!!    structure.  Note that although the diffusion solver does not alter the
!!    material system volume fractions (materials are conserved), the phase
!!    mixture of a multi-phase material will depend on the current value of
!!    the state variables and so the Truchas VoF data, which is really the
!!    volume fraction of phases, will need to be updated, but only when
!!    phase change is possible.
!!

#include "f90_assert.fpp"

!#define EXTRA_VOF_DIAGNOSTICS

module mesh_interop

  use kinds, only: r8
  use base_mesh_class
  use matl_mesh_func_type
  use material_model_driver, only: matl_model
  use truchas_logging_services, only: TLS_info
  implicit none
  private

  public :: mmf_init, update_mmf_from_matl, update_matl_from_mmf
  public :: void_is_present

contains

  logical function void_is_present ()
    void_is_present = matl_model%have_void
  end function void_is_present

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MMF_INIT
 !!
 !! This routine creates the material mesh function MMF and initializes its
 !! volume fraction data using the VoF data from Truchas' MATL data structure.
 !! The MMF describes the distribution of material systems across the mesh and
 !! forms a basis for defining the property mesh functions (e.g. conductivity)
 !! passed to the heat transfer/species diffusion solver.
 !!
 !! Although a MMF is capable of describing situations where material systems
 !! are restricted to certain mesh regions, such information is not currently
 !! available within the Truchas framework, and so this routine creates a MMF
 !! with a single region comprising the whole mesh that may contain any of
 !! the material systems that have been defined.  One consequence of this
 !! choice is exploited by UPDATE_MMF_FROM_MATL and UPDATE_MATL_FROM_MMF:
 !!
 !!   The list of region cells returned by MMF_REG_CELL will be the identity
 !!   mapping; i.e., the first region cell is cell 1, the seccond is cell 2,
 !!   and so forth.
 !!
 !! These two lists (MMF%REG_CELL, MMF%REG_MATL) are normally needed to
 !! associate the elements of the rank-2 array returned by MMF_REG_VOL_FRAC
 !! to a cell and material system, but in this case the array indices are
 !! themselves the cell number and material system ID.
 !!
 !! Prerequisites:
 !! A. The material table is fully populated and no materials deleted so that
 !!    the material IDs are contiguous from 1.  If any subsequent changes are
 !!    made to the material table, they will not be reflected in the MMF that
 !!    is initialized here.
 !! B. GENERATE_MATERIAL_MAPPINGS from the MATERIAL_INTEROP module has been
 !!    called to establish the mappings between Truchas material numbers and
 !!    material system IDs that are used in the MMF <--> MATL update routines.
 !!

  subroutine mmf_init (mesh, mmf, stat, errmsg)

    class(base_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(out) :: mmf
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    integer, allocatable :: all_materials(:)
    integer, pointer :: all_cell_sets(:)

    !! Single region consisting of the entire mesh allowing all materials.
    n = matl_model%nmatl_real
    if (matl_model%have_void) then
      allocate(all_materials(0:n))
      all_materials(0) = 0  ! Special ID for void -- no corresponding material
    else
      allocate(all_materials(1:n))
    end if
    all_materials(1:n) = [(j, j=1,n)]
    all_cell_sets => mesh%cell_set_id
    call mmf%init(mesh)
    call mmf%define_region(all_cell_sets, all_materials, stat, errmsg)
    deallocate(all_materials)
    if (stat /= 0) then
      errmsg = 'error defining MMF region: ' // trim(errmsg)
      return
    end if
    call mmf%define_complete(stat, errmsg)
    if (stat /= 0) then
      errmsg = 'error creating MMF: ' // trim(errmsg)
      return
    end if

    !! Assumed condition; see note.
    INSIST(all(mmf%reg_cell(1) == (/(j,j=1,mesh%ncell)/)))

    !! Finally initialize the MMF volume fraction data from the MATL VoF data.
    call update_mmf_from_matl (mmf)

    stat = 0

  end subroutine mmf_init

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! UPDATE_MMF_FROM_MATL
 !!
 !! This routine copies the Truchas material VoF data stored in the MATL
 !! structure into the material system volume fraction data stored in the
 !! MMF structure.  Note that though a material may transform from one phase
 !! to another within the material system, the material itself is conserved.
 !! Thus in the absence of flow, which could transport material between cells,
 !! the material system volume fractions remain fixed, and the work done by
 !! this routine only needs to be done once.
 !!
 !! The routine assumes a single-region MMF having an identity mapping between
 !! region cells and mesh cells, and uses the mapping MATERIAL_TO_SYSTEM from
 !! the MATERIAL_INTEROP module.
 !!

  subroutine update_mmf_from_matl (mmf)

!    use material_class
    use matl_module, only: gather_vof
    use legacy_mesh_api, only: ncells
    use index_partitioning, only: gather_boundary
    !use material_interop, only: material_to_system
#ifdef EXTRA_VOF_DIAGNOSTICS
    use parallel_communication, only: global_minval, global_maxval
#endif

    type(matl_mesh_func), intent(inout) :: mmf

    integer :: i, m, p1, p2
    integer, pointer :: material_id(:)
    real(r8), allocatable :: vf(:), vofm(:)
    real(r8), pointer :: vfrac(:,:)
    class(base_mesh), pointer :: mesh
#ifdef EXTRA_VOF_DIAGNOSTICS
    character(len=90) :: string
#endif
    integer :: material_to_system(matl_model%nphase)
!    class(material), pointer :: matl

    ASSERT(mmf%num_reg() == 1)

    mesh  => mmf%mesh_ptr()
    vfrac => mmf%reg_vol_frac(1)

    if (.not.associated(vfrac)) return  ! single-material region; vfrac=1

    !! We assume an identity mapping from region cells to mesh cells.
    ASSERT(all(mmf%reg_cell(1) == (/(m,m=1,mesh%ncell)/)))

! All materials plus void exist in MMF
! non-void MMF material IDs are same as in MATL_MODEL
! MMF material ID 0 is void and corresponds to phase NMAT
! if void is present it is the first material in MMF

    material_id => mmf%reg_matl(1) ! array of material system IDs

    material_to_system(matl_model%nphase) = 0  ! in case last is void
    do m = 1, matl_model%nmatl_real
      call matl_model%get_matl_phase_index_range(m, p1, p2)
      material_to_system(p1:p2) = m
    end do

    allocate(vofm(ncells), vf(mesh%ncell))
    vfrac = 0.0_r8
    do m = 1, matl_model%nphase
      call gather_vof (m, vofm)
      vf(:mesh%ncell_onP) = vofm(:mesh%ncell_onP)
      call gather_boundary (mesh%cell_ip, vf)
      do i = size(material_id), 1, -1 ! find the destination dim
        if (material_id(i) == material_to_system(m)) exit
      end do
      INSIST(i > 0)
      vfrac(:,i) = vfrac(:,i) + vf
    end do
    deallocate(vofm, vf)

#ifdef EXTRA_VOF_DIAGNOSTICS
    INSIST(all(vfrac >= 0.0_r8 .and. vfrac <= 1.0_r8))
    write(string,'(3(a,f18.16))') 'UPDATE_MMF_FROM_MATL: vfrac sum interval: [', &
       global_minval(sum(vfrac,dim=2)), ', ', global_maxval(sum(vfrac,dim=2)) , ']'
    call TLS_info (string)
#endif

  end subroutine update_mmf_from_matl

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! UPDATE_MATL_FROM_MMF
 !!
 !! This routine updates the Truchas MATL structure with new VoF data derived
 !! from phase mixtures at the current solution state and material system
 !! volume fractions stored in the MMF structure.  Note that although the
 !! diffusion solver does not alter the material volume fractions (materials
 !! are conserved), the phase mixture of a multi-phase material will depend on
 !! the current value of the state variables and so the Truchas VoF data, which
 !! is really the volume fraction of phases, will need to be updated, but only
 !! when phase change is possible.
 !!
 !! It is very difficult to directly manipulate the compressed MATL structure
 !! to only update those VoFs involved in phase change.  So instead we use the
 !! provided UPDATE_MATL routine that essentially recreates MATL given the full
 !! VoF array generated here -- potentially a lot of needless effort.
 !!
 !! The routine assumes a single-region MMF having an identity mapping between
 !! region cells and mesh cells, and uses the mapping PHASE_TO_MATERIAL from
 !! the MATERIAL_INTEROP module.
 !!

  subroutine update_matl_from_mmf (mmf, state)

    use matl_utilities, only: update_matl, matl_get_cell_vof
    use legacy_mesh_api, only: ncells
#ifdef EXTRA_VOF_DIAGNOSTICS
    use parallel_communication, only: global_minval, global_maxval
#endif

    type(matl_mesh_func), intent(inout) :: mmf
    real(r8), intent(in) :: state(:,:)

    integer :: i, j, m, p1, p2
    integer,  pointer :: material_id(:)
    real(r8), pointer :: vfrac(:,:)
    real(r8), allocatable :: vofm(:), vof(:,:)
    class(base_mesh),  pointer :: mesh
#ifdef EXTRA_VOF_DIAGNOSTICS
    character(len=90) :: string
#endif

    mesh => mmf%mesh_ptr()

    ASSERT(mmf%num_reg() == 1)
    ASSERT(size(state,dim=1) == mesh%ncell_onP)

    material_id => mmf%reg_matl(1)

    mesh  => mmf%mesh_ptr()
    vfrac => mmf%reg_vol_frac(1)

    !! We assume an identity mapping from region cells to mesh cells.
    ASSERT(all(mmf%reg_cell(1) == (/(m,m=1,mesh%ncell)/)))

    allocate(vof(0:matl_model%nphase,ncells), vofm(ncells))
    vof = 0.0_r8

    !! Copy void volume fraction, if any, into VOF
    if (material_id(1) == 0) then
      ASSERT(matl_model%have_void)
      if (associated(vfrac)) then
        vofm(:mesh%ncell_onP) = vfrac(:mesh%ncell_onP,1)
        vof(matl_model%void_index,:) = vofm
      else  ! single-material region, just void!
        vof(matl_model%void_index,:) = 1.0_r8
      end if
    end if

    ! phase index == vof material index for non-void
    do i = 1, size(material_id) ! loop over MMF materials
      if (material_id(i) == 0) cycle  ! void -- no corresponding material
      call matl_model%get_matl_phase_index_range(material_id(i), p1, p2)
      select case (p2-p1+1)
      case (1)! Single-phase material system

        !! Copy its volume fraction into the Truchas material VoF array.
        m = p1
        if (associated(vfrac)) then
          vofm(:mesh%ncell_onP) = vfrac(:mesh%ncell_onP,i)
          vof(m,:) = vofm
        else  ! single-material region; volume fraction is 1
          vof(m,:) = 1.0_r8
        end if

      case default  ! Multi-phase material system

        !! Compute the state-dependent phase volume fractions for the material.
        if (associated(vfrac)) then
          do j = 1, mesh%ncell_onP
            if (vfrac(j,i) > 0.0_r8) then
              call matl_model%get_matl_phase_frac(material_id(i), state(j,1), vof(p1:p2,j))
              vof(p1:p2,j) = vfrac(j,i)*vof(p1:p2,j)
            else
              vof(p1:p2,j) = 0.0_r8
            end if
          end do
        else  ! single-material region; volume fraction is 1
          do j = 1, mesh%ncell_onP
            call matl_model%get_matl_phase_frac(material_id(i), state(j,1), vof(p1:p2,j))
          end do
        end if

      end select
    end do

    !! Copy the volume fraction data for gap elements from MATL into VOF.
    do j = mesh%ncell_onP+1, ncells
      call matl_get_cell_vof (j, vof(1:,j))
    end do

#ifdef EXTRA_VOF_DIAGNOSTICS
    INSIST(all(vof(1:,:) >= 0.0_r8 .and. vof(1:,:) <= 1.0_r8))
    write(string,'(3(a,f18.16))') 'UPDATE_MATL_FROM_MMF: vof sum interval: [', &
       global_minval(sum(vof(1:,:),dim=1)), ', ', global_maxval(sum(vof(1:,:),dim=1)) , ']'
    call TLS_info (string)
#endif

    !! Update the MATL structure using the uncompressed VoF data.
    call update_matl (vof)

  end subroutine update_matl_from_mmf

end module mesh_interop
