!!
!! MESH_INTEROP
!!
!! This module provides some data and procedures that facilitate the exchange
!! of mesh-based data between the Truchas mesh structure and the new mesh data
!! structure used by the diffusion solver.
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL GENERATE_MESH_MAPPINGS (MESH) generates the mappings that
!!    are needed to move cell based fields between the diffusion solver MESH
!!    and the old Truchas mesh structure.  The two meshes are geometrically
!!    equivalent but may differ in the cell numbering and parallel partitioning.
!!    In addition the Truchas mesh may include so-called "gap elements" that
!!    are not included in the diffusion solver mesh.  The following public
!!    module variables provide access to this information:
!!
!!    T_GAP_ELEMENTS(:)
!!        An array of the process-local Truchas mesh cell indices that are
!!        gap elements.
!!
!!    PCELL_T_TO_DS
!!        The parallel mapping object that describes the mapping from Truchas
!!        mesh cells to diffusion solver mesh cells.  It is used with the
!!        REARRANGE subroutine to pull-back a cell-based field defined on the
!!        diffusion solver mesh to the Truchas mesh, as in
!!            CALL REARRANGE (PCELL_T_TO_DS, DEST_T_FIELD, SRC_DS_FIELD).
!!        All elements of the destination field are defined except for those
!!        corresponding to the gap elements given by the T_GAP_ELEMENTS array.
!!
!!    PCELL_DS_TO_T
!!        The analogous parallel mapping object that describes the inverse
!!        mapping from diffusion solver mesh cells to Truchas mesh cells.
!!        It is used in exactly the same way to pull-back a cell-based field
!!        defined on the Truchas mesh to the diffusion solver mesh.  In this
!!        case source elements corresponding to gap elements are ignored and
!!        all elements of the destination field are defined.
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
  use parallel_permutations
  use base_mesh_class
  use material_mesh_function
  use truchas_logging_services, only: TLS_info
  implicit none
  private
  
  public :: generate_mesh_mappings, delete_mesh_mappings
  public :: mmf_init, update_mmf_from_matl, update_matl_from_mmf
  public :: void_is_present

  !! Permutation structures connecting the DS and T mesh cell labelings.
  type(par_perm), allocatable, public, save :: pcell_t_to_ds, pcell_ds_to_t
  integer, pointer, public, save :: t_gap_elements(:) => null()
  
contains

  subroutine generate_mesh_mappings (mesh)
    use mesh_module, only: unpermute_mesh_vector
    class(base_mesh), intent(in) :: mesh
    integer, pointer :: dummy(:) => null()
    allocate(pcell_t_to_ds, pcell_ds_to_t)
    call create_par_perm (unpermute_mesh_vector, mesh%xcell(:mesh%ncell_onP), &
                          pcell_t_to_ds, t_gap_elements, pcell_ds_to_t, dummy)
    INSIST(size(dummy) == 0)
    deallocate(dummy)
    INSIST(are_gap_elements(t_gap_elements))
  contains
    logical function are_gap_elements (list)
      use mesh_module, only: mesh, GAP_ELEMENT_1
      integer, intent(in) :: list(:)
      integer :: j
      are_gap_elements = .false.
      do j = 1, size(list)
        if (mesh(list(j))%cell_shape < GAP_ELEMENT_1) return
      end do
      are_gap_elements = .true.
    end function are_gap_elements
  end subroutine generate_mesh_mappings
  
  subroutine delete_mesh_mappings ()
    if (allocated(pcell_t_to_ds)) deallocate(pcell_t_to_ds)
    if (allocated(pcell_ds_to_t)) deallocate(pcell_ds_to_t)
    if (associated(t_gap_elements)) deallocate(t_gap_elements)
  end subroutine delete_mesh_mappings
  
  logical function void_is_present ()
    use material_interop, only: void_material_index
    void_is_present = (void_material_index > 0)
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
 !! These two lists (MMF_REG_CELL, MMF_REG_MATID) are normally needed to
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
 !! C. GENERATE_MESH_MAPPINGS has been called to establish the mesh mappings.
 !!

  subroutine mmf_init (mesh, mmf, stat, errmsg)

    use material_table
    use material_interop, only: void_material_index

    class(base_mesh), intent(in), target :: mesh
    type(mat_mf), intent(out) :: mmf
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j
    integer, allocatable :: all_materials(:)
    integer, pointer :: all_cell_sets(:)

    !! Single region consisting of the entire mesh allowing all materials.
    if (void_material_index > 0) then
      allocate(all_materials(0:mt_num_material()))
      all_materials(0) = 0  ! Special ID for void -- no corresponding material
    else
      allocate(all_materials(1:mt_num_material()))
    end if
    call mt_get_material_ids (all_materials(1:))
    all_cell_sets => mesh%cell_set_id
    call mmf_prep (mmf, mesh)
    call mmf_define_region (mmf, all_cell_sets, all_materials, stat, errmsg)
    deallocate(all_materials)
    if (stat /= 0) then
      errmsg = 'error defining MMF region: ' // trim(errmsg)
      return
    end if
    call mmf_done (mmf, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'error creating MMF: ' // trim(errmsg)
      return
    end if

    !! Assumed condition; see note.
    INSIST(all(mmf_reg_cell(mmf,1) == (/(j,j=1,mesh%ncell)/)))

    !! Finally initialize the MMF volume fraction data from the MATL VoF data.
    call update_mmf_from_matl (mmf)

    stat = 0
    errmsg = ''

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
 !! the MATERIAL_INTEROP module, and the mesh mapping PCELL_DS_TO_T.
 !!

  subroutine update_mmf_from_matl (mmf)

    use matl_module, only: gather_vof
    use parameter_module, only: nmat, ncells
    use index_partitioning, only: gather_boundary
    use material_interop, only: material_to_system
#ifdef EXTRA_VOF_DIAGNOSTICS
    use parallel_communication, only: global_minval, global_maxval
#endif

    type(mat_mf), intent(inout) :: mmf

    integer :: i, m
    integer, pointer :: material_id(:)
    real(r8), allocatable :: vf(:), vofm(:)
    real(r8), pointer :: vfrac(:,:)
    class(base_mesh), pointer :: mesh
#ifdef EXTRA_VOF_DIAGNOSTICS
    character(len=90) :: string
#endif

    ASSERT(mmf_num_reg(mmf) == 1)

    mesh  => mmf_mesh(mmf)
    vfrac => mmf_reg_vol_frac(mmf,1)

    if (.not.associated(vfrac)) return  ! single-material region; vfrac=1

    !! We assume an identity mapping from region cells to mesh cells.
    ASSERT(all(mmf_reg_cell(mmf,1) == (/(m,m=1,mesh%ncell)/)))

    material_id => mmf_reg_matid(mmf,1) ! array of material system IDs

    allocate(vofm(ncells), vf(mesh%ncell))
    vfrac = 0.0_r8
    do m = 1, nmat  ! loop over Truchas material numbers
      call gather_vof (m, vofm)
      call rearrange (pcell_ds_to_t, vf(:mesh%ncell_onP), vofm)
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

    use material_system
    use material_table
    use matl_utilities, only: update_matl, matl_get_cell_vof
    use parameter_module, only: nmat, ncells
    use material_interop, only: phase_to_material, void_material_index
#ifdef EXTRA_VOF_DIAGNOSTICS
    use parallel_communication, only: global_minval, global_maxval
#endif

    type(mat_mf), intent(inout) :: mmf
    real(r8), intent(in) :: state(:,:)

    integer :: i, j, k, m
    integer,  pointer :: material_id(:), phase_id(:)
    real(r8), pointer :: vfrac(:,:)
    real(r8), allocatable :: pfrac(:,:), vofm(:), vof(:,:)
    type(mat_system), pointer :: ms
    class(base_mesh),  pointer :: mesh
#ifdef EXTRA_VOF_DIAGNOSTICS
    character(len=90) :: string
#endif

    mesh => mmf_mesh(mmf)

    ASSERT(mmf_num_reg(mmf) == 1)
    ASSERT(size(state,dim=1) == mesh%ncell_onP)

    material_id => mmf_reg_matid(mmf,1)

    mesh  => mmf_mesh(mmf)
    vfrac => mmf_reg_vol_frac(mmf,1)

    !! We assume an identity mapping from region cells to mesh cells.
    ASSERT(all(mmf_reg_cell(mmf,1) == (/(m,m=1,mesh%ncell)/)))

    allocate(vof(0:nmat,ncells), vofm(ncells))
    vof = 0.0_r8
    
    !! Copy void volume fraction, if any, into VOF
    if (material_id(1) == 0) then
      ASSERT(void_material_index > 0)
      if (associated(vfrac)) then
        call rearrange (pcell_t_to_ds, vofm, vfrac(:mesh%ncell_onP,1))
        vof(void_material_index,:) = vofm
      else  ! single-material region, just void!
        vof(void_material_index,:) = 1.0_r8
      end if
    end if

    do i = 1, size(material_id) ! loop over MMF materials
      if (material_id(i) == 0) cycle  ! void -- no corresponding material
      ms => mt_get_material(material_id(i))
      ASSERT(associated(ms))
      call ms_get_phase_id (ms, phase_id)
      select case (size(phase_id))
      case (1)! Single-phase material system

        !! Copy its volume fraction into the Truchas material VoF array.
        m = phase_to_material(phase_id(1))
        if (associated(vfrac)) then
          call rearrange (pcell_t_to_ds, vofm, vfrac(:mesh%ncell_onP,i))
          vof(m,:) = vofm
        else  ! single-material region; volume fraction is 1
          vof(m,:) = 1.0_r8
        end if

      case default  ! Multi-phase material system

        !! Compute the state-dependent phase volume fractions for the material.
        allocate(pfrac(mesh%ncell_onP,size(phase_id)))
        if (associated(vfrac)) then
          do j = 1, mesh%ncell_onP
            if (vfrac(j,i) > 0.0_r8) then
              call ms_phase_mixture (ms, state(j,:), pfrac(j,:))
              pfrac(j,:) = vfrac(j,i)*pfrac(j,:)
            else
              pfrac(j,:) = 0.0_r8
            end if
          end do
        else  ! single-material region; volume fraction is 1
          do j = 1, mesh%ncell_onP
            call ms_phase_mixture (ms, state(j,:), pfrac(j,:))
          end do
        end if

        !! Copy into the Truchas material VoF array.
        do k = 1, size(phase_id)
          m = phase_to_material(phase_id(k))
          call rearrange (pcell_t_to_ds, vofm, pfrac(:,k))
          vof(m,:) = vofm
        end do
        deallocate(pfrac)

      end select
      deallocate(phase_id)
    end do
    deallocate(vofm)
    
    !! Copy the volume fraction data for gap elements from MATL into VOF.
    do i = 1, size(t_gap_elements)
      j = t_gap_elements(i)
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
    deallocate(vof)

  end subroutine update_matl_from_mmf
  
end module mesh_interop
