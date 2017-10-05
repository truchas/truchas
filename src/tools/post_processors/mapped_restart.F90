!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module mapped_restart

  use grid_mapping_module, only: grid_int_vols
  use kinds, only: r8
  implicit none
  private
  
  public :: write_mapped_restart
  
  type :: mesh_map
    private
    type(grid_int_vols) :: map_data
  end type mesh_map
  
  private :: mesh_map_init, mesh_map_delete
  private :: mesh_map_map_field
  
  integer, parameter, public :: MAP_FORM_DEFAULT  = 0
  integer, parameter, public :: MAP_FORM_CONSTANT = 1
  integer, parameter, public :: MAP_FORM_INTEGRAL = 2
  
  type :: base_mesh
    integer :: ncell, nnode
    real(r8), pointer :: coord(:,:)   => null() ! node coordinates
    integer,  pointer :: connect(:,:) => null() ! cell connectivity
    integer,  pointer :: blockid(:)   => null() ! cell block IDs
    logical :: is_tet_mesh = .false.
  end type base_mesh
  
  private :: base_mesh_init, base_mesh_delete
  
  interface base_mesh_init
    module procedure base_mesh_init_danu, base_mesh_init_exo
  end interface

contains

  subroutine write_mapped_restart (ofile, seq_num, exomesh, unit)
  
    use h5out_type
    use exodus_mesh_type
    use danu_module
    use,intrinsic :: iso_c_binding, only: c_ptr
    
    type(h5out), intent(inout) :: ofile
    integer, intent(in) :: seq_num
    type(exodus_mesh), intent(in) :: exomesh
    integer, intent(in) :: unit
    
    integer :: stat
    type(c_ptr) :: seq_id
    type(base_mesh) :: src_mesh, dest_mesh
    type(mesh_map) :: map

    call base_mesh_init (src_mesh, ofile)
    call base_mesh_init (dest_mesh, exomesh)
    call mesh_map_init (map, src_mesh, dest_mesh)
    
    call sequence_get_id (ofile%sim_id, seq_num, seq_id, stat)
    INSIST(stat == DANU_SUCCESS)

    call write_header_segment
    call write_mesh_segment
    call write_core_data_segment
    !if (has_feature(seq_id,'solid_mechanics')) call write_solid_mech_data_segment
    if (has_feature(seq_id,'species'))         call write_species_data_segment
    !if (has_feature(seq_id,'joule_heat'))      call write_joule_heat_data_segment

  contains
  
    subroutine write_header_segment

      real(r8) :: t, dt
      integer  :: n, stat, cycle_num
      character(8), parameter :: fileid = 'TRF-3'
      character(32), allocatable :: feature_list(:)

      !! Item 1 -- file format magic number
      write(unit) fileid

      !! Item 2 -- file features
      !! The only one retained for a mapped restart is species.
      if (has_feature(seq_id, 'species')) then
        allocate(feature_list(1))
        feature_list(1) = 'species'
      else
        allocate(feature_list(0))
      end if
      write(unit) size(feature_list)
      do n = 1, size(feature_list)
        write(unit) feature_list(n)
      end do
      deallocate(feature_list)

      !! Item 3 -- simulation specification (free use -- not using any now)
      write(unit) 0

      !! Item 4 -- global data
      call attribute_read (seq_id, 'time', t, stat)
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (seq_id, 'time step', dt, stat)
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (seq_id, 'cycle', cycle_num, stat)
      INSIST(stat == DANU_SUCCESS)

      write(unit) t
      write(unit) dt
      write(unit) cycle_num
      write(unit) dest_mesh%ncell
      write(unit) dest_mesh%nnode

    end subroutine write_header_segment

    subroutine write_mesh_segment

      integer :: k

      do k = 1, size(dest_mesh%connect,1)
        write(unit) dest_mesh%connect(k,:)
      end do

      write(unit) 1
      write(unit) dest_mesh%blockid

      write(unit) dest_mesh%coord(1,:)
      write(unit) dest_mesh%coord(2,:)
      write(unit) dest_mesh%coord(3,:)

    end subroutine write_mesh_segment
    
    subroutine write_core_data_segment

      integer :: stat, k, nmat, vof_shape(2)
      real(r8) :: src_array(ofile%ncells), dest_array(dest_mesh%ncell)
      real(r8), allocatable :: array2(:,:)

      !! RHO -- cell densities
      call simulation_data_read (seq_id, 'Z_RHO', src_array, stat)
      INSIST(stat == DANU_SUCCESS)
      call mesh_map_map_field (map, MAP_FORM_CONSTANT, 0.0_r8, src_array, dest_array)
      write(unit) dest_array

      !! TEMP -- cell temperatures
      call simulation_data_read (seq_id, 'Z_TEMP', src_array, stat)
      INSIST(stat == DANU_SUCCESS)
      call mesh_map_map_field (map, MAP_FORM_CONSTANT, 0.0_r8, src_array, dest_array)
      write(unit) dest_array

      !! ENTHALPY -- cell enthalpies
      call simulation_data_read (seq_id, 'Z_ENTHALPY', src_array, stat)
      INSIST(stat == DANU_SUCCESS)
      call mesh_map_map_field (map, MAP_FORM_CONSTANT, 0.0_r8, src_array, dest_array)
      write(unit) dest_array

      !! Dummy flow data.
      dest_array = 0.0_r8
      do k = 1, 10
        write(unit) dest_array
      end do

      !! VF_* -- phase volume fractions
      !! N.B.  The volume fractions are supposed to sum to one.  If we use constant
      !! preserving mapping, then under linearity of the mapping(?), the sum of the
      !! mapped volume fractions should also sum to one.
      if (simulation_data_exists(seq_id, 'VOF')) then
        call simulation_data_dimensions (seq_id, 'VOF', vof_shape, stat)
        INSIST(stat == DANU_SUCCESS)
        nmat = vof_shape(1)
        allocate(array2(nmat,ofile%ncells))
        call simulation_data_read (seq_id, 'VOF', array2, stat)
        INSIST(stat == DANU_SUCCESS)
        write(unit) nmat
        do k = 1, nmat
          call mesh_map_map_field (map, MAP_FORM_CONSTANT, 0.0_r8, array2(k,:), dest_array)
          write(unit) dest_array
        end do
        deallocate(array2)
      else  ! single phase problem -- volume fraction is implicit
        dest_array = 1.0_r8
        write(unit) 1
        write(unit) dest_array
      end if

    end subroutine write_core_data_segment

    subroutine write_species_data_segment

      integer :: stat, n, num_species
      character(8) :: name
      real(r8) :: src_array(src_mesh%ncell), dest_array(dest_mesh%ncell)

      call attribute_read (ofile%sim_id, 'NUM_SPECIES', num_species, stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(num_species > 0)
      write(unit) num_species

      do n = 1, num_species
        write(name,'(a,i0)') 'phi', n
        call simulation_data_read (seq_id, name, src_array, stat)
        INSIST(stat == DANU_SUCCESS)
        call mesh_map_map_field (map, MAP_FORM_CONSTANT, 0.0_r8, src_array, dest_array)
        write(unit) dest_array
      end do

    end subroutine write_species_data_segment

  end subroutine write_mapped_restart

  subroutine mesh_map_init (this, src, dest)
    use grid_mapping_module, only: gm_mesh, compute_int_volumes
    type(mesh_map), intent(out) :: this
    type(base_mesh), intent(in) :: src, dest
    type(gm_mesh) :: gm_src, gm_dest
    !! Wire the mapping source mesh.
    gm_src%nnod = src%nnode
    gm_src%nelt = src%ncell
    gm_src%pos_node  = src%coord
    if (src%is_tet_mesh) then
      gm_src%node_elt  = src%connect(2:5,:)
    else
      gm_src%node_elt  = src%connect
    end if
    gm_src%block_elt = src%blockid
    !! Wire the mapping destination mesh.
    gm_dest%nnod = dest%nnode
    gm_dest%nelt = dest%ncell
    gm_dest%pos_node  = dest%coord
    if (dest%is_tet_mesh) then
      gm_dest%node_elt = dest%connect(2:5,:)
    else
      gm_dest%node_elt  = dest%connect
    end if
    gm_dest%block_elt = dest%blockid
    !! Compute the mesh mapping data.
    call compute_int_volumes (gm_src, gm_dest, this%map_data)
  end subroutine mesh_map_init
  
  subroutine mesh_map_delete (this)
    use grid_mapping_module, only: destroy_grid_int_vols
    type(mesh_map), intent(inout) :: this
    call destroy_grid_int_vols (this%map_data)
  end subroutine mesh_map_delete
  
  subroutine mesh_map_map_field (this, form, defval, src, dest)
  
    use grid_mapping_module, only: map_cell_field
  
    type(mesh_map), intent(in) :: this
    integer, intent(in) :: form
    real(r8), intent(in) :: defval
    real(r8), intent(in) :: src(:)
    real(r8), intent(out) :: dest(:)
    
    select case (form)
    case (MAP_FORM_CONSTANT)  ! THIS PRESERVES A CONSTANT-VALUED FIELD
      call map_cell_field (src, dest, this%map_data, strict=.true., defval=defval, preserve_constants=.true.)
          
    case (MAP_FORM_INTEGRAL)  ! THIS PRESERVES THE INTEGRAL OF THE FIELD
      call map_cell_field (src, dest, this%map_data, strict=.true., defval=defval, exactly_conservative=.true.)
          
    case (MAP_FORM_DEFAULT)   ! I DON'T ACTUALLY KNOW WHAT THIS DOES
      call map_cell_field (src, dest, this%map_data, strict=.true., defval=defval)
          
    case default
      INSIST(.false.)
    end select

  end subroutine mesh_map_map_field

  !! Initialize a BASE_MESH from a Danu HDF output file

  subroutine base_mesh_init_danu (this, ofile)
    use h5out_type
    type(base_mesh), intent(out) :: this
    type(h5out), intent(in) :: ofile
    this%nnode = h5out_num_nodes(ofile)
    this%ncell = h5out_num_cells(ofile)
    allocate(this%coord(3,this%nnode))
    call h5out_get_coordinates (ofile, this%coord)
    allocate(this%connect(8,this%ncell))
    call h5out_get_connectivity (ofile, this%connect)
    allocate(this%blockid(this%ncell))
    call h5out_get_block_ids (ofile, this%blockid)
    this%is_tet_mesh = all(this%connect(1,:) == this%connect(2,:))
  end subroutine base_mesh_init_danu
  
  !! Initialize a BASE_MESH from an ExodusII mesh

  subroutine base_mesh_init_exo (this, mesh)

    use exodus_mesh_type
    use,intrinsic :: iso_fortran_env, only: error_unit

    type(base_mesh),   intent(out) :: this
    type(exodus_mesh), intent(in)  :: mesh

    !! Truchas degenerate hex node numbering to ExodusII element node numbering.
    integer, parameter ::   TET_NODE_MAP(8) = (/1,1,2,3,4,4,4,4/)
    integer, parameter :: WEDGE_NODE_MAP(8) = (/1,4,5,2,3,6,6,3/)

    !! ExodusII element types identified by their number of nodes.
    integer, parameter :: TET=4, WEDGE=6, HEX=8

    integer :: n, i, j, nodes_per_elem
    logical :: non_hex_mesh

    this%nnode = mesh%num_node
    this%ncell = mesh%num_elem
    
    if (mesh%num_dim /= 3) then
      write(error_unit,'(a,i0)') 'Error: target mesh is not 3D: num_dim=', mesh%num_dim
      stop
    end if
    allocate(this%coord(3,this%nnode))
    this%coord = mesh%coord

    allocate(this%connect(8,this%ncell), this%blockid(this%ncell))

    !! Translate Exodus element connectivity to Truchas convention.
    n = 0
    non_hex_mesh = .false.
    this%is_tet_mesh = .true.
    do i = 1, mesh%num_eblk
      nodes_per_elem = size(mesh%eblk(i)%connect,dim=1)
      !! Translate Exodus element connectivity to Truchas convention.  See the ExodusII
      !! section in MESH_READ from MESH_INPUT_MODULE and the later section that converts
      !! non-hex elements into degenerate hexes as the reference for what must be done here.
      select case (nodes_per_elem)
      case (HEX)
        this%is_tet_mesh = .false.
        do j = 1, mesh%eblk(i)%num_elem
          n = n + 1
          this%blockid(n) = mesh%eblk(i)%ID
          this%connect(:,n) = mesh%eblk(i)%connect(:,j)
        end do
      case (TET)
        non_hex_mesh = .true.
        do j = 1, mesh%eblk(i)%num_elem
          n = n + 1
          this%blockid(n) = mesh%eblk(i)%ID
          this%connect(:,n) = mesh%eblk(i)%connect(TET_NODE_MAP,j)
        end do
      case (WEDGE)
        this%is_tet_mesh = .false.
        non_hex_mesh = .true.
        do j = 1, mesh%eblk(i)%num_elem
          n = n + 1
          this%blockid(n) = mesh%eblk(i)%ID
          this%connect(:,n) = mesh%eblk(i)%connect(WEDGE_NODE_MAP,j)
        end do
      case default
        write(error_unit,'(a,i0,a)') 'Error: unknown element type in target mesh: ', nodes_per_elem, '-node element'
        stop
      end select
    end do
    INSIST(n == this%ncell)

    if (non_hex_mesh .and. .not.this%is_tet_mesh) then
      write(error_unit,'(a)') 'Warning: target mesh contains degenerate hex elements'
    end if

  end subroutine base_mesh_init_exo
  
  subroutine base_mesh_delete (this)
    type(base_mesh), intent(inout) :: this
    if (associated(this%coord)) deallocate(this%coord)
    if (associated(this%connect)) deallocate(this%connect)
    if (associated(this%blockid)) deallocate(this%blockid)
  end subroutine base_mesh_delete

end module mapped_restart
