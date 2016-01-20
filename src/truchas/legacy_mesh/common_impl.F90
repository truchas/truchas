!!
!! COMMON_IMPL
!!
!! Provides the core data describing the relationship between the legacy mesh
!! and new mesh object.  This includes various mappings connecting the different
!! partitioning of cells and nodes, different conventions for numbering the
!! vertices and sides of a cell, and relation of gap cells to new mesh links.
!!

#include "f90_assert.fpp"

module common_impl

  use unstr_mesh_type, only: unstr_mesh
  use parallel_permutations, only: par_perm
  use gap_node_map_type
  implicit none
  private

  public :: init_common_impl

  integer, public, allocatable, protected :: unpermute_mesh_vector(:), unpermute_vertex_vector(:)

  integer, public, protected :: ncells, nnodes
  integer, public, protected :: ncells_tot, nnodes_tot

  type(unstr_mesh), pointer, public :: new_mesh => null()
  type(par_perm),   public :: pcell_new_to_old, pcell_old_to_new
  type(par_perm),   public :: pnode_new_to_old, pnode_old_to_new
  integer, pointer, public :: gap_cells(:) => null()  ! exist in old but not new
  integer, pointer, public :: gap_nodes(:) => null()  ! exist in new but not old

  logical, allocatable, public :: gap_link_mask(:)
  type(par_perm), public :: pgap_old_to_new, pgap_new_to_old
  type(gap_node_map), public :: gnmap

  !! Mappings from old cell sides to new. These imagine new cell side data
  !! stored in the leading part of an array of size 6 with dummy data in the
  !! unused part. These mappings point to that dummy data for degenerate sides.
  integer, parameter, public :: OLD_TET_SIDE_MAP(6) = [2,5,3,1,4,6]  ! old sides 2 and 6 are degenerate
  integer, parameter, public :: OLD_PYR_SIDE_MAP(6) = [3,1,4,2,5,6]  ! old side 6 is degenerate
  integer, parameter, public :: OLD_PRI_SIDE_MAP(6) = [2,3,4,5,1,6]  ! old side 6 is degenerate
  integer, parameter, public :: OLD_HEX_SIDE_MAP(6) = [3,1,4,2,5,6]

  !! And the inverse mappings from new cell sides to old.
  integer, parameter, public :: NEW_TET_SIDE_MAP(6) = [4,1,3,5,2,6] ! dummy sides 5, 6 point to degenerate
  integer, parameter, public :: NEW_PYR_SIDE_MAP(6) = [2,4,1,3,5,6] ! dummy side 6 points to degenerate
  integer, parameter, public :: NEW_PRI_SIDE_MAP(6) = [5,1,2,3,4,6] ! dummy side 6 points to degenerate
  integer, parameter, public :: NEW_HEX_SIDE_MAP(6) = [2,4,1,3,5,6]

  !! Mappings from old cell nodes to new.
  integer, parameter, public :: OLD_TET_NODE_MAP(8) = [1,1,2,3,4,4,4,4]
  integer, parameter, public :: OLD_PYR_NODE_MAP(8) = [1,2,3,4,5,5,5,5]
  integer, parameter, public :: OLD_PRI_NODE_MAP(8) = [1,4,5,2,3,6,6,3]

  !! Mappings from new cell nodes to old.
  integer, parameter, public :: NEW_TET_NODE_MAP(4) = [2,3,4,5]
  integer, parameter, public :: NEW_PYR_NODE_MAP(5) = [1,2,3,4,5]
  integer, parameter, public :: NEW_PRI_NODE_MAP(6) = [1,4,5,2,3,6]

contains

  subroutine init_common_impl (old_cell_perm, old_node_perm)

    use mesh_manager, only: unstr_mesh_ptr
    use parallel_permutations, only: create_par_perm
    use parallel_communication, only: global_sum

    integer, intent(in) :: old_cell_perm(:), old_node_perm(:)
    integer, pointer :: dummy(:)

    new_mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(new_mesh))

    unpermute_mesh_vector = old_cell_perm
    unpermute_vertex_vector = old_node_perm

    ncells = size(unpermute_mesh_vector)
    nnodes = size(unpermute_vertex_vector)

    ncells_tot = global_sum(ncells)
    nnodes_tot = global_sum(nnodes)

    !! Cell mappings between the old and new meshes.
    call create_par_perm (unpermute_mesh_vector, new_mesh%xcell(:new_mesh%ncell_onP), &
                          pcell_old_to_new, gap_cells, pcell_new_to_old, dummy)
    INSIST(size(dummy) == 0)  ! every new cell is an old cell
    deallocate(dummy)

    !! Node mappings between the old and new meshes.
    call create_par_perm (unpermute_vertex_vector, new_mesh%xnode(:new_mesh%nnode_onP), &
                          pnode_old_to_new, dummy, pnode_new_to_old, gap_nodes)
    INSIST(size(dummy) == 0)  ! every old node is a new node
    deallocate(dummy)

    call init_gap_maps

  end subroutine init_common_impl

  !! Initializes the various maps involving gap cells and nodes.
  !! Input is GAP_CELLS

  subroutine init_gap_maps

    use parallel_permutations, only: create_par_perm

    integer, allocatable :: gap_cell_id(:), gap_link_id(:)
    integer, pointer :: dummy1(:), dummy2(:)
    integer :: unit, bndry_size

    !! On-process link mask identifying links coming from gap cells.
    gap_link_mask = (new_mesh%link_cell_id(:new_mesh%nlink_onP) > 0)

    !! Permutations between old mesh gap cells and new mesh links.
    gap_cell_id = unpermute_mesh_vector(gap_cells)
    gap_link_id = pack(new_mesh%link_cell_id(:new_mesh%nlink_onP), mask=gap_link_mask)
    call create_par_perm (gap_cell_id, gap_link_id, pgap_old_to_new, dummy1, &
                                                    pgap_new_to_old, dummy2)
    INSIST(size(dummy1) == 0)
    deallocate(dummy1)
    INSIST(size(dummy2) == 0)
    deallocate(dummy2)

    !! Initialize the mapping between gap nodes and their parents.
    call gnmap%init (new_mesh)

  end subroutine init_gap_maps

end module common_impl
