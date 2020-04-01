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
  implicit none
  private

  public :: init_common_impl

  type(unstr_mesh), pointer, public :: new_mesh => null()
  integer, public, protected :: ncells, nnodes, ncells_tot, nnodes_tot, ncells_real
  integer, public, allocatable, protected :: unpermute_mesh_vector(:), unpermute_vertex_vector(:)

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

  subroutine init_common_impl

    use mesh_manager, only: unstr_mesh_ptr
    use parallel_communication, only: global_sum

    integer :: j
    integer, allocatable :: gap_links(:)

    new_mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(new_mesh))

    nnodes = new_mesh%nnode_onP
    nnodes_tot = global_sum(nnodes)
    unpermute_vertex_vector = new_mesh%xnode(:nnodes)

    !! Mesh links generated from elements in the Exodus mesh are presented as
    !! regular cells by the legacy API.  These are appended to the list of mesh
    !! cells locally.
    gap_links = pack(array=[(j, j=1, new_mesh%nlink_onP)], &
                     mask=(new_mesh%link_cell_id(:new_mesh%nlink_onP)>0))
    ncells_real = new_mesh%ncell_onP
    ncells = new_mesh%ncell_onP + size(gap_links)
    ncells_tot = global_sum(ncells)
    allocate(unpermute_mesh_vector(ncells))
    unpermute_mesh_vector(:new_mesh%ncell_onP) = new_mesh%xcell(:new_mesh%ncell_onP)
    unpermute_mesh_vector(new_mesh%ncell_onP+1:) = new_mesh%link_cell_id(gap_links)

  end subroutine init_common_impl

end module common_impl
