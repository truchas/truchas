!!
!! DEGEN_HEX_TOPOLOGY
!!
!! This module provides mappings between the labeling of cell nodes and sides
!! used by the mesh data structure and their labeling as degenerate hexahedra,
!! which is still used for output and restart input.
!!
!! These parameters were originally part of the legacy mesh API. Here "new"
!! refers to the labeling currently used in Truchas, and "old" to the labeling
!! as degenerate hexahedra. The labeling of hex cell nodes is the same in both
!! schemes.
!!

#include "f90_assert.fpp"

module degen_hex_topology

  implicit none
  private

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

end module degen_hex_topology
