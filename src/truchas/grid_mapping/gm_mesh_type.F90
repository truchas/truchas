!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gm_mesh_type
  !=======================================================================
  ! Purpose(s):
  !
  !    Provide GM_MESH derived type for GRID_MAPPING_MODULE
  !
  !    Public Interface(s):
  !
  !      type gm_mesh
  !
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov)
  !
  !=======================================================================
  implicit none
  private 

  integer, parameter :: dp = kind(1.0d0)

  type, public :: gm_mesh
     integer :: nnod = 0
     integer :: nelt = 0
     real(dp), allocatable :: pos_node(:,:)
     integer, allocatable :: node_elt(:,:)
     integer, allocatable :: block_elt(:)
  end type gm_mesh

end module gm_mesh_type

