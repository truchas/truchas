!!
!! EXT_EXODUS_MESH_TYPE
!!
!! This module defines an extension of the EXODUS_MESH type that adds data
!! components to store mesh interface link info which describes a connection
!! of element sides across an interface.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ext_exodus_mesh_type

  use exodus_mesh_type
  implicit none
  private
  
  public :: exodus_mesh, elem_blk, node_set, side_set ! re-export

  type, extends(exodus_mesh), public :: ext_exodus_mesh
    integer :: nlink=0, nlblock=0
    integer, allocatable :: xlnode(:), lnode(:)
    integer, allocatable :: link_block(:)
    integer, allocatable :: link_block_id(:)
    !! Additional link data aiding transition from old mesh (TEMPORARY)
    integer, allocatable :: link_cell_id(:)  ! index of cell that became link (or 0)
    integer, allocatable :: parent_node(:)   ! index of the parent node (or 0)
  contains
    procedure :: set_no_links
    procedure :: delete_links
  end type ext_exodus_mesh

contains

  !! Sets the number of links to 0 and allocates the link arrays accordingly.
  subroutine set_no_links (this)
    class(ext_exodus_mesh), intent(inout) :: this
    this%nlink = 0
    this%nlblock = 0
    this%xlnode = [1]
    this%lnode = [integer::]  ! 0-sized array
    this%link_block = [integer::]  ! 0-sized array
    this%link_block_id = [integer::]  ! 0-sized array
    this%link_cell_id = [integer::]   ! 0-sized array
  end subroutine set_no_links

  !! Sets the number of links to 0 and deallocates the link arrays if allocated.
  subroutine delete_links (this)
    class(ext_exodus_mesh), intent(inout) :: this
    this%nlink = 0
    this%nlblock = 0
    if (allocated(this%xlnode)) deallocate(this%xlnode)
    if (allocated(this%lnode)) deallocate(this%lnode)
    if (allocated(this%link_block)) deallocate(this%link_block)
    if (allocated(this%link_block_id)) deallocate(this%link_block_id)
    if (allocated(this%link_cell_id)) deallocate(this%link_cell_id)
    if (allocated(this%parent_node)) deallocate(this%parent_node)
  end subroutine delete_links

end module ext_exodus_mesh_type
