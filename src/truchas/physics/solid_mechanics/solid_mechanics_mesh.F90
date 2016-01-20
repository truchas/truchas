!!
!! Solid Mechanics Mesh Module
!!      
!! This module holds a copy of the mesh pointer and defines subroutines
!! and functions related to the mesh required throughout the solid mechanics physics
!! package.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PUBLIC Data 
!!
!!    SM_MESH
!!    A pointer to solid mechanics dist_mesh 
!!
!!    SM_MESH_TO_OLD_CELL, SM_OLD_TO_MESH_CELL
!!    SM_MESH_TO_OLD_NODE, SM_OLD_TO_MESH_NODE
!!    Cell and node based data structures that permute the old mesh to dist_mesh
!!
!!    SM_GAP_ELEMENT
!!    Integer pointer. I have no idea what this is. Indices map to gap elements
!!    old mesh? I was mimicking what I saw in the heat transfer module.
!!
!! PUBLIC Parameters
!!
!!    These are duplicate parameters that originally came from the
!!    mesh module. Preserved here to maintain the package capability to
!!    handle meshes that are not hex meshes. 
!!    SM_CELL_TET
!!    SM_CELL_PYRAMID
!!    SM_CELL_PRISM
!!    SM_CELL_HEX
!!    SM_GAP_ELEMENT_1      
!!    SM_GAP_ELEMENT_3      
!!    SM_GAP_ELEMENT_5
!!
!!    These are duplicate parameters from the parameter module. Solid mechanics
!!    modules should these parameters instead of the PARAMETER_MODULE parameters
!!    These should be removed one the change to DIST_MESH is complete.
!!    NDIM
!!    NFC
!!    NVC
!!    NVF
!!    NEC
!!    NROT
!!    NCOMPS
!!
!! PUBLIC Subroutine/Functions
!!
!!    SM_MESH_INIT
!!    Initialize (copy) SM_MESH pointer and other public data 
!!
!!    SM_CELL_SHAPE
!!    A function that return cell shape parameter. Only checks the number of
!!    nodes of the first cell element. THIS WILL FAIL with mixed element meshes
!!
!!    SM_FACE_CENTROID
!!    Compute the physical face centroid. Based on FACE_CENTROID_PHYSICAL in the cell geometry
!!    module.
!!
!!    SM_CELL_CENTROID
!!    Compute the physical cell centroid. Based on CELL_CENTROID_PHYSICAL in the cell geometry
!!    module.
!!
!!    SM_CELL_CENTROID_LOGICAL
!!    Based on CELL_CENTROID_LOGICAL in the cell geometry module.

#include "f90_assert.fpp"
      

module solid_mechanics_mesh

  use unstr_mesh_type
  use parallel_permutations, only: par_perm

  implicit none
  private

  !! Solid mechanics DIST_MESH Pointer
  type(unstr_mesh), pointer, public :: sm_mesh => null()

  !! Permutation structures to transfer data to and from old mesh
  type(par_perm), target, save, public :: sm_mesh_to_old_cell, sm_old_to_mesh_cell
  type(par_perm), target, save, public :: sm_mesh_to_old_node, sm_old_to_mesh_node

  !! I don't know if this should be a module variable or local in sm_mesh_init
  integer, pointer, public, save :: sm_gap_elements(:) => null()

  !!
  public sm_mesh_init!, sm_cell_shape

!NNC: not yet used -- disabling
!NNC  !! Imported Cell geometry subroutines
!NNC  public :: sm_face_centroid, sm_cell_centroid, sm_face_centroid_logical

  !! From the old mesh module a data block defining the edge ordering of a HEX
  !! element.
  integer, public :: Cell_Edge(2,12)
  data Cell_Edge/1,2, 2,3, 3,4, 4,1, 2,6, 3,7, 4,8, 1,5, 5,6, 6,7, 7,8, 8,5/

  
  !! Parameter definitions originally from the old mesh module
  integer, parameter, public :: SM_CELL_TET     = 1
  integer, parameter, public :: SM_CELL_PYRAMID = 2
  integer, parameter, public :: SM_CELL_PRISM   = 3
  integer, parameter, public :: SM_CELL_HEX     = 4
  !! define another cell shape for "zero volume" gap elements
  !! The gap element number corresponds to the first non-degenerate cell face
  integer, parameter, public :: SM_GAP_ELEMENT_1  = 5
  integer, parameter, public :: SM_GAP_ELEMENT_3  = 6
  integer, parameter, public :: SM_GAP_ELEMENT_5  = 7

  !! These are parameters orignally from the parameter module
  !! Eventually these will be removed and replaced with 
  !! Number of physical dimensions
  integer, parameter, public :: ndim = 3

  !! Number of faces per cell
  integer, parameter, public :: nfc = 2*ndim

  !! Number of vertices per cell
  integer, parameter, public :: nvc = 2**ndim

  !! Number of vertices per face
  integer, parameter, public :: nvf = 2**(ndim - 1)

  !! Number of edges per cell
  integer, parameter, public :: nec = ndim*2**(ndim - 1)

  !! Number of rotation axes
  integer, parameter, public :: nrot = 2*ndim - 3

  !! Number of stress/strain components
  integer, parameter, public :: ncomps = (ndim*(ndim + 1))/2
  
  !! Namelist mesh name, used in sm_mesh_init
  character(4) :: sm_mesh_name = 'MAIN'


contains

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 

  ! Purpose: Initialize the solid mechanics mesh pointer and the permutation
  !          structures.
  subroutine sm_mesh_init
    use mesh_manager, only: unstr_mesh_ptr
    use legacy_mesh_api, only: unpermute_mesh_vector, unpermute_vertex_vector
    use parallel_permutations, only: create_par_perm

    integer, pointer :: dummy(:) => null()
    logical :: found

    ! Initialize the mesh pointer
    sm_mesh => unstr_mesh_ptr(sm_mesh_name)
    INSIST(associated(sm_mesh))

    ! Initialize the permutation structures -- I've copied
    ! the generate_mesh_mappings routine from  mesh_interop.F90
    call create_par_perm(unpermute_mesh_vector, &
                         sm_mesh%xcell(:sm_mesh%ncell_onp), &
                         sm_old_to_mesh_cell, sm_gap_elements, &
                         sm_mesh_to_old_cell, dummy)
    INSIST(size(dummy) == 0)
    deallocate(dummy)
    INSIST(are_gap_elements(sm_gap_elements))

    call create_par_perm(unpermute_vertex_vector, &
                         sm_mesh%xnode(:sm_mesh%nnode_onp), &
                         sm_old_to_mesh_node, sm_gap_elements, &
                         sm_mesh_to_old_node,dummy)
    INSIST(size(dummy) == 0)
    deallocate(dummy)
    INSIST(are_gap_elements(sm_gap_elements))

    contains
      logical function are_gap_elements (list)
      use legacy_mesh_api, only: mesh, GAP_ELEMENT_1
      integer, intent(in) :: list(:)
      integer :: j
      are_gap_elements = .false.
      do j = 1, size(list)
        if (mesh(list(j))%cell_shape < GAP_ELEMENT_1) return
      end do
      are_gap_elements = .true.
    end function are_gap_elements
  end subroutine sm_mesh_init

End module solid_mechanics_mesh 
