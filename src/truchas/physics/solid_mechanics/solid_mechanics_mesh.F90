!!
!! Solid Mechanics Mesh Module
!!      
!! This module holds a copy of the mesh pointer and defines subroutines
!! and functions related to the mesh required throughout the solid mechanics physics
!! package.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"
      
module solid_mechanics_mesh

  implicit none
  private

  public sm_mesh_init

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

contains

  subroutine sm_mesh_init
  end subroutine sm_mesh_init

end module solid_mechanics_mesh 
