!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MESH_GEN_DATA
  !=======================================================================
  ! Purpose(s):
  !
  !   The data required to determine how to generate or read a mesh.
  !
  !   **Note that the use of this module is being phased in.**
  !
  !   Public Interface:
  !
  !     * call set_Generated_Mesh()
  !
  !        Pass in a .TRUE. (if generating mesh internally) or 
  !                  .FALSE. (if not generating a mesh internally) flag
  !
  !     *  logical :: Generated_Mesh()
  !
  !        Returns .TRUE. if generating a mesh internally, .FALSE. otherwise.
  !
  !	
  ! Author(s): Robert Ferrell (ferrell.cpca.com)
  !
  !=======================================================================
  implicit none
  private

  public :: set_Generated_Mesh, Generated_Mesh

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  integer, save, public :: partitions_total

  ! Private parameters
  logical, SAVE :: Generated_Mesh_Flag

CONTAINS

  ! Set the generated_mesh flag
  subroutine set_generated_mesh (flag)
    logical :: flag
    generated_mesh_flag = flag
  end subroutine set_generated_mesh
  
  ! Return the Generated_Mesh flag
  logical function generated_mesh()
    generated_mesh = generated_mesh_flag
  end function generated_mesh

END MODULE MESH_GEN_DATA
  
  
  
