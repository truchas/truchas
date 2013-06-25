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
  use kind_module, only : int_kind, &
                          log_kind
  use constants_module, only: IPRESET

  ! Private Module
  private

  ! Public Subroutines and functions
  public :: set_Generated_Mesh, &
            Generated_Mesh

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  integer (INT_KIND), save, public :: partitions_total
  integer (INT_KIND), save, public :: partitions_per_process

  integer (INT_KIND), parameter, public :: partitions_total_DEFAULT       = IPRESET
  integer (INT_KIND), parameter, public :: partitions_per_process_DEFAULT = IPRESET

  ! Private parameters
  logical, SAVE :: Generated_Mesh_Flag

CONTAINS

  SUBROUTINE set_Generated_Mesh(Flag)
    !=======================================================================
    ! Purpose(s):
    !
    !   Set the Generated_Mesh flag
    !
    !   INPUT:
    !        logical :: Flag
    !
    !   OUTPUT:
    !
    !=======================================================================
    implicit none
    
    ! Subroutine arguments
    logical(log_kind) :: Flag

    Generated_Mesh_Flag = Flag
    return
  end SUBROUTINE set_Generated_Mesh
  
  
  FUNCTION Generated_Mesh()
    !=======================================================================
    ! Purpose(s):
    !
    !   Return the Generated_Mesh flag
    !
    !   INPUT:
    !
    !   OUTPUT:
    !         .TRUE. if the mesh is generated internally, otherwise .FALSE.
    !
    !=======================================================================
    implicit none
    logical (log_kind) :: Generated_Mesh
    
    Generated_Mesh = Generated_Mesh_Flag

    return
  end FUNCTION Generated_Mesh

END MODULE MESH_GEN_DATA
  
  
  
