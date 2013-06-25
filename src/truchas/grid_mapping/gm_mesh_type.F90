Module GM_MESH_TYPE
  !=======================================================================
  ! Purpose(s):
  !
  !    Provide GM_MESH derived type for GRID_MAPPING_MODULE
  !
  !    Public Interface(s):
  !
  !      destroy_gm_mesh
  !      type gm_mesh
  !
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov)
  !
  !=======================================================================
  implicit none

  private 

  PUBLIC :: destroy_gm_mesh

  integer, parameter :: dp = KIND(1.0d0)

  type, public :: gm_mesh
     integer :: nnod
     integer :: nelt
     real(dp), dimension(:,:), pointer :: pos_node => null()
     integer, dimension(:,:), pointer :: node_elt => null()
     integer, dimension(:), pointer :: block_elt => null()
  end type gm_mesh

contains

  subroutine destroy_gm_mesh(mesh)
    type(gm_mesh), intent(inout) :: mesh
    mesh%nnod=0
    mesh%nelt=0
    if (associated(mesh%pos_node)) deallocate(mesh%pos_node)
    if (associated(mesh%node_elt)) deallocate(mesh%node_elt)
    if (associated(mesh%block_elt)) deallocate(mesh%block_elt)
    return
  end subroutine destroy_gm_mesh

end module GM_MESH_TYPE

