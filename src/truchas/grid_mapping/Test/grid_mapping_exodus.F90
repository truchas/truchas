module grid_mapping_exodus
  !=======================================================================
  ! Purpose(s):
  !
  !    Convert an EXODUS mesh object to a GM_MESH object.
  !
  !    Public Interface(s):
  !
  !      copy_exodus_mesh_to_gm_mesh
  !
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov)
  !
  !=======================================================================

  use gm_mesh_type, only : gm_mesh
  implicit none

  private

  public :: copy_exodus_mesh_to_gm_mesh

contains

  subroutine copy_exodus_mesh_to_gm_mesh(e_mesh,g_mesh)

    use exodus_mesh_type

    implicit none 

    type(exodus_mesh), intent(in) :: e_mesh
    type(gm_mesh), intent(out) :: g_mesh
    integer :: i,j,ielt,nodes_per_elt
    logical :: tetmesh

    g_mesh%nnod=e_mesh%num_node
    g_mesh%nelt=e_mesh%num_elem
    allocate(g_mesh%pos_node(3,g_mesh%nnod))
    g_mesh%pos_node = e_mesh%coord
    ! If all element blocks have 4 nodes per element, we make
    ! g_mesh a tet mesh.  Otherwise, we make it a degenerate
    ! hex mesh
    tetmesh=.true.
    do i=1,e_mesh%num_eblk
       if (size(e_mesh%eblk(i)%connect,1).ne.4) then
          tetmesh=.false.
          exit
       endif
    enddo
    if (tetmesh) then
       allocate(g_mesh%node_elt(4,g_mesh%nelt))
       allocate(g_mesh%block_elt(g_mesh%nelt))
       ielt=0
       do i=1,e_mesh%num_eblk
          do j=1,e_mesh%eblk(i)%num_elem
             ielt=ielt+1
             g_mesh%node_elt(:,ielt)=e_mesh%eblk(i)%connect(:,j)
             g_mesh%block_elt(ielt)=e_mesh%eblk(i)%ID
          enddo
       enddo
       if (ielt.ne.g_mesh%nelt) then
          print*,'Number of elements in mesh does not add up!'
          stop
       endif
    else
       allocate(g_mesh%node_elt(8,g_mesh%nelt))
       allocate(g_mesh%block_elt(g_mesh%nelt))
       ielt=0
       do i=1,e_mesh%num_eblk
          nodes_per_elt=size(e_mesh%eblk(i)%connect,1)
          if (nodes_per_elt.eq.8) then
             ! Nondegenerate hexes
             do j=1,e_mesh%eblk(i)%num_elem
                ielt=ielt+1
                g_mesh%node_elt(:,ielt)=e_mesh%eblk(i)%connect(:,j)
                g_mesh%block_elt(ielt)=e_mesh%eblk(i)%ID
             enddo
          elseif (nodes_per_elt.eq.6) then
             ! prisms
             do j=1,e_mesh%eblk(i)%num_elem
                ielt=ielt+1
                g_mesh%node_elt(1,ielt)=e_mesh%eblk(i)%connect(1,j)
                g_mesh%node_elt(2,ielt)=e_mesh%eblk(i)%connect(4,j)
                g_mesh%node_elt(3,ielt)=e_mesh%eblk(i)%connect(5,j)
                g_mesh%node_elt(4,ielt)=e_mesh%eblk(i)%connect(2,j)
                g_mesh%node_elt(5,ielt)=e_mesh%eblk(i)%connect(3,j)
                g_mesh%node_elt(6,ielt)=e_mesh%eblk(i)%connect(6,j)
                g_mesh%node_elt(7,ielt)=e_mesh%eblk(i)%connect(6,j)
                g_mesh%node_elt(8,ielt)=e_mesh%eblk(i)%connect(3,j)
                
                g_mesh%block_elt(ielt)=e_mesh%eblk(i)%ID
             enddo
          elseif (nodes_per_elt.eq.5) then
             ! pyramids
             do j=1,e_mesh%eblk(i)%num_elem
                ielt=ielt+1
                g_mesh%node_elt(1:5,ielt)=e_mesh%eblk(i)%connect(:,j)
                g_mesh%node_elt(6:8,ielt)=e_mesh%eblk(i)%connect(5,j)
                
                g_mesh%block_elt(ielt)=e_mesh%eblk(i)%ID
             enddo
          elseif (nodes_per_elt.eq.4) then
             ! tets
             do j=1,e_mesh%eblk(i)%num_elem
                ielt=ielt+1
                g_mesh%node_elt(1,ielt)=e_mesh%eblk(i)%connect(1,j)
                g_mesh%node_elt(2,ielt)=e_mesh%eblk(i)%connect(1,j)
                g_mesh%node_elt(3,ielt)=e_mesh%eblk(i)%connect(2,j)
                g_mesh%node_elt(4,ielt)=e_mesh%eblk(i)%connect(3,j)
                g_mesh%node_elt(5,ielt)=e_mesh%eblk(i)%connect(4,j)
                g_mesh%node_elt(6,ielt)=e_mesh%eblk(i)%connect(4,j)
                g_mesh%node_elt(7,ielt)=e_mesh%eblk(i)%connect(4,j)
                g_mesh%node_elt(8,ielt)=e_mesh%eblk(i)%connect(4,j)
                
                g_mesh%block_elt(ielt)=e_mesh%eblk(i)%ID
             enddo
          endif
       enddo
       if (ielt.ne.g_mesh%nelt) then
          print*,'Number of elements in mesh does not add up!'
          stop
       endif
    endif
    
  end subroutine copy_exodus_mesh_to_gm_mesh

end module grid_mapping_exodus
