!!
!! SIMPL_MESH_TYPE
!!
!! This module provides a derived type that encapsulates the data describing a
!! distributed unstructured mesh that is comprised either entirely of hex cells
!! or entirely of tet cells.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Revised April 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type SIMPL_MESH for encapsulating the data
!!  that describes a distributed unstructured mesh.  Objects of this type are
!!  intended to be used by trusted code, and so its data components are public.
!!  However, the components should be regarded as read-only because an object
!!  may be shared amongst multiple clients.  The following components exist,
!!  though depending on the method used to initialize the object, not all will
!!  necessarily be defined.  On each process the object describes a complete
!!  mesh of some subdomain that references only local entities (cells, faces,
!!  edges, nodes), and can rightly be understood as a serial mesh for that
!!  subdomain.  Globally the mesh-conforming subdomains will overlap, perhaps
!!  only along a boundary but more generally on some collection of cells.
!!  Additional mesh data describes this overlap and provides for communication
!!  between overlapping entities.
!!
!!  Basic mesh data:
!!    cell_type - the type of cells that the mesh is composed of;
!!                either 'HEX' or 'TET'
!!        nnode - the number of nodes in the (subdomain) mesh
!!        nedge - the number of edges in the (subdomain) mesh
!!        nface - the number of faces in the (subdomain) mesh
!!        ncell - the number of cells in the (subdomain) mesh
!!        cnode - the rank-2 integer cell-node array: cnode(k,j) is the node
!!                index of vertex k of cell j.  Its shape is [4,ncell] for a
!!                tet mesh, and [8,ncell] for a hex mesh.
!!        cface - the rank-2 integer cell-face array: cface(k,j) is the face
!!                index of face k of cell j.  Its shape is [4,ncell] for a tet
!!                mesh, and [6,ncell] for a hex mesh.
!!        cfpar - the rank-1 integer cell-face parity bitmask array;
!!                btest(cfpar(j),k) returns true if face cface(k,j) is oriented
!!                inward with respect to cell j.
!!        cedge - the rank-2 integer cell-edge array: cedge(k,j) is the edge
!!                index of edge k of cell j.  Its shape is [6,ncell] for a tet
!!                mesh, and [12,ncell] for a hex mesh.
!!        cepar - the rank-1 integer cell-edge parity bitmask array;
!!                btest(cepar(j),k) returns true if edge cedge(k,j) is oriented
!!                opposite to the local orientation of that edge with respect
!!                to cell j.
!!        fnode - the rank-2 integer face-node array; fnode(k,j) is the node
!!                index of vertex k of face j.  Its shape is [3,nface] for a
!!                tet mesh, and [4,nface] for a hex mesh.
!!        fedge - the rank-2 integer face-edge array; fedge(k,j) is the edge
!!                index of edge k of face j.  Its shape is [3,nface] for a
!!                tet mesh, and [4,nface] for a hex mesh.
!!        fepar - the rank-1 integer face-edge parity bitmask array;
!!                btest(fepar(j),k) returns true if edge fedge(k,j) is oriented
!!                opposite to the local orientation of that edge with respect
!!                to face j.
!!        enode - the rank-2 integer edge-node array; enode(k,j) is the node
!!                index of vertex k of edge j.  Its shape is [2,nedge] for both
!!                tet and hex meshes.
!!        xnode - a rank-1 integer array giving the mapping from local node
!!                indices to their external (global) index (as defined in the
!!                mesh file, for example).  Its shape is [nnode].
!!        xcell - a rank-1 integer array giving the mapping from local cell
!!                indices to their external (global) index (as defined in the
!!                mesh file, for example).  Its shape is [ncell].
!!
!!  Mesh entity sets:
!!       cblock - a rank-1 integer array; cblock(j) is the block ID assigned
!!                to cell j.
!!     block_id - a rank-1 integer array storing the unique block IDs
!!                (globally).  This array is the same on all processes.
!!    cell_set_id - a rank-1 integer array storing the unique cell set IDs.
!!    cell_set_mask - a rank-1 bitmask array; btest(cell_set_mask(j),k)
!!                returns true if cell j belongs to the cell set with ID
!!                cell_set_id(k).
!!    face_set_id - a rank-1 integer array storing the unique face set IDs.
!!    face_set_mask - a rank-1 bitmask array; btest(face_set_mask(j),k)
!!                returns true if face j belongs to the face set with ID
!!                face_set_id(k).  Btest(face_set_mask(j),0) returns true
!!                if face j is a boundary face (global mesh).
!!    node_set_id - a rank-1 integer array storing the unique node set IDs.
!!    node_set_mask - a rank-1 bitmask array; btest(node_set_mask(j),k)
!!                returns true if node j belongs to the node set with ID
!!                node_set_id(k).
!!
!!  Geometry data:
!!            x - the rank-2 real array of node coordinates; x(:,j) is the
!!                position in R^3 of node j.  Its shape is [3,nnode].
!!       length - the rank-1 real array of edge lengths; length(j) is the
!!                length of edge j. Its shape is [nedge].
!!         area - the rank-1 real array of face areas; area(j) is the area
!!                of face j.  Its shape is [nface].
!!       volume - the rank-1 real array of cell volumes; volume(j) is the
!!                volume of cell j.  Its shape is [ncell].
!!       normal - the rank-2 real array of oriented face areas; normal(:,j)
!!                is the oriented area of face j.  Its shape is [3,nface].
!!
!!  Parallel data:
!!    nnode_onP, nedge_onP, nface_onP, ncell_onP - the number of local nodes,
!!        edges, faces, and cells that that are uniquely owned by the process
!!        (on-process).
!!    node_ip, edge_ip, face_ip, cell_ip - derived types that describe the
!!        partitioning and overlap of nodes, edges, faces, and cells, including
!!        information necessary to communicate off-process data between processes.
!!
!!  Link data:
!!    nlink - the number of links in the (subdomain) mesh
!!    lface - a rank-2 integer array; faces lface(1,j) and lface(2,j) are
!!            linked.  Its size is [2,nlink].
!!    link_set_id - a rank-1 integer array storing the unique link set IDs.
!!    link_set_mask - a rank-1 bitmask array; btest(link_set_mask(j),k)
!!                returns true if link j belongs to the link set with ID
!!                link_set_id(k).
!!    nlink_onP -- the number of links that are uniquely owned by the process.
!!    link_ip - a derived type component that describes the partitioning and
!!              overlap of links, and includes info necessary to communicate
!!              off-process link data between processes.
!!

#include "f90_assert.fpp"

module simpl_mesh_type

  use kinds, only: r8
  use base_mesh_class
  use parallel_communication
  use index_partitioning
  use bitfield_type
  implicit none
  private

  type, extends(base_mesh), public :: simpl_mesh
    integer :: nedge=0
    character(:), allocatable :: cell_type
    !! Primary indexing arrays which define the mesh topology (see Note 1)
    integer, pointer :: cnode(:,:) => null()  ! cell nodes
    integer, pointer :: cedge(:,:) => null()  ! cell edges
    integer, pointer :: cface(:,:) => null()  ! cell faces

    !! Secondary indexing arrays derivable from the primary indexing arrays.
    integer, allocatable :: fnode(:,:)  ! face nodes
    integer, allocatable :: fedge(:,:)  ! face edges
    integer, allocatable :: enode(:,:)  ! edge nodes

    !! Cell block ID arrays.
    integer, allocatable :: block_id(:) ! user-assigned ID for each cell block.
    integer, allocatable :: cblock(:)   ! cell block index.

    real(r8), allocatable :: length(:)

    integer :: nedge_onP=0

    !! Partitioning and inter-process communication data.
    type(ip_desc) :: edge_ip
  contains
    procedure :: get_global_cnode_array
    procedure :: get_global_cedge_array
    procedure :: get_global_cface_array
    procedure :: get_global_cblock_array
    procedure :: compute_geometry
    procedure :: write_profile
    final :: simpl_mesh_delete
  end type simpl_mesh

contains

  !! Final subroutine for SIMPL_MESH objects.
  subroutine simpl_mesh_delete (this)
    type(simpl_mesh), intent(inout) :: this
    if (associated(this%cnode)) deallocate(this%cnode)
    if (associated(this%cedge)) deallocate(this%cedge)
    if (associated(this%cface)) deallocate(this%cface)
    call destroy (this%node_ip)
    call destroy (this%edge_ip)
    call destroy (this%face_ip)
    call destroy (this%cell_ip)
  end subroutine simpl_mesh_delete

  !! Compute the geometric data components from the node coordinates.
  subroutine compute_geometry (this)
    use simplex_geometry, only: edge_length, tri_area, tet_volume
    class(simpl_mesh), intent(inout) :: this
    integer :: j
    ASSERT(allocated(this%length))
    ASSERT(allocated(this%area))
    ASSERT(allocated(this%volume))
    do j = 1, this%nedge
      this%length(j) = edge_length(this%x(:,this%enode(:,j)))
    end do
    do j = 1, this%nface
      this%area(j) = tri_area(this%length(this%fedge(:,j)))
    end do
    do j = 1, this%ncell
      this%volume(j) = tet_volume(this%x(:,this%cnode(:,j)))
    end do
  end subroutine compute_geometry

 !! Writes to the tty and output file a profile of the distributed mesh:
 !! numbers of noded, edges, faces, and cells assigned to each processor;
 !! numbers of  on-process and off-process objects.
  subroutine write_profile (this)

    class(simpl_mesh), intent(in) :: this

    integer :: n
    character(len=80) :: line
    integer, dimension(nPE) :: nnode_vec, nedge_vec, nface_vec, ncell_vec
    integer, dimension(2,nPE) :: nvec, evec, fvec, cvec

    call collate (nnode_vec, this%nnode)
    call collate (nedge_vec, this%nedge)
    call collate (nface_vec, this%nface)
    call collate (ncell_vec, this%ncell)

    call broadcast (nnode_vec)
    call broadcast (nedge_vec)
    call broadcast (nface_vec)
    call broadcast (ncell_vec)

    call wline ('  Distributed Mesh Profile:')
    write(line,fmt='(4x,a3,a,4a9)') 'PE', '|', 'nnode', 'nedge', 'nface', 'ncell'
    call wline (line)
    call wline ('    ---+'//repeat('-',36))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,4i9)') n, '|', nnode_vec(n), nedge_vec(n), nface_vec(n), ncell_vec(n)
      call wline (line)
    end do

    if (defined(this%node_ip)) then
      call collate (nvec(1,:), this%node_ip%offP_size())
      call collate (nvec(2,:), this%node_ip%onP_size())
      call broadcast (nvec)
    else
      nvec = 0
    end if

    if (defined(this%edge_ip)) then
      call collate (evec(1,:), this%edge_ip%offP_size())
      call collate (evec(2,:), this%edge_ip%onP_size())
      call broadcast (evec)
    else
      evec = 0
    end if

    if (defined(this%face_ip)) then
      call collate (fvec(1,:), this%face_ip%offP_size())
      call collate (fvec(2,:), this%face_ip%onP_size())
      call broadcast (fvec)
    else
      fvec = 0
    end if

    if (defined(this%cell_ip)) then
      call collate (cvec(1,:), this%cell_ip%offP_size())
      call collate (cvec(2,:), this%cell_ip%onP_size())
      call broadcast (cvec)
    else
      cvec = 0
    end if

    call wline ('  Mesh Communication Profile:')
    write(line,fmt='(4x,3x,1x,a11,3a16)')  'Nodes', 'Edges', 'Faces', 'Cells'
    call wline (line)
    write(line,fmt='(4x,a3,a,4a16)') 'PE', '|', ('off-PE   on-PE', n=1,4)
    call wline (line)
    call wline ('    ---+'//repeat('-',64))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,4(i7,i9))') n, '|', nvec(:,n), evec(:,n), fvec(:,n), cvec(:,n)
      call wline (line)
    end do

  contains

    subroutine wline (line)
      use truchas_logging_services
      character(len=*), intent(in) :: line
      call TLS_info (line)
    end subroutine wline

  end subroutine write_profile

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_GLOBAL_CNODE_ARRAY
 !! GET_GLOBAL_CEDGE_ARRAY
 !! GET_GLOBAL_CFACE_ARRAY
 !!
 !! These routines return the so-named global index array that is associated
 !! with a distributed mesh.  The collated global index array is returned on
 !! the IO processor and a 0-sized array on all others.  The returned array
 !! pointer is allocated by these procedures.
 !!
 !! N.B.: Any storage the pointer may have been associated with at entry
 !! is _not_ deallocated before allocating the pointer anew.
 !!

  subroutine get_global_cnode_array (this, cnode)
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cnode(:,:)
    ASSERT(associated(this%cnode))
    ASSERT(defined(this%cell_ip))
    ASSERT(defined(this%node_ip))
    ASSERT(size(this%cnode,2) == this%cell_ip%local_size())
    ASSERT(minval(this%cnode) >= 1)
    ASSERT(maxval(this%cnode) <= this%node_ip%local_size())
    allocate(cnode(size(this%cnode,1),merge(this%cell_ip%global_size(),0,is_IOP)))
    call collate (cnode, this%node_ip%global_index(this%cnode(:,:this%ncell_onP)))
    ASSERT(minval(cnode) >= 1 .and. maxval(cnode) <= this%node_ip%global_size())
  end subroutine get_global_cnode_array

  subroutine get_global_cedge_array (this, cedge)
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cedge(:,:)
    ASSERT(associated(this%cedge))
    ASSERT(defined(this%cell_ip))
    ASSERT(defined(this%edge_ip))
    ASSERT(size(this%cedge,2) == this%cell_ip%local_size())
    ASSERT(minval(this%cedge) >= 1)
    ASSERT(maxval(this%cedge) <= this%edge_ip%local_size())
    allocate(cedge(size(this%cedge,1),merge(this%cell_ip%global_size(),0,is_IOP)))
    call collate (cedge, this%edge_ip%global_index(this%cedge(:,:this%ncell_onP)))
    ASSERT(minval(cedge) >= 1 .and. maxval(cedge) <= this%edge_ip%global_size())
  end subroutine get_global_cedge_array

  subroutine get_global_cface_array (this, cface)
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cface(:,:)
    ASSERT(associated(this%cface))
    ASSERT(defined(this%cell_ip))
    ASSERT(defined(this%face_ip))
    ASSERT(size(this%cface,2) == this%cell_ip%local_size())
    ASSERT(minval(this%cface) >= 1)
    ASSERT(maxval(this%cface) <= this%face_ip%local_size())
    allocate(cface(size(this%cface,1),merge(this%cell_ip%global_size(),0,is_IOP)))
    call collate (cface, this%face_ip%global_index(this%cface(:,:this%ncell_onP)))
    ASSERT(minval(cface) >= 1 .and. maxval(cface) <= this%face_ip%global_size())
  end subroutine get_global_cface_array

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_GLOBAL_CBLOCK_ARRAY
 !!
 !! This routine returns the global cell block ID array that is associated
 !! with a distributed mesh.  The collated global array is returned on the IO
 !! processor and a 0-sized array on all others.  The returned array pointer
 !! is allocated by this procedure, and any storage the pointer may have been
 !! associated with at entry is _not_ deallocated before allocating the pointer
 !! anew.
 !!

  subroutine get_global_cblock_array (this, cblock)
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cblock(:)
    ASSERT(allocated(this%cblock))
    ASSERT(defined(this%cell_ip))
    allocate(cblock(merge(this%cell_ip%global_size(),0,is_IOP)))
    call collate (cblock, this%cblock(:this%ncell_onP))
  end subroutine get_global_cblock_array

end module simpl_mesh_type
