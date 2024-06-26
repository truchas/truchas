!!
!! SIMPL_MESH_TYPE
!!
!! This module provides a derived type that encapsulates the data describing
!! a distributed tetrahedral mesh.  The implementation views the mesh as a
!! simplicial complex, and is specialized for mimetic discretizations on the
!! complex involving degrees of freedom centered at nodes, edges, faces, and
!! tets of the mesh, such as that used for electromagnetics.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Revised August 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type SIMPL_MESH that encapsulates the data
!!  describing a distributed unstructured tetrahedral mesh.  On each process
!!  the object describes a complete mesh of some subdomain that references only
!!  local entities (cells, faces, edges, nodes), and can rightly be considered
!!  as a serial mesh for that subdomain.  Globally, the mesh-conforming
!!  subdomains will overlap, perhaps only along a boundary but more generally
!!  on some collection of cells.  Additional mesh data describes this overlap
!!  and provides for communication between overlapping entities.  See also the
!!  comments for the SIMPLEX_TOPOLOGY module which describe the conventions
!!  and properties of the mesh.
!!
!!  Objects of this type are intended to be used by trusted code, and so its
!!  data components are public.  However, the components must be treated as
!!  read-only because an object may be shared amongst multiple clients.  The
!!  following data components are accessible:
!!
!!    nnode, nedge, nface, ncell - the number of nodes, edges, faces, and
!!        cells in the (subdomain) mesh.
!!
!!    cnode - the rank-2 integer cell-node array: cnode(k,j) is the node index
!!        of vertex k of cell j.  Its shape is [4,ncell].
!!
!!    cface - the rank-2 integer cell-face array: cface(k,j) is the face index
!!         of face k of cell j.  Its shape is [4,ncell].
!!
!!    cedge - the rank-2 integer cell-edge array: cedge(k,j) is the edge index
!!        of edge k of cell j.  Its shape is [6,ncell].
!!
!!    fnode - the rank-2 integer face-node array; fnode(k,j) is the node index
!!        of vertex k of face j.  Its shape is [3,nface].
!!
!!    fedge - the rank-2 integer face-edge array; fedge(k,j) is the edge index
!!        of edge k of face j.  Its shape is [3,nface].
!!
!!    enode - the rank-2 integer edge-node array; enode(k,j) is the node index
!!        of vertex k of edge j.  Its shape is [2,nedge].
!!
!!    xnode - a rank-1 integer array giving the mapping from local node indices
!!        to their external (global) index (as defined in the mesh file, for
!!        example).  Its shape is [nnode].
!!
!!    xcell - a rank-1 integer array giving the mapping from local cell indices
!!        to their external (global) index (as defined in the mesh file, for
!!        example).  Its shape is [ncell].
!!
!!  PARALLEL DATA:
!!
!!    nnode_onP, nedge_onP, nface_onP, ncell_onP - the number of local nodes,
!!        edges, faces, and cells that that are uniquely owned (on-process).
!!
!!    node_imap, edge_imap, face_imap, cell_imap - derived types that describe the
!!        partitioning and overlap of nodes, edges, faces, and cells, including
!!        info necessary to communicate off-process data between processes.
!!
!!  MESH ENTITY DATA:
!!
!!    cell_set_id - a rank-1 integer array storing the unique cell set IDs.
!!        This data is replicated on each process.
!!
!!    cell_set_mask - a rank-1 bitmask array: btest(cell_set_mask(j),k)
!!        returns true if cell j belongs to the cell set with ID cell_set_id(k).
!!
!!    face_set_id - a rank-1 integer array storing the unique face set IDs.
!!        This data is replicated on each process.
!!
!!    face_set_mask - a rank-1 bitmask array: btest(face_set_mask(j),k) returns
!!        true if face j belongs to the face set with ID face_set_id(k).
!!        Btest(face_set_mask(j),0) returns true if face j is a boundary face (global mesh).
!!
!!    node_set_id - a rank-1 integer array storing the unique node set IDs.
!!        This data is replicated on each process.
!!
!!    node_set_mask - a rank-1 bitmask array: btest(node_set_mask(j),k) returns
!!        true if node j belongs to the node set with ID node_set_id(k).
!!
!!    block_id - a rank-1 integer array storing the unique element block IDs.
!!        This data is replicated on each process.
!!
!!    cblock - a rank-1 integer array: cblock(j) is the ID of the element block
!!        cell j belongs to.  Its shape is [ncell].
!!
!!  GEOMETRY DATA:
!!
!!    x - the rank-2 real array of node coordinates; x(:,j) is the position in
!!        R^3 of node j.  Its shape is [3,nnode].
!!
!!    length - the rank-1 real array of edge lengths; length(j) is the length
!!        of edge j. Its shape is [nedge].
!!
!!    area - the rank-1 real array of face areas; area(j) is the area of face j.
!!        Its shape is [nface].
!!
!!    volume - the rank-1 real array of signed cell volumes; volume(j) is the
!!        signed volume of cell j.  Its shape is [ncell].
!!
!!    normal - the rank-2 real array of oriented face areas; normal(:,j) is the
!!        oriented area of face j.  Its shape is [3,nface].
!!

#include "f90_assert.fpp"

module simpl_mesh_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use index_map_type
  implicit none
  private

  type, extends(base_mesh), public :: simpl_mesh
    integer :: nedge=0, nedge_onP=0
    character(:), allocatable :: cell_type
    !! Primary indexing arrays which define the mesh topology.
    integer, allocatable :: cnode(:,:), cedge(:,:), cface(:,:)
    !! Secondary indexing arrays derivable from the primary indexing arrays.
    integer, allocatable :: fnode(:,:), fedge(:,:), enode(:,:), fcell(:,:)
    !! Cell block ID arrays.
    integer, allocatable :: block_id(:) ! user-assigned ID for each cell block
    integer, allocatable :: cblock(:)   ! cell block index
    real(r8), allocatable :: length(:)
    !! Partitioning and inter-process communication data.
    type(index_map) :: edge_imap
  contains
    procedure :: get_global_cnode_array
    procedure :: get_global_cedge_array
    procedure :: get_global_cface_array
    procedure :: get_global_fnode_array
    procedure :: get_global_cblock_array
    procedure :: compute_geometry
    procedure :: write_profile
    procedure :: write_faces_vtk
  end type simpl_mesh

contains

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

  !! Creates the global CNODE array on the IO process, 0-sized array on others.
  subroutine get_global_cnode_array (this, cnode)
    use parallel_communication, only: is_IOP, gather
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cnode(:,:)
    allocate(cnode(size(this%cnode,1),merge(this%cell_imap%global_size,0,is_IOP)))
    call gather (this%node_imap%global_index(this%cnode(:,:this%ncell_onP)), cnode)
  end subroutine get_global_cnode_array

  !! Creates the global CEDGE array on the IO process, 0-sized array on others.
  subroutine get_global_cedge_array (this, cedge)
    use parallel_communication, only: is_IOP, gather
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cedge(:,:)
    allocate(cedge(size(this%cedge,1),merge(this%cell_imap%global_size,0,is_IOP)))
    call gather (this%edge_imap%global_index(this%cedge(:,:this%ncell_onP)), cedge)
  end subroutine get_global_cedge_array

  !! Creates the global CFACE array on the IO process, 0-sized array on others.
  subroutine get_global_cface_array (this, cface)
    use parallel_communication, only: is_IOP, gather
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cface(:,:)
    allocate(cface(size(this%cface,1),merge(this%cell_imap%global_size,0,is_IOP)))
    call gather (this%face_imap%global_index(this%cface(:,:this%ncell_onP)), cface)
  end subroutine get_global_cface_array

  !! Creates the global FNODE array on the IO process, 0-sized array on others.
  subroutine get_global_fnode_array (this, fnode)
    use parallel_communication, only: is_IOP, gather
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: fnode(:,:)
    allocate(fnode(size(this%fnode,1),merge(this%face_imap%global_size,0,is_IOP)))
    call gather (this%node_imap%global_index(this%fnode(:,:this%ncell_onP)), fnode)
  end subroutine get_global_fnode_array

  !! Creates the global CBLOCK array on the IO process; 0-sized array on others.
  subroutine get_global_cblock_array (this, cblock)
    use parallel_communication, only: is_IOP, gather
    class(simpl_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cblock(:)
    allocate(cblock(merge(this%cell_imap%global_size,0,is_IOP)))
    call gather (this%cblock(:this%ncell_onP), cblock)
  end subroutine get_global_cblock_array

  !! Writes to the tty and output file a profile of the distributed mesh:
  !! numbers of noded, edges, faces, and cells assigned to each processor;
  !! numbers of  on-process and off-process objects.

  subroutine write_profile (this)

    use parallel_communication, only: nPE, broadcast, gather
    use truchas_logging_services

    class(simpl_mesh), intent(in) :: this

    integer :: n
    character(80) :: line
    integer, dimension(nPE) :: nnode_vec, nedge_vec, nface_vec, ncell_vec
    integer, dimension(2,nPE) :: nvec, evec, fvec, cvec

    call gather (this%nnode, nnode_vec)
    call gather (this%nedge, nedge_vec)
    call gather (this%nface, nface_vec)
    call gather (this%ncell, ncell_vec)

    call broadcast (nnode_vec)
    call broadcast (nedge_vec)
    call broadcast (nface_vec)
    call broadcast (ncell_vec)

    call TLS_info ('  Distributed Mesh Profile:')
    write(line,fmt='(4x,a3,a,4a9)') 'PE', '|', 'nnode', 'nedge', 'nface', 'ncell'
    call TLS_info (line)
    call TLS_info ('    ---+'//repeat('-',36))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,4i9)') n, '|', nnode_vec(n), nedge_vec(n), nface_vec(n), ncell_vec(n)
      call TLS_info (line)
    end do

    call gather (this%node_imap%offp_size, nvec(1,:))
    call gather (this%node_imap%onp_size, nvec(2,:))
    call broadcast (nvec)

    call gather (this%edge_imap%offp_size, evec(1,:))
    call gather (this%edge_imap%onp_size, evec(2,:))
    call broadcast (evec)

    call gather (this%face_imap%offp_size, fvec(1,:))
    call gather (this%face_imap%onp_size, fvec(2,:))
    call broadcast (fvec)

    call gather (this%cell_imap%offp_size, cvec(1,:))
    call gather (this%cell_imap%onp_size, cvec(2,:))
    call broadcast (cvec)

    call TLS_info ('  Mesh Communication Profile:')
    write(line,fmt='(4x,3x,1x,a11,3a16)')  'Nodes', 'Edges', 'Faces', 'Cells'
    call TLS_info (line)
    write(line,fmt='(4x,a3,a,4a16)') 'PE', '|', ('off-PE   on-PE', n=1,4)
    call TLS_info (line)
    call TLS_info ('    ---+'//repeat('-',64))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,4(i7,i9))') n, '|', nvec(:,n), evec(:,n), fvec(:,n), cvec(:,n)
      call TLS_info (line)
    end do

  end subroutine write_profile

  subroutine write_faces_vtk(this, mask, file, comment)

    use parallel_communication, only: nPE, is_IOP, this_PE, broadcast, gather, global_sum

    class(simpl_mesh), intent(in) :: this
    logical,      intent(in) :: mask(:) ! faces to be written
    character(*), intent(in) :: file    ! file path
    character(*), intent(in) :: comment ! arbitrary vtk description string

    integer :: lun, j, n, sizes(nPE), offset, ntot
    integer, allocatable :: map(:), fnode_loc(:,:), fnode_all(:,:)
    real(r8), allocatable :: x_loc(:,:), x_all(:,:)

    if (is_IOP) then
      open(newunit=lun,file=file,action='write',position='rewind',status='replace')
      write(lun,'("# vtk DataFile Version 3.0")')
      write(lun,'(a)') comment
      write(lun,'("ASCII")')
      write(lun,'("DATASET POLYDATA")')
    end if

    !! Mark nodes that belong to one of the specified faces.
    allocate(map(this%nnode))
    map = 0
    do j = 1, this%nface_onP
      if (mask(j)) map(this%fnode(:,j)) = 1
    end do
    call this%node_imap%scatter_offp_max(map)

    n = sum(map(:this%nnode_onP))
    call gather(n, sizes)
    call broadcast(sizes)
    ntot = sum(sizes)
    offset = sum(sizes(:this_PE-1)) ! node numbering offset

    !! Collate the nodes to be written onto the IO process,
    !! and generate the global numbering of that node subset.
    allocate(x_loc(3,n))
    n = 0
    do j = 1, this%nnode_onP
      if (map(j) == 1) then
        n = n + 1
        x_loc(:,n) = this%x(:,j)
        map(j) = n + offset
      end if
    end do
    call this%node_imap%gather_offp(map)
    allocate(x_all(3,merge(ntot,0,is_IOP)))
    call gather(x_loc, x_all)

    !! Write the node coordinate data.
    if (is_IOP) then
      write(lun,'("POINTS ",i0," double")') ntot
      write(lun,'((3es13.5))') x_all
    end if
    deallocate(x_all, x_loc)

    !! Collate the faces to be written onto the IO process.
    !! The face connectivity data is remapped to the node subset numbering.
    n = count(mask(:this%nface_onP))
    ntot = global_sum(n)
    allocate(fnode_loc(3,n))
    n = 0
    do j = 1, this%nface_onP
      if (mask(j)) then
        n = n + 1
        fnode_loc(:,n) = map(this%fnode(:,j)) - 1 ! vtk uses 0-based numbering
      end if
    end do
    allocate(fnode_all(3,merge(ntot,0,is_IOP)))
    call gather(fnode_loc, fnode_all)

    !! Write the face connectivity data.
    if (is_IOP) then
      write(lun,'("POLYGONS ",i0,1x,i0)') ntot, 4*ntot
      write(lun,'(("3",3(1x,i0)))') fnode_all
    end if

    if (is_IOP) close(lun)

  end subroutine write_faces_vtk

end module simpl_mesh_type
