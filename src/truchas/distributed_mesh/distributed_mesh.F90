!!
!! DISTRIBUTED_MESH
!!
!! This module provides a parallel data structure that enscapsulates the
!! description of a distributed mesh, and a few general procedures that operate
!! on instances of the structure.  The data structure is currently designed to
!! handle the requirements of mimetic discretizations over tetrahedral and
!! hexahedral meshes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 10 Apr 2007.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type DIST_MESH has the following public components.  Depending
!!  on the method used to instantiate the mesh object, not all components will
!!  necessarily be defined.
!!
!!  o NNODE, NEDGE, NFACE, NCELL are the number of nodes, edges faces, and
!!    cells in the mesh.
!!
!!  o CNODE(:,:) is the cell node array; CNODE(k,j) is the node number of local
!!    node k of cell j.  CNODE is dimensioned (4,NCELL) for a tetrahedral mesh
!!    and (8,NCELL) for a hexahedral mesh.
!!
!!  o CFACE(:,:) is the cell face array; CFACE(k,j) is the face number of local
!!    face k of cell j.  CFACE is dimensioned (4,NCELL) for a tetrahedral mesh
!!    and (6,NCELL) for a hexahedral mesh.
!!
!!  o CFPAR(:) is the cell face parity bitmask; BTEST(CFPAR(j),k) returns true
!!    if the orientation of face CFACE(k,j) is inward with respect to cell j.
!!
!!  o CEDGE(:,:) is the cell edge array; CEDGE(k,j) is the edge number of local
!!    edge k of cell j.  CEDGE is dimensioned (6,NCELL) for a tetrahedral mesh
!!    and (12,NCELL) for a hexahedral mesh.
!!
!!  o CEPAR(:) is the cell edge parity bitmask; BTEST(CEPAR(j),k) returns true
!!    if the orientation of edge CEDGE(k,j) is opposite to the local orientation
!!    of that edge with respect to cell j.
!!
!!  o FNODE(:,:) is the face node array; FNODE(k,j) is the node number of local
!!    node k of face j.  FNODE is dimensioned (3,NFACE) for a tetrahedral mesh,
!!    which has triangular faces, and (4,NFACE) for a hexahedral mesh, which
!!    has quad faces.
!!
!!  o FEDGE(:,:) is the face edge array; FEDGE(k,j) is the edge number of local
!!    edge k of face j.  FEDGE is dimensioned (3,NFACE) for a tetrahedral mesh,
!!    which has triangular faces, and (4,NFACE) for a hexahedral mesh, which
!!    has quad faces.
!!
!!  o FEPAR(:) is the face edge parity bitmask; BTEST(FEPAR(j),k) returns true
!!    if the orientation of edge FEDGE(k,j) is opposite to the local orientation
!!    of that edge with respect to face j.
!!
!!  o ENODE(:,:) is the edge node array; ENODE(k,j) is the node number of local
!!    node k of edge j.  ENODE is dimensioned (2,NEDGE).
!!
!!  o XNODE(:) is the external node numbering array; XNODE(j) is the external
!!    number (as defined in the mesh file, e.g.) of local node j.
!!
!!  o XCELL(:) is the external cell numbering array; XCELL(j) is the external
!!    number (as defined in the mesh file, e.g.) of local cell j.
!!
!!  o CBLOCK(:) is the cell block ID array; CBLOCK(j) is the block ID assigned
!!    to cell j.  BLOCK_ID(:) is the list of unique IDs in no particular order.
!!    IDs must be positive integers.  These are legacy arrays tied to Cubit's
!!    element blocks; their use is deprecated in favor of the CELL_SET_* arrays.
!!
!!  o CELL_SET_ID(:) is the list of unique cell set IDs (positive integers), and
!!    CELL_SET_MASK(:) is the cell set bitmask array; BTEST(CELL_SET_MASK(j),k)
!!    returns true if cell j belongs to the cell set with ID CELL_SET_ID(k).
!!
!!  o FACE_SET_ID(:) is the list of unique face set IDs (positive integers), and
!!    FACE_SET_MASK(:) is the face set bitmask array; BTEST(FACE_SET_MASK(j),k)
!!    returns true if face j belongs to the face set with ID FACE_SET_ID(k).
!!    BTEST(FACE_SET_MASK(J),0) returns true if face j is a boundary face.
!!
!!  o X(:,:) is the array of node coordinates, dimensioned (3,NNODE).
!!
!!  o LENGTH(:) is the vector of edge lengths, dimensioned (NEDGE).
!!
!!  o AREA(:) is the vector of face areas, dimensioned (NFACE).
!!
!!  o VOLUME(:) is the vector of cell volumes, dimensioned (NCELL).  This may
!!    be a signed volume that depends on the intrinsic orientation of the cell
!!    given by the local node order of the cell.
!!
!!  o NORMAL(:,:) is the array of oriented face areas; NORMAL(:,j) is the
!!    oriented face area of face j.
!!
!!  o CORNER_VOLUME(:,:) is the array of cell corner volumes; CORNER_VOLUME(k,j)
!!    is the volume of the 'corner' of cell j adjacent to the local vertex k.
!!    This is only relevant to hexahedral meshes, where the 'corner' is the
!!    tetrahedron subtended by the vertex k and the 3 adjacent vertices of the
!!    cell.
!!
!!  o NNODE_ONP, NEDGE_ONP, NFACE_ONP are the number of nodes, edges, and
!!    faces that are owned by the underlying process (on-process).
!!
!!  o NODE_IP, EDGE_IP, FACE_IP, and CELL_IP are derived types that describe
!!    the partitioning of the node, edge, face and cell index sets, and they 
!!    include information necessary to communication certain off-process
!!    node, edge, face and cell data between processes.
!!
!! UTILITY PROCEDURES
!!
!!  CALL GET_FACE_SET_IDS (MESH, FACES, SETIDS)
!!    TYPE(DIST_MESH), INTENT(IN) :: MESH
!!    INTEGER, INTENT(IN) :: FACES(:)
!!    INTEGER, POINTER :: SETIDS(:)
!!
!!    This procedure returns a list of face set IDs in the array SETIDS.
!!    A set ID is included in the list if and only if some face in the
!!    specified list of faces belongs to the set.  The procedure allocates
!!    the return array SETIDS.  This is a parallel procedure returning the
!!    same result on all processes; the list of (local) faces may be
!!    process-dependent, of course.
!!

#include "f90_assert.fpp"

module distributed_mesh

  use kinds, only: r8
  use parallel_communication
  use index_partitioning
  use bitfield_type

  implicit none
  private

  public :: destroy, write_dist_mesh_profile
  public :: get_global_cnode_array, get_global_cedge_array, get_global_cface_array
  public :: get_global_cblock_array, get_global_x_array, get_global_volume_array
  public :: write_mesh
  public :: get_face_set_IDs

  type, public :: dist_mesh
    integer :: nnode=0, nedge=0, nface=0, ncell=0
    !! Primary indexing arrays which define the mesh topology.
    integer, pointer :: cnode(:,:) => null()  ! cell nodes
    integer, pointer :: cedge(:,:) => null()  ! cell edges
    integer, pointer :: cface(:,:) => null()  ! cell faces

    integer, pointer :: cfpar(:) => null() ! relative cell face orientation (bitfield)
    integer, pointer :: cepar(:) => null() ! relative cell edge orientation (bitfield)

    !! Secondary indexing arrays derivable from the primary indexing arrays.
    integer, pointer :: fnode(:,:) => null()  ! face nodes
    integer, pointer :: fedge(:,:) => null()  ! face edges
    integer, pointer :: enode(:,:) => null()  ! edge nodes

    integer, pointer :: fepar(:) => null() ! relative face edge orientation (bitfield)

    !! Relationship to external numbering.
    integer, pointer :: xnode(:) => null()  ! external node number
    integer, pointer :: xcell(:) => null()  ! external cell number

    !! Cell block ID arrays.
    integer, pointer :: block_id(:) => null() ! user-assigned ID for each cell block.
    integer, pointer :: cblock(:) => null()   ! cell block index.

    !! Cell set arrays.
    integer, pointer :: cell_set_id(:) => null()
    integer, pointer :: cell_set_mask(:) => null()

    !! Face set arrays.
    integer, pointer :: face_set_id(:) => null()
    type(bitfield), pointer :: face_set_mask(:) => null()

    !! Node set arrays; not going to use a BITFIELD until proved necessary.
    integer, pointer :: node_set_id(:) => null()
    integer, pointer :: node_set_mask(:) => null()
    
    !! Mesh interface links.
    integer :: nlink = 0, nlink_onP = 0
    integer, pointer :: lface(:,:) => null()
    integer, pointer :: link_set_id(:) => null()    ! user-assigned ID for each link block
    type(bitfield), pointer :: link_set_mask(:) => null()  ! link block index
    type(ip_desc) :: link_ip

    real(kind=r8), pointer :: x(:,:)    => null()
    real(kind=r8), pointer :: length(:) => null()
    real(kind=r8), pointer :: area(:)   => null()
    real(kind=r8), pointer :: volume(:) => null()

    real(kind=r8), pointer :: normal(:,:) => null()
    real(kind=r8), pointer :: corner_volume(:,:) => null()

    integer :: nnode_onP=0, nedge_onP=0, nface_onP=0, ncell_onP=0

    !! Partitioning and inter-process communication data.
    type(ip_desc) :: node_ip, edge_ip, face_ip, cell_ip
  contains
    procedure :: get_cell_set_bitmask
    procedure :: get_face_set_bitmask
  end type dist_mesh

  interface destroy
    module procedure destroy_dist_mesh
  end interface

contains

  !! Returns a scalar bit mask for use in bit operations with the cell_set_mask
  !! array component.  The corresponding bit is set for each cell set ID given
  !! in the array SETIDS.  STAT returns a non-zero value if an unknown cell set
  !! ID is specified, and the optional allocatable deferred-length character
  !! ERRMSG is assigned an explanatory message if present.

  subroutine get_cell_set_bitmask (this, setids, bitmask, stat, errmsg)
    use string_utilities, only: i_to_c
    class(dist_mesh), intent(in) :: this
    integer, intent(in) :: setids(:)
    integer(kind(this%cell_set_mask)), intent(out) :: bitmask
    integer, intent(out) :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: i, j
    bitmask = 0
    do i = 1, size(setids)
      do j = size(this%cell_set_ID), 1, -1
        if (setids(i) == this%cell_set_ID(j)) exit
      end do
      if (j == 0) then
        stat = 1
        if (present(errmsg)) errmsg = 'unknown cell set ID: ' // i_to_c(setids(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do
    stat = 0
  end subroutine get_cell_set_bitmask

  !! Returns a scalar bit mask for use in bit operations with the face_set_mask
  !! array component.  The corresponding bit is set for each face set ID given
  !! in the array SETIDS.  STAT returns a non-zero value if an unknown face set
  !! ID is specified, and the optional allocatable deferred-length character
  !! ERRMSG is assigned an explanatory message if present.

  subroutine get_face_set_bitmask (this, setids, bitmask, stat, errmsg)
    use bitfield_type
    use string_utilities, only: i_to_c
    class(dist_mesh), intent(in) :: this
    integer, intent(in) :: setids(:)
    type(bitfield), intent(out) :: bitmask
    integer, intent(out) :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: i, j
    bitmask = ZERO_BITFIELD
    do i = 1, size(setids)
      do j = size(this%face_set_ID), 1, -1
        if (setids(i) == this%face_set_ID(j)) exit
      end do
      if (j == 0) then
        stat = 1
        if (present(errmsg)) errmsg = 'unknown face set ID: ' // i_to_c(setids(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do
    stat = 0
  end subroutine get_face_set_bitmask

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DESTROY_DIST_MESH
 !!

  subroutine destroy_dist_mesh (this)

    type(dist_mesh), intent(inout) :: this

    if (associated(this%cnode)) deallocate(this%cnode)
    if (associated(this%cedge)) deallocate(this%cedge)
    if (associated(this%cface)) deallocate(this%cface)
    if (associated(this%fnode)) deallocate(this%fnode)
    if (associated(this%fedge)) deallocate(this%fedge)
    if (associated(this%enode)) deallocate(this%enode)
    if (associated(this%cfpar)) deallocate(this%cfpar)
    if (associated(this%cepar)) deallocate(this%cepar)
    if (associated(this%fepar)) deallocate(this%fepar)

    if (associated(this%xnode)) deallocate(this%xnode)
    if (associated(this%xcell)) deallocate(this%xcell)

    if (associated(this%block_id)) deallocate(this%block_id)
    if (associated(this%cblock)) deallocate(this%cblock)
    if (associated(this%cell_set_id)) deallocate(this%cell_set_id)
    if (associated(this%cell_set_mask)) deallocate(this%cell_set_mask)
    if (associated(this%face_set_ID)) deallocate(this%face_set_ID)
    if (associated(this%face_set_mask)) deallocate(this%face_set_mask)
    if (associated(this%node_set_ID)) deallocate(this%node_set_ID)
    if (associated(this%node_set_mask)) deallocate(this%node_set_mask)

    if (associated(this%x)) deallocate(this%x)
    if (associated(this%length)) deallocate(this%length)
    if (associated(this%area)) deallocate(this%area)
    if (associated(this%volume)) deallocate(this%volume)

    if (associated(this%normal)) deallocate(this%normal)
    if (associated(this%corner_volume)) deallocate(this%corner_volume)

    call destroy (this%node_ip)
    call destroy (this%edge_ip)
    call destroy (this%face_ip)
    call destroy (this%cell_ip)

  end subroutine destroy_dist_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_DIST_MESH_PROFILE
 !!
 !! Writes to the tty and output file a profile of the distributed mesh:
 !! numbers of noded, edges, faces, and cells assigned to each processor;
 !! numbers of  on-process and off-process objects.
 !!

  subroutine write_dist_mesh_profile (mesh)

    type(dist_mesh), intent(in) :: mesh

    integer :: n
    character(len=80) :: line
    integer, dimension(nPE) :: nnode_vec, nedge_vec, nface_vec, ncell_vec
    integer, dimension(2,nPE) :: nvec, evec, fvec, cvec

    call collate (nnode_vec, mesh%nnode)
    call collate (nedge_vec, mesh%nedge)
    call collate (nface_vec, mesh%nface)
    call collate (ncell_vec, mesh%ncell)

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

    if (defined(mesh%node_ip)) then
      call collate (nvec(1,:), offP_size(mesh%node_ip))
      call collate (nvec(2,:), onP_size(mesh%node_ip))
      call broadcast (nvec)
    else
      nvec = 0
    end if

    if (defined(mesh%edge_ip)) then
      call collate (evec(1,:), offP_size(mesh%edge_ip))
      call collate (evec(2,:), onP_size(mesh%edge_ip))
      call broadcast (evec)
    else
      evec = 0
    end if

    if (defined(mesh%face_ip)) then
      call collate (fvec(1,:), offP_size(mesh%face_ip))
      call collate (fvec(2,:), onP_size(mesh%face_ip))
      call broadcast (fvec)
    else
      fvec = 0
    end if

    if (defined(mesh%cell_ip)) then
      call collate (cvec(1,:), offP_size(mesh%cell_ip))
      call collate (cvec(2,:), onP_size(mesh%cell_ip))
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

  end subroutine write_dist_mesh_profile

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

  subroutine get_global_cnode_array (mesh, cnode)
    type(dist_mesh), intent(in) :: mesh
    integer, pointer :: cnode(:,:)
    ASSERT( associated(mesh%cnode) )
    ASSERT( defined(mesh%cell_ip) )
    ASSERT( defined(mesh%node_ip) )
    ASSERT( size(mesh%cnode,2) == local_size(mesh%cell_ip) )
    ASSERT( minval(mesh%cnode) >= 1 )
    ASSERT( maxval(mesh%cnode) <= local_size(mesh%node_ip) )
    call allocate_collated_array (cnode, size(mesh%cnode,1), global_size(mesh%cell_ip))
    call collate (cnode, global_index(mesh%node_ip, mesh%cnode(:,:mesh%ncell_onP)))
    ASSERT( minval(cnode) >= 1 .and. maxval(cnode) <= global_size(mesh%node_ip) )
  end subroutine get_global_cnode_array

  subroutine get_global_cedge_array (mesh, cedge)
    type(dist_mesh), intent(in) :: mesh
    integer, pointer :: cedge(:,:)
    ASSERT( associated(mesh%cedge) )
    ASSERT( defined(mesh%cell_ip) )
    ASSERT( defined(mesh%edge_ip) )
    ASSERT( size(mesh%cedge,2) == local_size(mesh%cell_ip) )
    ASSERT( minval(mesh%cedge) >= 1 )
    ASSERT( maxval(mesh%cedge) <= local_size(mesh%edge_ip) )
    call allocate_collated_array (cedge, size(mesh%cedge,1), global_size(mesh%cell_ip))
    call collate (cedge, global_index(mesh%edge_ip, mesh%cedge(:,:mesh%ncell_onP)))
    ASSERT( minval(cedge) >= 1 .and. maxval(cedge) <= global_size(mesh%edge_ip) )
  end subroutine get_global_cedge_array

  subroutine get_global_cface_array (mesh, cface)
    type(dist_mesh), intent(in) :: mesh
    integer, pointer :: cface(:,:)
    ASSERT( associated(mesh%cface) )
    ASSERT( defined(mesh%cell_ip) )
    ASSERT( defined(mesh%face_ip) )
    ASSERT( size(mesh%cface,2) == local_size(mesh%cell_ip) )
    ASSERT( minval(mesh%cface) >= 1 )
    ASSERT( maxval(mesh%cface) <= local_size(mesh%face_ip) )
    call allocate_collated_array (cface, size(mesh%cface,1), global_size(mesh%cell_ip))
    call collate (cface, global_index(mesh%face_ip, mesh%cface(:,:mesh%ncell_onP)))
    ASSERT( minval(cface) >= 1 .and. maxval(cface) <= global_size(mesh%face_ip) )
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

  subroutine get_global_cblock_array (mesh, cblock)
    type(dist_mesh), intent(in) :: mesh
    integer, pointer :: cblock(:)
    ASSERT( associated(mesh%cblock) )
    ASSERT( defined(mesh%cell_ip) )
    call allocate_collated_array (cblock, global_size(mesh%cell_ip))
    call collate (cblock, mesh%cblock(:mesh%ncell_onP))
  end subroutine get_global_cblock_array

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_GLOBAL_X_ARRAY
 !!
 !! This routine returns the global node position array X that is associated
 !! with a distributed mesh.  The collated position array is returned on
 !! the IO processor and a 0-sized array on all others.  The returned array
 !! pointer is allocated by this procedure, and any storage the pointer may
 !! have been associated with at entry is _not_ deallocated before allocating
 !! the pointer anew.
 !!

  subroutine get_global_x_array (mesh, x)
    type(dist_mesh), intent(in) :: mesh
    real(kind=r8), pointer :: x(:,:)
    ASSERT( associated(mesh%x) )
    ASSERT( defined(mesh%node_ip) )
    call allocate_collated_array (x, size(mesh%x,1), global_size(mesh%node_ip))
    call collate (x, mesh%x(:,:mesh%nnode_onP))
  end subroutine get_global_x_array

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_GLOBAL_VOLUME_ARRAY
 !!
 !! This routine returns the global tet volume array VOLUME that is associated
 !! with a distributed mesh.  The collated volume array is returned on
 !! the IO processor and a 0-sized array on all others.  The returned array
 !! pointer is allocated by this procedure, and any storage the pointer may
 !! have been associated with at entry is _not_ deallocated before allocating
 !! the pointer anew.
 !!

  subroutine get_global_volume_array (mesh, volume)
    type(dist_mesh), intent(in) :: mesh
    real(kind=r8), pointer :: volume(:)
    ASSERT( associated(mesh%volume) )
    ASSERT( defined(mesh%cell_ip) )
    call allocate_collated_array (volume, global_size(mesh%cell_ip))
    call collate (volume, mesh%volume(:mesh%ncell_onP))
  end subroutine get_global_volume_array

  subroutine write_mesh (mesh, label)

    type(dist_mesh),  intent(in) :: mesh
    character(len=*), intent(in) :: label

    !! Stub routine.

  end subroutine write_mesh

  subroutine get_face_set_IDs (mesh, faces, setids)
  
    type(dist_mesh), intent(in) :: mesh
    integer, intent(in) :: faces(:)
    integer, pointer :: setids(:)

    integer :: j, n
    type(bitfield) :: bitmask

    bitmask = ZERO_BITFIELD
    do j = 1, size(faces)
      bitmask = ior(bitmask, mesh%face_set_mask(faces(j)))
    end do
    bitmask = ibclr(bitmask, pos=0) ! clear the boundary flag
    bitmask = global_ior(bitmask)

    !! Create the list of involved side set IDS.
    n = 0 ! count first to allocate
    do j = 1, size(mesh%face_set_id)
      if (btest(bitmask,j)) n = n + 1
    end do
    allocate(setids(n))
    n = 0 ! now store the data
    do j = 1, size(mesh%face_set_id)
      if (btest(bitmask,j)) then
        n = n + 1
        setids(n) = mesh%face_set_id(j)
      end if
    end do

  end subroutine get_face_set_IDs

end module distributed_mesh
