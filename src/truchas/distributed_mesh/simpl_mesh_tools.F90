!!
!! SIMPL_MESH_TOOLS
!!
!! This module provides a collection of procedures for generating higher-order
!! topological information for unstructured simplicial meshes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Refactored August 2015
!!

#include "f90_assert.fpp"

module simpl_mesh_tools

  implicit none
  private

  public :: get_tri_neighbor_array, get_tet_neighbor_array
  public :: label_tri_mesh_edges, label_tet_mesh_faces, label_tet_mesh_edges

contains

  !! Returns the cell neighbor array for a 2D triangular mesh. CNODE(:,j)
  !! are the three nodes of cell j in the CCW orientation.  CNHBR(k,j) is
  !! the index of the cell across edge k of cell j (edge k is the edge
  !! opposite vertex k).  The procedure validates the mesh topology, and
  !! if faulty topology is detected, STAT returns a non-zero value, and
  !! the corresponding elements of CNHBR assigned a negative value.

  subroutine get_tri_neighbor_array (cnode, cnhbr, stat)
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out) :: stat
    call get_cell_neighbor_array (cnode, get_tri_edge, seg_parity, cnhbr, stat)
  end subroutine get_tri_neighbor_array

  !! NB: Ideally, the following two procedures should be internal procedures
  !! of the preceding subroutine (F2008).

  !! Return the nodes on the specified edge of a triangle.  The edge must
  !! be oriented CCW with respect to the triangle.

  pure subroutine get_tri_edge (cell, face_index, face)
    integer, intent(in) :: cell(:), face_index
    integer, allocatable, intent(inout) :: face(:)
    integer, parameter :: TRI_FACE_VERT(2,3) = reshape(source=[2,3, 3,1, 1,2], shape=[2,3])
    face = cell(TRI_FACE_VERT(:,face_index))
  end subroutine get_tri_edge

  !! Returns the relative parity of oriented segments A and B: 1 if A and B
  !! are the same segment with the same orientation; -1 if A and B are the
  !! same segment with opposite orientations; and 0 otherwise.

  pure integer function seg_parity (a, b) result (parity)
    integer, intent(in) :: a(:), b(:)
    parity = 0
    if (a(1) == b(1)) then
      if (a(2) == b(2)) parity = 1
    else if (a(1) == b(2)) then
      if (a(2) == b(1)) parity = -1
    end if
  end function seg_parity

  !! Returns the cell neighbor array for a 3D tetrahedral mesh. CNODE(:,j)
  !! are the four nodes defining cell j.  CNHBR(k,j) is the index of the cell
  !! across face k of cell j (face k is the face opposite vertex k).  The
  !! procedure validates the mesh topology, and if faulty topology is detected,
  !! STAT returns a non-zero value, and the corresponding elements of CNHBR
  !! assigned a negative value.
  !! NB: The cells should be oriented such that their volume is positive

  subroutine get_tet_neighbor_array (cnode, cnhbr, stat)
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out) :: stat
    call get_cell_neighbor_array (cnode, get_tet_face, tri_parity, cnhbr, stat)
  end subroutine get_tet_neighbor_array

  !! NB: Ideally, the following two procedures should be internal procedures
  !! of the preceding subroutine (F2008).

  !! Return the nodes on the specified face of a tetrahedron.  The face must
  !! be oriented outward with respect to the tetrahedron.

  pure subroutine get_tet_face (cell, face_index, face)
    integer, intent(in) :: cell(:), face_index
    integer, allocatable, intent(inout) :: face(:)
    integer, parameter :: TET_FACE_VERT(3,4) = &
        reshape(source=[2,3,4, 1,4,3, 1,2,4, 1,3,2], shape=[3,4])
    face = cell(TET_FACE_VERT(:,face_index))
  end subroutine get_tet_face

  !! Returns the relative parity of oriented triangles A and B: 1 if A and B
  !! are the same triangle with the same orientation; -1 if A and B are the
  !! same triangle with opposite orientations; and 0 otherwise.

  pure integer function tri_parity (a, b) result (parity)
    integer, intent(in) :: a(:), b(:)
    parity = 0
    if (a(1) == b(1)) then
      if (a(2) == b(2)) then
        if (a(3) == b(3)) parity = 1
      else if (a(2) == b(3)) then
        if (a(3) == b(2)) parity = -1
      end if
    else if (a(1) == b(2)) then
      if (a(2) == b(3)) then
        if (a(3) == b(1)) parity = 1
      else if (a(2) == b(1)) then
        if (a(3) == b(3)) parity = -1
      end if
    else if (a(1) == b(3)) then
      if (a(2) == b(1)) then
        if (a(3) == b(2)) parity = 1
      else if (a(2) == b(2)) then
        if (a(3) == b(1)) parity = -1
      end if
    end if
  end function tri_parity

  !! This auxiliary procedure does the heavy lifting for computing the array
  !! of cell neighbors across cell faces.  The procedure is actually quite
  !! general and applies to meshes comprised of a single cell type.  CNODE
  !! stores the cell node connectivity, and the procedure argument GET_FACE
  !! provides the means of getting the faces of a cell.  The procedure checks
  !! the validity of the mesh topology: each cell face is shared with at most
  !! one other cell, and in the latter case the outward orientation of the
  !! face with respect to one cell is opposite to the outward orientation
  !! with respect to the other cell.  If faulty topology is detected, STAT
  !! returns a non-zero value.

  subroutine get_cell_neighbor_array (cnode, get_face, face_parity, cnhbr, stat)

    use facet_hash_type

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out) :: stat

    interface
      !! Must return the specified oriented face of the given cell as an
      !! ordered list of node indices in the allocatable array FACE.  The face
      !! must be oriented outward with respect to the cell.  The returned size
      !! of FACE must be consistent with the size of the face.
      pure subroutine get_face (cell, face_index, face)
        integer, intent(in) :: cell(:), face_index
        integer, allocatable, intent(inout) :: face(:)
      end subroutine
      !! Must return 1 if the oriented faces A and B are the same with the same
      !! orientation; -1 if they are the same but with opposite orientations;
      !! and 0 otherwise.
      pure integer function face_parity (a, b) result (parity)
        integer, intent(in) :: a(:), b(:)
      end function
    end interface

    integer :: i, j, k, n, jj, kk, max_bin_size, nmatch, bad_faces
    integer, allocatable :: face(:), xbin(:), p(:)
    type(facet_hash) :: hpar

    type :: table_entry
      integer :: j, k                 ! cell and side indices
      integer, allocatable :: face(:) ! face node list (outward and normalized)
    end type table_entry
    type(table_entry), allocatable :: bin_table(:)

    ASSERT(size(cnode,dim=2) == size(cnhbr,dim=2))

    n = size(cnhbr)
    allocate(bin_table(n))

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of faces, but is generally much closer to the
    !! lower bound.  Setting N to the number of faces is generous.
    call hpar%init (n, maxval(cnode))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    allocate(xbin(0:n))
    xbin = 0
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        call get_face (cnode(:,j), k, face)
        call hpar%hash (face, n)
        xbin(n+1) = 1 + xbin(n+1)
      end do
    end do
    max_bin_size = maxval(xbin)

    !! Prepare XBIN: bin J will be BIN_TABLE(XBIN(J):XBIN(J+1)-1)
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j-1) + xbin(j)
    end do

    !! Fill the bin table; use XBIN as a temporary to hold the next free
    !! location for each bin.
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        call get_face (cnode(:,j), k, face)
        call hpar%hash (face, n)
        i = xbin(n)
        bin_table(i)%j = j
        bin_table(i)%k = k
        call move_alloc (face, bin_table(i)%face)
        xbin(n) = i + 1
      end do
    end do

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(xbin,1), 1, -1
      xbin(j) = xbin(j-1)
    end do
    xbin(0) = 1

    cnhbr = 0
    bad_faces = 0
    stat = 0
    allocate(p(max_bin_size)) ! temp that stores parity info while searching a bin
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        if (cnhbr(k,j) /= 0) cycle  ! info already assigned
        !! Get a face and its corresponding bin.
        call get_face (cnode(:,j), k, face)
        call hpar%hash (face, n)
        associate (bin => bin_table(xbin(n):xbin(n+1)-1))
          !! Scan bin for *all* matching faces.
          jj = 0
          kk = 0
          nmatch = 0
          do i = 1, size(bin)
            if (bin(i)%j == j) then ! found myself
              p(i) = 1
            else
              p(i) = face_parity(bin(i)%face, face)
              select case (p(i))
              case (-1) ! a good match (if only one)
                nmatch = 1 + nmatch
                jj = bin(i)%j
                kk = bin(i)%k
              case (1)  ! a bad match (wrong orientation)
                nmatch = 1 + nmatch
              end select
            end if
          end do
          !! Store the neighbor information.
          if (nmatch == 1 .and. jj /= 0) then
            !! Found a unique neighbor; assign the neighbor data.
            cnhbr(k,j) = jj   ! my neighbor, and
            cnhbr(kk,jj) = j  ! my neighbor's neighbor (me!)
          else if (nmatch /= 0) then
            !! Bad mesh topology; tag the faces involved.
            bad_faces = 1 + bad_faces
            do i = 1, size(bin)
              if (p(i) /= 0) cnhbr(bin(i)%k,bin(i)%j) = -bad_faces
            end do
            stat = -1
          end if
        end associate
      end do
    end do

  end subroutine get_cell_neighbor_array

  !! Uniquely enumerates the edges of a tri mesh and initializes the CEDGE
  !! array: CEDGE(k,j) is the edge index of edge k of tet j.  NEDGE returns
  !! the number of edges.  The cell node lists CNODE(:,j) must be sorted.

  subroutine label_tri_mesh_edges (cnode, cedge, nedge)
    use simplex_topology, only: TRI_EDGE_VERT
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cedge(:,:)
    integer, intent(out) :: nedge
    integer :: max_edge
    ASSERT(size(cnode,dim=1) == 3)
    ASSERT(size(cedge,dim=1) == 3)
    ASSERT(size(cnode,dim=2) == size(cedge,dim=2))
    max_edge = size(cedge)  ! worst case; realistically much closer to half this
    call label_mesh_facets (cnode, TRI_EDGE_VERT, max_edge, cedge, nedge)
    ASSERT(all(cedge /= 0))
    ASSERT(all(cedge > 0))
  end subroutine label_tri_mesh_edges

  !! Uniquely enumerates the faces of a tet mesh and initializes the CFACE
  !! array: CFACE(k,j) is the face index of face k of tet j.  NFACE returns
  !! the number of faces.  The cell node lists CNODE(:,j) must be sorted.

  subroutine label_tet_mesh_faces (cnode, cface, nface)
    use simplex_topology, only: TET_FACE_VERT
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cface(:,:)
    integer, intent(out) :: nface
    integer :: max_face
    ASSERT(size(cnode,dim=1) == 4)
    ASSERT(size(cface,dim=1) == 4)
    ASSERT(size(cnode,dim=2) == size(cface,dim=2))
    max_face = size(cface)  ! worst case; realistically much closer to half this
    call label_mesh_facets (cnode, TET_FACE_VERT, max_face, cface, nface)
    ASSERT(all(cface /= 0))
    ASSERT(all(cface > 0))
  end subroutine label_tet_mesh_faces

  !! Uniquely enumerates the edges of a tet mesh and initializes the CEDGE
  !! array: CEDGE(k,j) is the edge index of edge k of tet j.  NEDGE returns
  !! the number of edges.  The cell node lists CNODE(:,j) must be sorted.

  subroutine label_tet_mesh_edges (cnode, cedge, nedge)
    use simplex_topology, only: TET_EDGE_VERT
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cedge(:,:)
    integer, intent(out) :: nedge
    integer :: max_edge
    ASSERT(size(cnode,dim=1) == 4)
    ASSERT(size(cedge,dim=1) == 6)
    ASSERT(size(cnode,dim=2) == size(cedge,dim=2))
    max_edge = 2.5*size(cnode,2)
    call label_mesh_facets (cnode, TET_EDGE_VERT, max_edge, cedge, nedge)
    ASSERT(all(cedge /= 0))
    ASSERT(all(cedge > 0))
  end subroutine label_tet_mesh_edges

  !! This auxiliary subroutine enumerates the facets (edges or faces) of a mesh
  !! comprised of cells of a single type.  CNODE(:,j) is the node connectivity
  !! of cell j and FACET_VERT(:,k) is the list of vertices defining oriented
  !! facet k for this type of cell.  MAX_FACET should be an upper bound on the
  !! number of facets in the mesh.  It determines the size of the hash table
  !! used internally.  NFACET returns the number of facets, and the CFACE array
  !! is defined: |CFACE(k,j)| is the facet index of facet k of cell j.  If
  !! CFACE(k,j) > 0 the facet is co-oriented with the local orientation of the
  !! facet, and if CFACE(k,j) < 0 the facet is contra-oriented.  Any zero value
  !! is an error due to a full hash table.
  !! NB: This procedure assumes that CNODE defines a valid mesh topology.
  !! However the procedure GET_CELL_NEIGHBOR_ARRAY does validate the topology
  !! and it is recommended to arrange that it be called first if possible.

  subroutine label_mesh_facets (cnode, facet_vert, max_facet, cfacet, nfacet)

    use facet_table_type

    integer, intent(in)  :: cnode(:,:)
    integer, intent(in)  :: facet_vert(:,:)
    integer, intent(in)  :: max_facet
    integer, intent(out) :: cfacet(:,:)
    integer, intent(out) :: nfacet

    integer :: j, k, node_max
    type(facet_table) :: table

    ASSERT(size(cnode,dim=2) == size(cfacet,dim=2))
    ASSERT(size(cfacet,dim=1) == size(facet_vert,dim=2))
    ASSERT(minval(facet_vert) >= 1)
    ASSERT(maxval(facet_vert) <= size(cnode,dim=1))
    ASSERT(minval(cnode) > 0)

    node_max = maxval(cnode)
    call table%init (max_facet, node_max)

    do j = 1, size(cfacet,dim=2)
      do k = 1, size(cfacet,dim=1)
        call table%get_facet_label (cnode(facet_vert(:,k),j), cfacet(k,j), insert=.true.)
      end do
    end do
    nfacet = table%number_of_facets()

    !call table%write_hash_performance (0)

  end subroutine label_mesh_facets

end module simpl_mesh_tools
