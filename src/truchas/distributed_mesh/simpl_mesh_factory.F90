!!
!! SIMPL_MESH_FACTORY
!!
!! Provides a procedure for instantiating a new SIMPL_MESH object that stores
!! a distributed unstructured tetrahedral mesh.  Information about the mesh
!! is passed using a PARAMETER_LIST object.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! July 2015 (refactored version of distributed_tet_mesh)
!!
!! PROGRAMMING INTERFACE
!!
!!  NEW_SIMPL_MESH(PARAMS, STAT, ERRMSG) returns a pointer to a newly allocated
!!    SIMPL_MESH object that has been initialized according to the information
!!    specified by the parameter list PARAMS.  If an error is encountered, the
!!    integer STAT is assigned a non-zero value and the alloctable deferred-
!!    length character string ERRMSG assigned an explanatory message.  The
!!    following parameter list parameters are recognized:
!!
!!    'mesh-file' -- path to an ExodusII mesh file (required)
!!    'coord-scale-factor' -- a multiplicative scaling factor applied to the
!!                            node coordinates (optional)
!!
!! IMPLEMENTATION NOTES
!!
!! 1. We exploit a property of suitably defined tetrahedral meshes that allows
!!    us to know, a priori without additional data, the relationship between
!!    the local orientation of a face or edge of a tetrahedron to the global
!!    orientation of the face or edge within the mesh.  An oriented edge in
!!    the mesh is an ordered list of its end points [n1 n2] with n1 < n2.
!!    Similarly an oriented triangular face in the mesh is an ordered list
!!    of its 3 vertices [n1 n2 n3] with n1 < n2 < n3.  We extend this to the
!!    mesh cells, describing them as oriented tetrahedra [n1 n2 n3 n4] with
!!    n1 < n2 < n3 < n4.  A tetrahedron is positively oriented if its volume
!!    with respect to the ordering of its vertices is positive, and negatively
!!    oriented if its volume is negative.  The faces of a tetrahedron are, in
!!    order, [n2 n3 n4], [n1 n3 n4], [n1 n2 n4], and [n1 n2 n3].  If the tet
!!    is positively oriented, the first and third faces are oriented outward
!!    with respect to the tet, and the second and fourth are oriented inward.
!!    If the tetrahedron is negatively oriented, the orientations are all
!!    reversed.  Similarly, the edges of a triangular face [n1 n2 n3] are,
!!    in order [n2 n3], [n1 n3], [n1 n2], with the first and third oriented
!!    ccw with respect to the orientation of the triangle, and the second
!!    oriented cw.
!!
!! 2. The original method of generating the cell adjacency data produced
!!    neighbor lists that were sorted.  A graph partitioner generally will
!!    be sensitive to reordering the lists, producing somewhat different,
!!    though equally good, partitions.  Here we have sorted the lists for
!!    no other reason than to be able to verify the new implementation
!!    against the original; this step can be eliminated.
!!

#include "f90_assert.fpp"

module simpl_mesh_factory

  use simpl_mesh_type
  use parameter_list_type
  implicit none
  private

  public :: simpl_mesh, new_simpl_mesh

contains

  function new_simpl_mesh (params, stat, errmsg) result (this)

    use kinds, only: r8
    use exodus_mesh_type
    use exodus_mesh_io
    use simpl_mesh_tools
    use index_partitioning
    use parallel_communication
    use permutations
    use truchas_logging_services

    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simpl_mesh), pointer :: this

    integer :: j, k, n, offset, ncell, nedge, nface, nnode
    integer, allocatable :: cnode(:,:), cedge(:,:), cface(:,:), node_perm(:), cell_perm(:)
    integer, dimension(nPE) :: node_bsize, edge_bsize, face_bsize, cell_bsize
    integer, allocatable :: side_map(:), perm(:), offP_size(:), offP_index(:), cblock(:)
    type(exodus_mesh), target :: mesh
    character(:), allocatable :: mesh_file
    real(r8) :: csf

    this => null()

    !! Read the Exodus mesh file into MESH.
    if (is_IOP) then
      call params%get('mesh-file', mesh_file)
      call TLS_info ('  Reading ExodusII mesh file "' // mesh_file // '"')
      call read_exodus_mesh (mesh_file, mesh, stat, errmsg)
      if (stat /= 0) errmsg = 'error reading mesh file: ' // errmsg
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    ncell = mesh%num_elem
    nnode = mesh%num_node

    !! Ensure the mesh comprises only tet cells.
    if (is_IOP) then
      do n = 1, mesh%num_eblk
        if (mesh%eblk(n)%elem_type(1:3) /= 'TET' .or. size(mesh%eblk(n)%connect,dim=1) /= 4) then
          stat = -1
          errmsg = 'ExodusII mesh includes non-tetrahedral elements'
          exit
        end if
      end do
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    !! Flatten the Exodus element block structure.
    allocate(cnode(4,ncell), cblock(ncell))
    call mesh%get_concat_elem_conn (cnode)
    offset = 0
    do n = 1, mesh%num_eblk
      cblock(offset+1:offset+mesh%eblk(n)%num_elem) = mesh%eblk(n)%id
      offset = offset + mesh%eblk(n)%num_elem
    end do

    !! Account for the different tet side labeling convention used by SIMPL_MESH.
    side_map = [3,1,2,4]
    do n = 1, mesh%num_sset
      do j = 1, mesh%sset(n)%num_side
        mesh%sset(n)%face(j) = side_map(mesh%sset(n)%face(j))
      end do
    end do

    !! Partition and order the cells.
    allocate(cell_perm(ncell))
    if (is_IOP) call partition_cells (cnode, nPE, cell_bsize, cell_perm, stat, errmsg)
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return
    if (is_IOP) then
      !! Reorder the cell-based arrays.
      call reorder (cnode,  cell_perm)
      call reorder (cblock, cell_perm)
      !! Map the values of cell-valued arrays.
      allocate(perm(size(cell_perm)))
      call invert_perm (cell_perm, perm)
      do n = 1, size(mesh%sset)
        mesh%sset(n)%elem = perm(mesh%sset(n)%elem)
      end do
      deallocate(perm)
    end if

    !! Partition and order the nodes.
    allocate(node_perm(nnode))
    if (is_IOP) then
      call partition_facets (cnode, cell_bsize, node_bsize, node_perm)
      !! Reorder the node-based arrays.
      call reorder (mesh%coord, node_perm)
      !! Map the values of node-valued arrays.
      allocate(perm(size(node_perm)))
      call invert_perm (node_perm, perm)
      do j = 1, size(cnode,2)
        do k = 1, size(cnode,1)
          cnode(k,j) = perm(cnode(k,j))
        end do
      end do
      do n = 1, size(mesh%nset)
        mesh%nset(n)%node = perm(mesh%nset(n)%node)
      end do
      deallocate(perm)
    end if

    !! Reorient the cells before enumerating the mesh faces and edges; see Note 1.
    if (is_IOP) call reorient_tets (cnode, mesh%sset)

    !! Enumerate and partition the mesh edges.
    allocate(cedge(6,ncell))
    if (is_IOP) then
      call label_tet_mesh_edges (cnode, cedge, nedge)
      allocate(perm(nedge))
      call partition_facets (cedge, cell_bsize, edge_bsize, perm)
      call invert_perm (perm)
      do j = 1, size(cedge,dim=2)
        do k = 1, size(cedge,dim=1)
          cedge(k,j) = perm(cedge(k,j))
        end do
      end do
      deallocate(perm)
    end if

    !! Enumerate and partition the mesh faces.
    allocate(cface(4,ncell))
    if (is_IOP) then
      call label_tet_mesh_faces (cnode, cface, nface)
      allocate(perm(nface))
      call partition_facets (cface, cell_bsize, face_bsize, perm)
      call invert_perm (perm)
      do j = 1, size(cface,dim=2)
        do k = 1, size(cface,dim=1)
          cface(k,j) = perm(cface(k,j))
        end do
      end do
      deallocate(perm)
    end if

    !! Identify off-process cells (ghosts) to include with each partition.
    call select_ghost_cells (cell_bsize, cnode, node_bsize, cedge, edge_bsize, &
                             cface, face_bsize, offP_size, offP_index)

    !! Begin initializing the SIMPL_MESH result.
    allocate(this)

    !! Create the cell index partition; include the off-process cells from above.
    call this%cell_ip%init (cell_bsize, offP_size, offP_index)
    deallocate(offP_size, offP_index)

    this%ncell = this%cell_ip%local_size()
    this%ncell_onP = this%cell_ip%onP_size()

    !! Distribute the cell permutation array; gives mapping to the external cell number.
    allocate(this%xcell(this%ncell))
    call distribute (this%xcell(:this%ncell_onP), cell_perm)
    call gather_boundary (this%cell_ip, this%xcell)
    deallocate(cell_perm)

    !! Create the node index partition and localize the global CNODE array,
    !! which identifies off-process nodes to augment the partition with.
    call this%node_ip%init (node_bsize)
    call localize_index_array (cnode, this%cell_ip, this%node_ip, this%cnode, offP_index)
    call this%node_ip%add_offP_index (offP_index)
    deallocate(offP_index)

    this%nnode = this%node_ip%local_size()
    this%nnode_onP = this%node_ip%onP_size()

    !! Distribute the node permutation array; gives mapping to the external node number.
    allocate(this%xnode(this%nnode))
    call distribute (this%xnode(:this%nnode_onP), node_perm)
    call gather_boundary (this%node_ip, this%xnode)
    deallocate(node_perm)

    !! Create the edge index partition and localize the global CEDGE array,
    !! which identifies off-process edges to augment the partition with.
    call this%edge_ip%init (edge_bsize)
    call localize_index_array (cedge, this%cell_ip, this%edge_ip, this%cedge, offP_index)
    call this%edge_ip%add_offP_index (offP_index)
    deallocate(cedge, offP_index)

    this%nedge = this%edge_ip%local_size()
    this%nedge_onP = this%edge_ip%onP_size()

    !! Create the face index partition and localize the global CFACE array,
    !! which identifies off-process faces to augment the partition with.
    call this%face_ip%init (face_bsize)
    call localize_index_array (cface, this%cell_ip, this%face_ip, this%cface, offP_index)
    call this%face_ip%add_offP_index (offP_index)
    deallocate(offP_index)

    this%nface = this%face_ip%local_size()
    this%nface_onP = this%face_ip%onP_size()

    !! Distribute the element block data
    if (is_IOP) n = mesh%num_eblk
    call broadcast (n)
    allocate(this%block_id(n))
    if (is_IOP) this%block_id = mesh%eblk%id
    call broadcast (this%block_id)

    allocate(this%cblock(this%ncell))
    call distribute (this%cblock(:this%ncell_onP), cblock)
    call gather_boundary (this%cell_ip, this%cblock)
    deallocate(cblock)

    !! Initialize the secondary indexing arrays.
    call init_face_node_data (this)
    call init_face_edge_data (this)
    call init_edge_node_data (this)

    !! Initialize the node, face, and cell set data.
    call init_face_set_data (this, mesh, cface)
    call init_node_set_data (this, mesh)
    call init_cell_set_data (this, mesh)
    deallocate(cface)

    !! Scale the node coordinates and distribute.
    if (is_IOP) then
      call params%get ('coord-scale-factor', csf, default=1.0_r8)
      if (csf /= 1.0_r8) mesh%coord = csf * mesh%coord
    else
      allocate(mesh%coord(3,0))
    end if
    allocate(this%x(3,this%nnode))
    call distribute (this%x(:,:this%nnode_onP), mesh%coord)
    call gather_boundary (this%node_ip, this%x)

    !! Initialize the mesh geometry data components.
    allocate(this%length(this%nedge), this%area(this%nface), this%volume(this%ncell))
    call this%compute_geometry

  end function new_simpl_mesh

  !! This auxiliary subroutine partitions the mesh cells and generates a new
  !! cell numbering for which the partition becomes a block partition and the
  !! cells within partitions are well-ordered.  This is a serial procedure
  !! that operates on the entire mesh.  CNODE(:,:) is the cell connectivity
  !! array, with the last index running over cells.  The cells are partitioned
  !! into NPART groups, and BSIZE(:) returns the number of cells in each of the
  !! partitions.  Under the new cell numbering, the first BSIZE(1) cells belong
  !! to the first partition, the next BSIZE(2) cells to the second, and so on.
  !! The new cell numbering is returned in the permutation array PERM, which
  !! maps new cell numbers to old.  STAT and ERRMSG return a nonzero value and
  !! explanatory message if an error occurs.
  !!
  !! NB: The caller is responsible for applying the renumbering.  That means
  !! permuting all cell-based arrays, including CNODE, and mapping the values
  !! of all cell-valued arrays.
  !!
  !! This implementation uses the cell adjacency graph (cells are adjacent if
  !! they share a face) as input to the partitioner, which seeks to minimize
  !! the number of graph edges connecting different partitions; this is a
  !! typical domain decomposition.
  !!
  !! Currently this is hardwired to use Chaco to do the partitioning (with
  !! hardwired Chaco parameters), but the framework is more or less in place
  !! to allow the input of partitioner choice and parameters.
  !!
  !! This implementation does not currently do anything about well-ordering
  !! the cells within a partition (FIXME!); RCM would be better than nothing.

  subroutine partition_cells (cnode, npart, bsize, perm, stat, errmsg)

    use permutations
    use simpl_mesh_tools, only: get_tet_neighbor_array
    use graph_partitioner_factory
    use parameter_list_type
    use string_utilities, only: i_to_c

    integer, intent(in)  :: cnode(:,:)    ! the cell node array
    integer, intent(in)  :: npart         ! number of partitions
    integer, intent(out) :: bsize(:)      ! partition block size
    integer, intent(out) :: perm(:)       ! permutation
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, ncell
    integer, allocatable :: cnhbr(:,:), xadj(:), adjncy(:)
    type(parameter_list) :: params
    class(graph_partitioner), allocatable :: gpart
    integer, allocatable :: part(:)
    real, allocatable :: ewgt(:)

    ASSERT(npart > 0)
    ASSERT(size(bsize) == npart)
    ASSERT(size(perm) == size(cnode,2))

    ncell = size(cnode,dim=2)

    if (npart == 1) then  ! single partition (serial mode)

      bsize(1) = ncell
      perm = identity_perm(ncell)

    else

      !! Generate the cell neighbor data array CNHBR.
      allocate(cnhbr, mold=cnode) ! valid for simplicial cells
      call get_tet_neighbor_array (cnode, cnhbr, stat)
      if (stat /= 0) then
        stat = -1
        errmsg = 'get_tet_neighbor_array: invalid mesh topology detected'
        return
      end if
      !! Sort the neighbor lists; see Note 2.
      do j = 1, size(cnhbr,dim=2)
        call insertion_sort (cnhbr(:,j))
      end do

      !! Pack the neighbor data into the adjacency structure required by the
      !! partitioner.  The CNHBR array is almost what we want, we just need to
      !! squeeze out 0 values marking "no neighbor" for cells on the boundary.
      allocate(xadj(ncell+1))
      adjncy = pack(cnhbr,mask=(cnhbr > 0))
      xadj(1) = 1
      xadj(2:) = count((cnhbr > 0), dim=1)
      do j = 2, size(xadj)
        xadj(j) = xadj(j) + xadj(j-1)
      end do

      !! Call Chaco to partition the cell adjacency graph.
      call params%set ('partitioner', 'chaco')
      call alloc_graph_partitioner (gpart, params)
      allocate(ewgt(size(adjncy)), part(ncell))
      ewgt = 1.0
      call gpart%compute (ncell, xadj, adjncy, ewgt, npart, part, stat)
      if (stat /= 0) then
        stat = -2
        errmsg = 'chaco graph partitioner error: ierr=' // i_to_c(stat)
        return
      end if

      !! Compute the permutation
      call blocked_partition (part, bsize, perm)

    end if

  end subroutine partition_cells

  !! This auxiliary subroutine partitions the mesh facets of one type (nodes,
  !! edges, or faces) and generates a new numbering of them for which the
  !! partition becomes a block partition and the facets within partitions are
  !! well-ordered.  This is a serial procedure that operates on facet data for
  !! the entire mesh.  FACET(:,:) is a cell based array: FACET(:,j) are the
  !! facets of one type belonging to cell j.  It is assumed that the cells
  !! have already been partitioned and renumbered, and that FACET reflects that
  !! new numbering.  BSIZE returns the number of facets in each partition, and
  !! the new facet numbering is returned in the permutation array PERM, which
  !! maps new facet numbers to old.
  !!
  !! NB: The caller is responsible for applying the renumbering.  That means
  !! permuting any facet-based arrays, and mapping the values of all facet-
  !! valued arrays, including FACET itself.
  !!
  !! The facet partition is based on the cell partition given by the cell
  !! partition sizes CELL_BSIZE.  A facet is assigned to the partition of one
  !! of the cells it is adjacent to.  We use a simple greedy algorithm: for
  !! n=1,2,... assign any unassigned facet adjacent to a cell in partition n
  !! to partition n.  NB: While this may lead to a poorly balanced partition
  !! for facets, it is not clear that this significantly impacts performance.
  !!
  !! To order the facets in a partition we rely on the well-orderedness of the
  !! cells.  We merely number them as they are encountered in the FACET array.

  subroutine partition_facets (facet, cell_bsize, bsize, perm)

    use permutations

    integer, intent(in)  :: facet(:,:)    ! cell facets
    integer, intent(in)  :: cell_bsize(:) ! cell partition block sizes
    integer, intent(out) :: bsize(:)      ! facet partition block sizes
    integer, intent(out) :: perm(:)       ! facet permutation

    integer :: j, k, n, offset, part(size(perm)), perm1(size(perm))

    ASSERT(size(bsize) == size(cell_bsize))
    ASSERT(sum(cell_bsize) == size(facet,2))
    ASSERT(minval(facet) == 1 .and. maxval(facet) == size(perm))

    !! Generate a good ordering of the facets: Number the facets consecutively
    !! as they are encountered in FACET.  PERM1 is old-to-good numbering.
    n = 0
    perm1 = 0
    do j = 1, size(facet,2)
      do k = 1, size(facet,1)
        if (perm1(facet(k,j)) /= 0) cycle  ! already numbered
        n = n + 1
        perm1(facet(k,j)) = n
      end do
    end do
    ASSERT(is_perm(perm1))

    call invert_perm (perm1)  ! PERM1 is now good-to-old numbering

    !! Assign facets to partitions: a simple greedy algorithm.
    part = 0
    offset = 0
    do n = 1, size(cell_bsize)
      do j = offset+1, offset+cell_bsize(n)
        do k = 1, size(facet,1)
          if (part(facet(k,j)) == 0) part(facet(k,j)) = n
        end do
      end do
      offset = offset + cell_bsize(n)
    end do

    !! Partition assignment relative to the good ordering.
    call reorder (part, perm1)

    !! Block partition permutation (new-to-good).
    call blocked_partition (part, bsize, perm)

    !! Total facet permutation (new-to-old).
    do j = 1, size(perm)
      perm(j) = perm1(perm(j))
    end do

  end subroutine partition_facets

  !! This auxiliary subroutine reorients each tet cell by sorting its list of
  !! node numbers into increasing order.  This is a critical step that allows
  !! us to establish a fixed relationship between the local orientation of the
  !! local facets of a simplex with their orientation as entities in the mesh.
  !! This would be trivial were it not for the side sets that must be updated.
  !! This is a serial procedure that operates on the entire mesh.

  subroutine reorient_tets (cnode, sset)

    use exodus_mesh_type, only: side_set

    integer, intent(inout) :: cnode(:,:)
    type(side_set), intent(inout) :: sset(:)

    integer :: n, j, k

    !! Temporarily modify the way side set sides are identified,
    !! using the global node number opposite the side.
    !! This is unaffected by the sort that follows.
    do n = 1, size(sset)
      do j = 1, sset(n)%num_side
        sset(n)%face(j) = cnode(sset(n)%face(j),sset(n)%elem(j))
      end do
    end do

    !! Sort the nodes defining each cell into increasing order.
    do j = 1, size(cnode,2)
      call insertion_sort (cnode(:,j))
    end do

    !! Restore the way side set sides are identified,
    !! using the local side index of the cell.
    do n = 1, size(sset)
      do j = 1, sset(n)%num_side
        do k = 1, size(cnode,dim=2)
          if (sset(n)%face(j) == cnode(k,sset(n)%elem(j))) then
            sset(n)%face(j) = k
            exit
          end if
        end do
      end do
    end do

  end subroutine reorient_tets

  !! Given a partition assignment array PART, this auxiliary subroutine
  !! computes the permutation array PERM that makes the partition a block
  !! partition.  The partition sizes are returned in BSIZE.  PERM maps the
  !! new numbering to the original, and preserves the relative order of
  !! elements in the same partition.  Specifically, 1) PART(PERM(:)) is
  !! sorted (non-decreasing), and 2) j < k whenever PERM(j) < PERM(k) and
  !! PART(PERM(j)) = PART(PERM(k)).

  subroutine blocked_partition (part, bsize, perm)

    use permutations

    integer, intent(in)  :: part(:)   ! partition assignment
    integer, intent(out) :: bsize(:)  ! partition block size
    integer, intent(out) :: perm(:)   ! permutation

    integer :: j, n, next(size(bsize))

    ASSERT(size(part) == size(perm))
    ASSERT(minval(part) >= 1 .and. maxval(part) <= size(bsize))

    !! Compute the block size of each partition.
    bsize = 0
    do j = 1, size(part)
      bsize(part(j)) = bsize(part(j)) + 1
    end do

    !! NEXT(j) is the next free cell number for partition j.
    next(1) = 1
    do n = 2, size(bsize)
      next(n) = next(n-1) + bsize(n-1)
    end do

    !! Generate the permutation.
    do j = 1, size(part)
      perm(next(part(j))) = j
      next(part(j)) = next(part(j)) + 1
    end do

    ASSERT(is_perm(perm))

  end subroutine blocked_partition

  !! This auxiliary subroutine identifies the off-process ghost cells that
  !! should be added to each subdomain.  The criterion is that if a cell
  !! contains a face, edge, or node that belongs to a subdomain and the
  !! cell itself does not belong to the subdomain, then it is added to the
  !! subdomain as a ghost cell.  This ensures that every face, edge, and
  !! node belonging to the subdomain will have complete cell support.  The
  !! motivation for stems from the finite element context in which the
  !! equation for a DoF located at a facet depends on calculations over
  !! all the cells that contain that facet.  This is a serial procedure.

  subroutine select_ghost_cells (cell_bsize, cnode, node_bsize, cedge, edge_bsize, &
                                 cface, face_bsize, offP_size, offP_index)

    use integer_set_type
    use parallel_communication, only: is_IOP, nPE

    integer, intent(in) :: cell_bsize(:), node_bsize(:), edge_bsize(:), face_bsize(:)
    integer, intent(in) :: cnode(:,:), cedge(:,:), cface(:,:)
    integer, allocatable, intent(out) :: offP_size(:), offP_index(:)

    integer :: n, offset
    type(integer_set), allocatable :: ghosts(:)

    if (is_IOP) then
      allocate(ghosts(nPE))
      call overlapping_cells (cnode, cell_bsize, node_bsize, ghosts)
      call overlapping_cells (cedge, cell_bsize, edge_bsize, ghosts)
      call overlapping_cells (cface, cell_bsize, face_bsize, ghosts)
      !! Copy the sets into packed array storage
#ifdef INTEL_DPD200362026
      allocate(offP_size(nPE))
      do n = 1, nPE
        offP_size(n) = ghosts(n)%size()
      end do
#else
      offP_size = ghosts%size()
#endif
      allocate(offP_index(sum(offP_size)))
      offset = 0
      do n = 1, nPE
        call ghosts(n)%copy_to_array (offP_index(offset+1:))
        offset = offset + offP_size(n)
      end do
      deallocate(ghosts)
    else
      allocate(offP_size(0), offP_index(0))
    end if

  end subroutine select_ghost_cells

  !! This auxiliary subroutine identifies for each partition those cells that
  !! contain a facet in the partition but that themselves belong to a different
  !! partition.  For a given partition the identitified cell set will be the
  !! minimal one that if added to the cells already belonging to the partition
  !! will ensure that every facet belonging to the partition has complete
  !! cell-support.  The motivation for this stems from the finite element
  !! context in which the equation for a DoF located at a facet depends on
  !! calculations over all the cells that contain that facet.
  !!
  !! The FACET array is a cell-based array: FACET(:,j) are the indices of the
  !! facets of one type (node, edge, or face) contained in cell j.  The cell
  !! partitioning is described by the cell partition block sizes CELL_BSIZE,
  !! and the facet partitioning by their block sizes BSIZE.  The cells
  !! identified for partition n are added to the set XCELLS(n).  Thus this
  !! routine can be called for different FACET arrays, each call adding to
  !! XCELLS the cells for that type of facet.  This is a serial procedure.

  subroutine overlapping_cells (facet, cell_bsize, bsize, xcells)

    use integer_set_type

    integer, intent(in) :: facet(:,:)     ! cell facets
    integer, intent(in) :: cell_bsize(:)  ! cell partition block sizes
    integer, intent(in) :: bsize(:)       ! facet partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, k, n, m, offset

    ASSERT(size(bsize) == size(cell_bsize))
    ASSERT(all(bsize >= 0))
    ASSERT(all(cell_bsize >= 0))
    ASSERT(size(xcells) == size(bsize))
    ASSERT(sum(cell_bsize) == size(facet,2))
    ASSERT(minval(facet) >= 1 .and. maxval(facet) <= sum(bsize))

    offset = 0
    do n = 1, size(cell_bsize)  ! loop over partitions
      do j = offset+1, offset+cell_bsize(n) ! loop over cells in partition N.
        do k = 1, size(facet,1) ! loop over the facets of cell J.
          m = fpart(facet(k,j))
          if (m /= n) call xcells(m)%add (j)
        end do
      end do
      offset = offset + cell_bsize(n)
    end do

  contains

    integer function fpart (n)
      integer, intent(in) :: n
      integer :: m
      m = n
      do fpart = 1, size(bsize)
        if (m <= bsize(fpart)) return
        m = m - bsize(fpart)
      end do
    end function fpart

  end subroutine overlapping_cells

  !! This subroutine initializes the node set data components.  NODE_SET_MASK
  !! is a node-based integer mask array, with bit n > 0 set if the node belongs
  !! to the nth node set.  Bit 0 is set if the node belongs to a boundary face.
  !! NODE_SET_ID stores the user-assigned integer IDs for the side sets, and is
  !! replicated on each process.

  subroutine init_node_set_data (this, exo_mesh)

    use exodus_mesh_type
    use bitfield_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary, scatter_boundary_or

    type(simpl_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: exo_mesh

    integer :: i, j, n, nnode_tot
    integer, allocatable :: node_set_mask(:)
    logical, allocatable :: bnode(:)

    ASSERT(allocated(this%fnode))
    ASSERT(allocated(this%face_set_mask))

    nnode_tot = this%node_ip%global_size()

    !! Initialize the node set data, %NODE_SET_MASK
    allocate(node_set_mask(merge(nnode_tot,0,is_IOP)))
    if (is_IOP) then
      INSIST(size(exo_mesh%nset)+1 <= bit_size(node_set_mask))
      node_set_mask = 0
      do n = 1, size(exo_mesh%nset)
        ASSERT(minval(exo_mesh%nset(n)%node) >= 1)
        ASSERT(maxval(exo_mesh%nset(n)%node) <= size(node_set_mask))
        do i = 1, exo_mesh%nset(n)%num_node
          j = exo_mesh%nset(n)%node(i)
          node_set_mask(j) = ibset(node_set_mask(j),pos=n)
        end do
      end do
    end if

    !! Initialize the distributed node set mask (%NODE_SET_MASK)
    allocate(this%node_set_mask(this%nnode))
    call distribute (this%node_set_mask(:this%nnode_onP), node_set_mask)
    call gather_boundary (this%node_ip, this%node_set_mask)
    deallocate(node_set_mask)

    !! Initialize the list of node set IDs (%NODE_SET_ID)
    if (is_IOP) n = size(exo_mesh%nset)
    call broadcast (n)
    allocate(this%node_set_id(n))
    if (is_IOP) this%node_set_id = exo_mesh%nset%id
    call broadcast (this%node_set_id)

    !! Tag boundary nodes in the node set mask (bit 0).
    allocate(bnode(this%nnode))
    bnode = .false.
    do j = 1, this%nface
      if (btest(this%face_set_mask(j),pos=0)) bnode(this%fnode(:,j)) = .true.
    end do
    call scatter_boundary_or (this%node_ip, bnode)
    where (bnode) this%node_set_mask = ibset(this%node_set_mask, pos=0)
    deallocate(bnode)

  end subroutine init_node_set_data

  !! This subroutine initializes the cell set data components.  CELL_SET_MASK
  !! is a cell-based integer mask array, with bit n > 0 set if the node belongs
  !! to the nth cell set.  Currently, the cells belonging to an Exodus element
  !! block are mapped to a cell set whose ID is the element block ID.
  !! CELL_SET_ID stores the user-assigned integer IDs for the cell sets, and is
  !! replicated on each process.

  subroutine init_cell_set_data (this, exo_mesh)

    use exodus_mesh_type
    use permutations, only: reorder
    use parallel_communication, only: is_IOP, distribute, broadcast, collate
    use index_partitioning, only: gather_boundary

    type(simpl_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: exo_mesh

    integer :: i, n, offset, ncell_tot
    integer, allocatable :: cell_set_mask(:), cell_perm(:)

    ncell_tot = this%cell_ip%global_size()

    !! Generate the global cell_set mask array (original cell ordering)
    allocate(cell_set_mask(merge(ncell_tot,0,is_IOP)))
    if (is_IOP) then
      INSIST(size(exo_mesh%eblk)+1 <= bit_size(cell_set_mask))
      offset = 0
      do n = 1, size(exo_mesh%eblk)
        do i = 1, exo_mesh%eblk(n)%num_elem
          cell_set_mask(offset+i) = ibset(0, pos=n)
        end do
        offset = offset + exo_mesh%eblk(n)%num_elem
      end do
    end if

    !! Reorder the global cell_set_mask to the internal cell ordering.
    allocate(cell_perm(merge(ncell_tot,0,is_IOP)))
    call collate (cell_perm, this%xcell(:this%ncell_onP))
    if (is_IOP) call reorder (cell_set_mask, cell_perm)
    deallocate(cell_perm)

    !! Initialize the distributed cell set mask (%CELL_SET_MASK)
    allocate(this%cell_set_mask(this%ncell))
    call distribute (this%cell_set_mask(:this%ncell_onP), cell_set_mask)
    call gather_boundary (this%cell_ip, this%cell_set_mask)
    deallocate(cell_set_mask)

    !! Initialize the list of cell set IDs (%CELL_SET_ID)
    if (is_IOP) n = size(exo_mesh%eblk)
    call broadcast (n)
    allocate(this%cell_set_id(n))
    if (is_IOP) this%cell_set_id = exo_mesh%eblk%id
    call broadcast (this%cell_set_id)

  end subroutine init_cell_set_data

  !! This subroutine initializes the face set data components.  FACE_SET_MASK
  !! is a face-based integer mask array, with bit n > 0 set if the face belongs
  !! to the nth face set.  CELL_SET_ID stores the user-assigned integer IDs for
  !! the cell sets, and is replicated on each process.  Currently, each Exodus
  !! side set is mapped to a face set whose ID is the side set ID.  An exodus
  !! "side" is a pair of indices identifying a cell and one of its sides.  This
  !! naturally identifies a unique mesh face, but any orientation information
  !! implicit with the side-of-a-cell description is lost.

  subroutine init_face_set_data (this, exo_mesh, cface)

    use bitfield_type
    use exodus_mesh_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary

    type(simpl_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: exo_mesh
    integer, intent(in) :: cface(:,:)

    integer :: i, j, k, n, nface_tot
    type(bitfield), allocatable :: face_set_mask(:)
    integer, allocatable :: tag(:)

    nface_tot = this%face_ip%global_size()

    !! Generate the global face set mask array.
    allocate(face_set_mask(merge(nface_tot,0,is_IOP)))
    if (is_IOP) then
      ASSERT(minval(cface) >= 1)
      ASSERT(maxval(cface) <= size(face_set_mask))
      INSIST(size(exo_mesh%sset)+1 <= bit_size(face_set_mask))
      face_set_mask = ZERO_BITFIELD
      do n = 1, size(exo_mesh%sset)
        ASSERT(minval(exo_mesh%sset(n)%face) >= 1)
        ASSERT(maxval(exo_mesh%sset(n)%face) <= size(cface,dim=1))
        ASSERT(minval(exo_mesh%sset(n)%elem) >= 1)
        ASSERT(maxval(exo_mesh%sset(n)%elem) <= size(cface,dim=2))
        do i = 1, exo_mesh%sset(n)%num_side
          j = cface(exo_mesh%sset(n)%face(i),exo_mesh%sset(n)%elem(i))
          face_set_mask(j) = ibset(face_set_mask(j), pos=n)
        end do
      end do

      !! Count the references to each face.  A count of 1 indicates a boundary
      !! face, 2 an interior face.  Anything else exposes a bad mesh topology.
      allocate(tag(size(face_set_mask)))
      tag = 0
      do j = 1, size(cface,dim=2)
        do k = 1, size(cface,dim=1)
          tag(cface(k,j)) = tag(cface(k,j)) + 1
        end do
      end do
      INSIST(all(tag >= 1 .and. tag <= 2))

      !! Tag boundary faces (bit 0)
      do j = 1, size(tag)
        if (tag(j) == 1) face_set_mask(j) = ibset(face_set_mask(j), pos=0)
      end do
      deallocate(tag)
    end if

    !! Initialize the distributed face set mask (%FACE_SET_MASK)
    allocate(this%face_set_mask(this%nface))
    call distribute (this%face_set_mask(:this%nface_onP), face_set_mask)
    call gather_boundary (this%face_ip, this%face_set_mask)
    deallocate(face_set_mask)

    !! Initialize the list of cell set IDs (%FACE_SET_ID)
    if (is_IOP) n = size(exo_mesh%sset)
    call broadcast (n)
    allocate(this%face_set_id(n))
    if (is_IOP) this%face_set_ID = exo_mesh%sset%id
    call broadcast (this%face_set_id)

  end subroutine init_face_set_data

  !! Initialize the face-node connectivity data component FNODE.
  subroutine init_face_node_data (this)
    use simplex_topology, only: TET_FACE_VERT
    class(simpl_mesh), intent(inout) :: this
    integer :: i, j, k
    allocate(this%fnode(3,this%nface))
    this%fnode = 0
    do j = 1, size(this%cface,dim=2)
      do k = 1, size(this%cface,dim=1)
        i = this%cface(k,j)
        if (this%fnode(1,i) == 0) this%fnode(:,i) = this%cnode(TET_FACE_VERT(:,k),j)
      end do
    end do
    ASSERT(all(this%fnode > 0))
  end subroutine init_face_node_data

  !! Initialize the face-edge connectivity data component FEDGE.
  subroutine init_face_edge_data (this)
    use simplex_topology, only: TET_FACE_EDGE
    class(simpl_mesh), intent(inout) :: this
    integer :: i, j, k
    allocate(this%fedge(3,this%nface))
    this%fedge = 0
    do j = 1, size(this%cface,dim=2)
      do k = 1, size(this%cface,dim=1)
        i = this%cface(k,j)
        if (this%fedge(1,i) == 0) this%fedge(:,i) = this%cedge(TET_FACE_EDGE(:,k),j)
      end do
    end do
    ASSERT(all(this%fedge > 0))
  end subroutine init_face_edge_data

  !! Initialize the edge-node connectivity data component ENODE.
  subroutine init_edge_node_data (this)
    use simplex_topology, only: TET_EDGE_VERT
    class(simpl_mesh), intent(inout) :: this
    integer :: i, j, k
    allocate(this%enode(2,this%nedge))
    this%enode = 0
    do j = 1, size(this%cedge,dim=2)
      do k = 1, size(this%cedge,dim=1)
        i = this%cedge(k,j)
        if (this%enode(1,i) == 0) this%enode(:,i) = this%cnode(TET_EDGE_VERT(:,k),j)
      end do
    end do
    ASSERT(all(this%enode > 0))
  end subroutine init_edge_node_data

  !! When a procedure executed only on the IO processor returns a status
  !! and possible error message, these need to be communicated to the other
  !! processes.  This auxiliary procedure performs that communication.
  !! The IO process value of STAT is broadcast, and if not 0, the IO process
  !! value of ERRMSG, which must be allocated, is broadcast also.

  subroutine broadcast_status (stat, errmsg)
    use parallel_communication, only: is_IOP, broadcast
    integer, intent(inout) :: stat
    character(:), allocatable, intent(inout) :: errmsg
    integer :: n
    call broadcast (stat)
    if (stat /= 0) then
      if (is_IOP) n = len(errmsg)
      call broadcast (n)
      if (.not.is_IOP) allocate(character(len=n)::errmsg)
      call broadcast (errmsg)
    end if
  end subroutine broadcast_status

  !! Standard insertion sort.  Should be moved to a "sort module".

  pure subroutine insertion_sort (array)
    integer, intent(inout) :: array(:)
    integer :: i, j, next
    do j = 2, size(array)
      i = j
      next = array(j)
      do while (i > 1)  ! find location i of next value
        if (array(i-1) <= next) exit
        array(i) = array(i-1)
        i = i - 1
      end do
      if (i /= j) array(i) = next
    end do
  end subroutine insertion_sort

end module simpl_mesh_factory

