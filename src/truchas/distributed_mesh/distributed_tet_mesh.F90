!!
!! DISTRIBUTED_TET_MESH
!!
!! This module provides a procedure that takes as input a bare, as-imported
!! tetrahedral mesh presented on the IO processor, and completely fleshed-out,
!! distributed mesh object.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!

#include "f90_assert.fpp"

module distributed_tet_mesh

  use kinds
  use parallel_communication
  use index_partitioning
  use dist_mesh_type

  implicit none
  private

  public :: create_dist_tet_mesh

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CREATE_DIST_TET_MESH
 !!
 !! This procedure takes the bare imported mesh MESH of TYPE(EXTERNAL_MESH) presented
 !! on the IO processor, and returns the distributed mesh object THIS, which
 !! is a fully fleshed-out distributed mesh containing all the basic topological
 !! and geometrical information required by a mimetic discretization scheme.
 !! It also returns the permutation vectors NODE_PERM and CELL_PERM which

  subroutine create_dist_tet_mesh (this, mesh)

    use mesh_importer, only: external_mesh
    use simplicial_mesh_support
    use parallel_communication
    use permutations
    use integer_set_type
    use bitfield_type

    type(dist_mesh), intent(out) :: this
    type(external_mesh), intent(inout), target :: mesh

    integer :: j, n, offset, ncell, nedge, nface, p(mesh%ncell)
    integer, pointer :: cnode(:,:), cedge(:,:), cface(:,:), node_perm(:), cell_perm(:)
    integer, dimension(nPE) :: node_bsize, edge_bsize, face_bsize, cell_bsize
    type(bitfield), pointer :: face_set_mask(:)
    integer, pointer :: offP_index(:), array(:)
    integer, allocatable :: edge_perm(:), face_perm(:), offP_size(:)
    type(integer_set), allocatable :: xcells(:)

    ncell = mesh%ncell
    cnode => mesh%cnode

    call allocate_collated_array (cedge, 6, ncell)
    call allocate_collated_array (cface, 4, ncell)

    call allocate_collated_array (node_perm, mesh%nnode)
    call allocate_collated_array (cell_perm, mesh%ncell)

    if (is_IOP) then
      !! Partition and order the cells.
      call organize_cells (nPE, 'CHACO', mesh%cnode, cell_bsize, cell_perm)
      call reorder (mesh%cell_block, cell_perm)  ! array is cell-based.

      !! Partition and order the nodes.
      call organize_facets (mesh%cnode, cell_bsize, node_bsize, node_perm)
      call reorder (mesh%x, node_perm)  ! array is node-based.

      !! Enumerate the edges and faces of the mesh.
      call label_cell_edges (mesh%cnode, cedge, nedge)
      call label_cell_faces (mesh%cnode, cface, nface)

      !! Partition and order the edges and faces.
      allocate(edge_perm(nedge), face_perm(nface))
      call organize_facets (cedge, cell_bsize, edge_bsize, edge_perm)
      call organize_facets (cface, cell_bsize, face_bsize, face_perm)
      deallocate(edge_perm, face_perm)
    end if

    !! Generate the face set mask array.
    call allocate_collated_array (face_set_mask, nface)
    if (is_IOP) then
      p = inverse_perm(cell_perm)
      do n = 1, size(mesh%sset)
        do j = 1, mesh%sset(n)%num_side
          mesh%sset(n)%elem(j) = p(mesh%sset(n)%elem(j))
        end do
      end do
      call generate_face_set_mask (face_set_mask, cface, mesh%sset)
    end if

    !!
    !! At this point, the mesh is fully fleshed out on the IO processor, and
    !! everything block partitioned.  We're now ready to build the distributed
    !! mesh that will be used by the application code.  The inputs to this
    !! section of the procedure are:
    !!  o global cell node, edge, and face arrays: CNODE, CEDGE, CFACE
    !!  o partition block sizes: NODE_BSIZE, EDGE_BSIZE, FACE_BSIZE, CELL_BSIZE
    !!  o node positions X and boundary face mask array BFACE
    !! All are collated arrays.
    !!

    !! Identify off-process cells to include with each partition.
    if (is_IOP) then

      !! Form the sets of off-process cells.
      allocate(xcells(nPE))
      call overlapping_cells (cnode, cell_bsize, node_bsize, xcells)
      call overlapping_cells (cedge, cell_bsize, edge_bsize, xcells)
      call overlapping_cells (cface, cell_bsize, face_bsize, xcells)

      !! Copy the cell sets into packed array storage.
      allocate(offP_size(nPE))
      do n = 1, nPE
        offP_size(n) = xcells(n)%size()
      end do

      n = sum(offP_size)
      allocate(offP_index(n))
      offset = 0
      do n = 1, nPE
        call xcells(n)%copy_to_array (offP_index(offset+1:))
        offset = offset + offP_size(n)
      end do
      deallocate(xcells)

    else
      allocate(offP_index(0), offP_size(0))
    end if

    !! Create the cell index partition; include the off-process cells from above.
    call create (this%cell_ip, cell_bsize, offP_size, offP_index)
    deallocate(offP_size, offP_index)

    !! Create the node index partition and localize the global CNODE array,
    !! which identifies off-process nodes to augment the partition with.
    call create (this%node_ip, node_bsize)
    call localize_index_array (cnode, this%cell_ip, this%node_ip, this%cnode, offP_index)
    call add_offP_index (this%node_ip, offP_index)
    deallocate(offP_index)

    !! Create the edge index partition and localize the global CEDGE array,
    !! which identifies off-process edges to augment the partition with.
    call create (this%edge_ip, edge_bsize)
    call localize_index_array (cedge, this%cell_ip, this%edge_ip, this%cedge, offP_index)
    call add_offP_index (this%edge_ip, offP_index)
    deallocate(cedge, offP_index)

    !! Create the face index partition and localize the global CFACE array,
    !! which identifies off-process faces to augment the partition with.
    call create (this%face_ip, face_bsize)
    call localize_index_array (cface, this%cell_ip, this%face_ip, this%cface, offP_index)
    call add_offP_index (this%face_ip, offP_index)
    deallocate(cface, offP_index)

    !! Local mesh sizes: on-process plus off-process.
    this%nnode = this%node_ip%local_size()
    this%nedge = this%edge_ip%local_size()
    this%nface = this%face_ip%local_size()
    this%ncell = this%cell_ip%local_size()

    !! On-process mesh sizes.
    this%nnode_onP = this%node_ip%onP_size()
    this%nedge_onP = this%edge_ip%onP_size()
    this%nface_onP = this%face_ip%onP_size()
    this%ncell_onP = this%cell_ip%onP_size()

    !! Distribute the cell block ID arrays.
    if (is_IOP) n = size(mesh%block_id)
    call broadcast (n)
    allocate(this%block_id(n))
    if (is_IOP) this%block_id = mesh%block_id
    call broadcast (this%block_id)
    allocate(this%cblock(this%ncell))
    call distribute (this%cblock(:this%ncell_onP), mesh%cell_block)
    call gather_boundary (this%cell_ip, this%cblock)

    !! Setup the cell set arrays.
    allocate(this%cell_set_id(size(this%block_id)))
    this%cell_set_id = this%block_id
    allocate(this%cell_set_mask(this%ncell))
    this%cell_set_mask = 0
    this%cell_set_mask = ibset(this%cell_set_mask, this%cblock)

    !! Distribute the boundary face mask array.
    allocate(this%face_set_mask(this%nface))
    call distribute (this%face_set_mask(:this%nface_onP), face_set_mask)
    call gather_boundary (this%face_ip, this%face_set_mask)
    deallocate(face_set_mask)

    !! Distribute the cell permutation array; gives mapping to the external cell number.
    allocate(this%xcell(this%ncell))
    call distribute (this%xcell(:this%ncell_onP), cell_perm)
    call gather_boundary (this%cell_ip, this%xcell)
    deallocate(cell_perm)

    !! Broadcast the side set IDs to all processors.
    if (is_IOP) n = size(mesh%sset)
    call broadcast (n)
    allocate(this%face_set_ID(n))
    if (is_IOP) this%face_set_ID = mesh%sset%ID
    call broadcast (this%face_set_ID)

    !! Distribute the node permutation array; gives mapping to the external node number.
    allocate(this%xnode(this%nnode))
    call distribute (this%xnode(:this%nnode_onP), node_perm)
    call gather_boundary (this%node_ip, this%xnode)
    deallocate(node_perm)

    !! Distribute the global node position array, and gather the off-PE node positions.
    allocate(this%x(3,this%nnode))
    call distribute (this%x(:,:this%nnode_onP), mesh%x)
    call gather_boundary (this%node_ip, this%x)

    !! Construct the secondary indexing arrays (THIS SUBDOMAIN MESH).
    allocate(this%fnode(3,this%nface), this%fedge(3,this%nface), this%enode(2,this%nedge))
    call assemble_face_node_list (this%cface, this%cnode, this%fnode)
    call assemble_face_edge_list (this%cface, this%cedge, this%fedge)
    call assemble_edge_node_list (this%cedge, this%cnode, this%enode)

    !! Edge lengths (THIS SUBDOMAIN MESH)
    allocate(this%length(this%nedge))
    do j = 1, this%nedge
      this%length(j) = edge_length(this%x(:,this%enode(:,j)))
    end do

    !! Face areas (THIS SUBDOMAIN MESH)
    allocate(this%area(this%nface))
    do j = 1, this%nface
      this%area(j) = tri_area(this%length(this%fedge(:,j)))
    end do

    !! Signed cell volumes (THIS SUBDOMAIN MESH)
    allocate(this%volume(this%ncell))
    do j = 1, this%ncell
      this%volume(j) = tet_volume(this%x(:,this%cnode(:,j)))
    end do

  end subroutine create_dist_tet_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ORGANIZE_CELLS
 !!
 !! This routine organizes the mesh cells in preparation for creating the
 !! distributed mesh data structure.  This involves partitioning the cells,
 !! and generating a new cell numbering in which the partition becomes a
 !! block partition and the cells within partitions are well-ordered.
 !! It takes the cell node array CNODE and the number of partitions NP, and
 !! returns the new order as a permutation vector PERM (new-to-old) and the
 !! partition block size vector BSIZE.
 !!
 !! N.B.: The CNODE array is permuted into the new ordering, but the caller
 !! must ensure that all other cell-based arrays are permuted using the
 !! returned permutation, and it will need to remap all cell-valued arrays
 !! with its inverse.
 !!
 !! The cell adjacency graph is used as the fuel for the partitioner, which
 !! seeks to minimize the number of (graph) edges joining different partitions;
 !! this is a typical domain decomposition.
 !!
 !! For cells belonging to a common partition, well-ordered simply means an
 !! order that gives good memory locality: adjacent cells--which share faces,
 !! edges, and nodes--should be numbered close together.  We use the heuristic
 !! reverse Cuthill-McKee (RCM) algorithm; there may well be better methods.
 !!

  subroutine organize_cells (np, pmeth, cnode, bsize, perm)

    use simplicial_mesh_support, only: assemble_cell_neighbor_list
    use graph_type
    use permutations

    integer, intent(in)    :: np            ! number of partitions
    character(len=*), intent(in) :: pmeth   ! partition method
    integer, intent(inout) :: cnode(:,:)    ! the cell node array
    integer, intent(out)   :: bsize(:)      ! partition block size
    integer, intent(out)   :: perm(:)       ! permutation

    integer :: j, k, ncell, stat
    integer :: pass(size(cnode,2)), cnhbr(size(cnode,1),size(cnode,2)), rcm_perm(size(perm))

    type(graph), allocatable :: g
    integer, allocatable :: xadj(:), adjncy(:)

    ASSERT( np > 0 )
    ASSERT( size(bsize) == np )
    ASSERT( size(perm) == size(cnode,2) )

    ncell = size(cnode,dim=2)

    !! Generate the cell neighbor array
    call assemble_cell_neighbor_list (cnode, cnhbr, stat)
    ASSERT( stat == 0 )

    !! Create the cell adjacency graph
    allocate(g)
    call g%init (ncell)
    do j = 1, ncell
      do k = 1, size(cnode,1)
        if (cnhbr(k,j) > j) call g%add_edge (j, cnhbr(k,j))
      end do
    end do
    call g%get_adjacency (xadj, adjncy)
    deallocate(g)

    if (np > 1) then  ! Assign cells to partitions.

      select case (pmeth)
      case ('CHACO')

        !! Use Chaco to assign the cells to partitions based on the cell
        !! adjacency graph.  Some Chaco control parameters are hard-wired
        !! within the wrapper.  In Chaco, graph nodes are numbered starting
        !! at one (as here), but its C-arrays use 0-based indexing, so we
        !! need to offset the values in XADJ, which point to ADJNCY, by 1.
        call chaco_f90_wrapper (ncell, np, xadj-1, adjncy, pass, stat)

        pass = pass + 1   ! we want 1-based numbering of the partitions.
        ASSERT( minval(pass) == 1 .and. maxval(pass) == np )

      case default ! shouldn't be here

        ASSERT( .false. )

      end select

      call blocked_partition (pass, bsize, perm)

    else  ! single partition (serial mode)

      bsize = 0
      bsize(1) = ncell
      perm = identity_perm(ncell)

    end if

    !! Generate a good ordering of common-partition cells.
    !! The plan is to use the cell neighbor graph and partition assignments
    !! to generate the RCM ordering of the cells in each subgraph defined by
    !! the partitioning.  For now we just leave the order as it is.
    rcm_perm = identity_perm(ncell)

    deallocate(xadj, adjncy)

    !! Total cell permutation (new-to-old)
    do j = 1, size(perm)
      perm(j) = rcm_perm(perm(j))
    end do

    !! Permute the cell node array accordingly.
    call reorder (cnode, perm)

  end subroutine organize_cells

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ORGANIZE_FACETS
 !!
 !! This routine organizes the mesh facets (i.e., nodes, edges, or faces) in
 !! preparation for creating the distributed mesh data structure. This involves
 !! partitioning the facets, and generating a new facet numbering in which the
 !! partition becomes a block partition and the facets within partitions are
 !! well-ordered.  The partition block sizes are returned in the vector BSIZE,
 !! and the permutation (new-to-old) in PERM.
 !!
 !! N.B.: The FACET array values are mapped to the new numbering, but the
 !! caller must remap any other facet-valued arrays, and permute any
 !! facet-based arrays, using the returned permutation.
 !!
 !! The FACET array is a cell-based array: FACET(:,j) are the facets of one
 !! type adjacent to cell j.  It is assumed that the cells have already been
 !! organized and FACET defined accordingly, and that all facets are referenced
 !! by the array.
 !!
 !! The partitioning of the facets is based on the existing cell partition
 !! described by the cell partition block size CELL_BSIZE.  A facet must be
 !! assigned to a partition of one of the cells it is adjacent to.  We use a
 !! simple greedy algorithm: for n=1,2,... assign any unassigned facet adjacent
 !! to a cell in partition n to partition n.  This may lead to a poorly
 !! balanced partition, but this can be addressed later if it turns out to be
 !! a significant issue.
 !!
 !! To order the facets in a common partition we rely on the well-ordered'ness
 !! of the cells.  We merely number them as they are encountered in the FACET
 !! array.
 !!

  subroutine organize_facets (facet, cell_bsize, bsize, perm)

    use permutations

    integer, intent(inout) :: facet(:,:)    ! cell facets
    integer, intent(in)    :: cell_bsize(:) ! cell partition block sizes
    integer, intent(out)   :: bsize(:)      ! facet partition block sizes
    integer, intent(out)   :: perm(:)       ! permutation

    integer :: j, k, n, offset, pass(size(perm)), perm1(size(perm))

    ASSERT( size(bsize) == size(cell_bsize) )
    ASSERT( sum(cell_bsize) == size(facet,2) )
    ASSERT( minval(facet) == 1 .and. maxval(facet) == size(perm) )

    !! Generate a good ordering of the facets: Number the facets consecutively
    !! as they are encountered in FACET.  PERM1 is old-to-good numbering.
    n = 0
    perm1 = 0
    do j = 1, size(facet,2)
      do k = 1, size(facet,1)
        if (perm1(facet(k,j)) /= 0) cycle  ! numbered this one already
        n = n + 1
        perm1(facet(k,j)) = n
      end do
    end do

    ASSERT( is_perm(perm1) )

    call invert_perm (perm1)  ! PERM1 is now good-to-old numbering

    !! Assign facets to partitions: a simple greedy algorithm.
    pass = 0
    offset = 0
    do n = 1, size(cell_bsize)
      do j = offset+1, offset+cell_bsize(n)
        where (pass(facet(:,j)) == 0) pass(facet(:,j)) = n
      end do
      offset = offset + cell_bsize(n)
    end do

    !! Partition assignment relative to the good ordering.
    call reorder (pass, perm1)

    !! Block partition permutation (new-to-good).
    call blocked_partition (pass, bsize, perm)

    !! Total facet permutation (new-to-old).
    do j = 1, size(perm)
      perm(j) = perm1(perm(j))
    end do

    !! Remap the facet array values; PERM1 is old-to-new.
    perm1 = inverse_perm(perm)
    do j = 1, size(facet,2)
      do k = 1, size(facet,1)
        facet(k,j) = perm1(facet(k,j))
      end do
    end do

  end subroutine organize_facets


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BLOCKED_PARTITION
 !!
 !! This is an auxillary routine that takes the partition assignment vector
 !! PASS and returns the permutation PERM that will make the partition a
 !! block partition.  The partition sizes are returned in the vector BSIZE.
 !! The permutation preserves the relative order of elements in the same
 !! partition.  PERM is a mapping from the new numbering to original.
 !!

  subroutine blocked_partition (pass, bsize, perm)

    use permutations

    integer, intent(in)  :: pass(:)   ! partition assignment
    integer, intent(out) :: bsize(:)  ! partition block size
    integer, intent(out) :: perm(:)   ! permutation

    integer :: j, n, next(size(bsize))

    ASSERT( size(pass) == size(perm) )
    ASSERT( minval(pass) >= 1 .and. maxval(pass) <= size(bsize) )

    !! Compute the block size of each partition.
    bsize = 0
    do j = 1, size(pass)
      bsize(pass(j)) = bsize(pass(j)) + 1
    end do

    !! NEXT(j) is the next free cell number for partition j.
    next(1) = 1
    do n = 2, size(bsize)
      next(n) = next(n-1) + bsize(n-1)
    end do

    !! Generate the permutation.
    do j = 1, size(pass)
      perm(next(pass(j))) = j
      next(pass(j)) = next(pass(j)) + 1
    end do

    ASSERT( is_perm(perm) )

  end subroutine blocked_partition

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GENERATE_FACE_SET_MASK
 !!
 !! This routine generates the face-based, bitmask vector FACE_SET_MASK.  Bit 0
 !! is set for faces lying on the boundary of the domain, and bit n, n>0, is
 !! set if the face belongs to the nth side set.  Note that the relationship
 !! between this implicit internal side set numbering and the user-assigned ID
 !! is stored separately.
 !!
 !! N.B.  The number of face sets is limited to 31 for the typical default
 !! 4-byte integer, but this can easily be amended by going to 8-byte integers
 !! or by using a general 'bit field' derived data type.
 !!

  subroutine generate_face_set_mask (face_set_mask, cface, sset)

    use mesh_importer, only: side_set
    use bitfield_type

    type(bitfield), intent(out) :: face_set_mask(:) ! face set bitmask array
    integer,        intent(in)  :: cface(:,:)       ! cell face array
    type(side_set), intent(in)  :: sset(:)          ! side set array

    integer :: i, j, k, tag(size(face_set_mask))

    ASSERT( minval(cface) >= 1 )
    ASSERT( maxval(cface) <= size(face_set_mask) )

    !! Count the references to each face.  A count of 1 indicates a boundary
    !! face, 2 an interior face.  Anything else exposes a bad mesh topology.
    tag = 0
    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        tag(cface(k,j)) = tag(cface(k,j)) + 1
      end do
    end do

    INSIST( all(tag >= 1) )
    INSIST( all(tag <= 2) )

    do j = 1, size(tag)
      if (tag(j) == 1) then
        face_set_mask(j) = ibset(ZERO_BITFIELD,pos=0)
      else
        face_set_mask(j) = ZERO_BITFIELD
      end if
    end do

    !! If a face belongs to a side set, set the corresponding bit.
    INSIST( size(sset) <= bit_size(face_set_mask)-1 )
    do k = 1, size(sset)
      ASSERT( minval(sset(k)%elem) >= 1 )
      ASSERT( maxval(sset(k)%elem) <= size(cface,2) )
      ASSERT( minval(sset(k)%face) >= 1 )
      ASSERT( maxval(sset(k)%face) <= size(cface,1) )
      do j = 1, sset(k)%num_side
        i = cface(sset(k)%face(j),sset(k)%elem(j))
        face_set_mask(i) = ibset(face_set_mask(i),pos=k)
      end do
    end do

  end subroutine generate_face_set_mask

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! This routine identifies for each partition those cells that contain a facet
 !! in the partition but which themselves belong to a different partition.  For
 !! a given partition the identitified cell set will be the minimal one that if
 !! added to the cells already belonging to the partition will ensure that
 !! every facet belonging to the partition has complete cell-support.  The
 !! motivation for this stems from the finite element context in which the
 !! equation for a DoF located at a facet depends on calculations over all the
 !! cells that contain that facet.
 !!
 !! The FACET array is a cell-based array: FACET(:,j) are the indices of the
 !! facets of one type (node, edge, or face) contained in cell j.  The cell
 !! partitioning is described by the cell partition block sizes CELL_BSIZE, and
 !! the facet partitioning by their block sizes BSIZE.  The cells identified
 !! for partition n are added to the set XCELLS(n).  Thus this routine can be
 !! called for different FACET arrays, each call adding to XCELLS the cells
 !! for that type of facet.
 !!
 !! This is a serial procedure.
 !!

  subroutine overlapping_cells (facet, cell_bsize, bsize, xcells)

    use integer_set_type

    integer, intent(in) :: facet(:,:)     ! cell facets
    integer, intent(in) :: cell_bsize(:)  ! cell partition block sizes
    integer, intent(in) :: bsize(:)       ! facet partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, k, n, m, offset

    ASSERT( size(bsize) == size(cell_bsize) )
    ASSERT( all(bsize >= 0) )
    ASSERT( all(cell_bsize >= 0) )
    ASSERT( size(xcells) == size(bsize) )
    ASSERT( sum(cell_bsize) == size(facet,2) )
    ASSERT( minval(facet) >= 1 .and. maxval(facet) <= sum(bsize) )

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

end module distributed_tet_mesh

