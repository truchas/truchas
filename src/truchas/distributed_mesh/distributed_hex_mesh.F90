!!
!! DISTRIBUTED_HEX_MESH
!!
!! This module provides a procedure that takes as input a bare, as-imported
!! hexahedral mesh presented on the IO processor, and returns a completely
!! fleshed-out, distributed mesh object.
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 12 Mar 2004
!!

#include "f90_assert.fpp"

module distributed_hex_mesh

  use kinds
  use parallel_communication
  use index_partitioning
  use distributed_mesh

  implicit none
  private

  public :: create_dist_hex_mesh

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CREATE_DIST_HEX_MESH
 !!
 !! This procedure takes the bare imported mesh MESH of TYPE(EXTERNAL_MESH) presented
 !! on the IO processor, and returns the distributed mesh object THIS, which
 !! is a fully fleshed-out distributed mesh containing all the basic topological
 !! and geometrical information required by a mimetic discretization scheme.
 !! It also returns the permutation vectors NODE_PERM and CELL_PERM which

  subroutine create_dist_hex_mesh (this, mesh)

    use mesh_importer
    use cell_geometry
    use hexahedral_mesh_support
    use parallel_communication
    use permutations
    use sets
    use bitfield_type

    type(dist_mesh), intent(out) :: this
    type(external_mesh), intent(inout), target :: mesh

    integer :: j, k, n, offset, ncell, nface, stat, p(mesh%ncell)
    integer, pointer :: cnode(:,:), cface(:,:), lface(:,:), node_perm(:), cell_perm(:)
    integer, dimension(nPE) :: node_bsize, face_bsize, cell_bsize, bsize
    type(bitfield), pointer :: face_set_mask(:)
    integer, pointer :: node_set_mask(:)
    integer, allocatable :: face_perm(:), offP_size(:), lnhbr(:,:), cnhbr(:,:), perm(:)
    integer, pointer :: cfpar(:), array(:), offP_index(:)
    type(integer_set), allocatable :: xcells(:)
    type(bitfield), pointer :: link_set_mask(:)

    ncell = mesh%ncell
    cnode => mesh%cnode
    
    !call allocate_collated_array (cedge, 12, ncell)
    call allocate_collated_array (cface, 6, ncell)
    call allocate_collated_array (lface, 2, mesh%nlink)
    call allocate_collated_array (cfpar, ncell)

    call allocate_collated_array (node_perm, mesh%nnode)
    call allocate_collated_array (cell_perm, mesh%ncell)

    if (is_IOP) then
    
      !! Partition and order the cells.
      allocate(cnhbr(6,ncell), lnhbr(2,mesh%nlink))
      call cell_neighbor_info (cnode, mesh%lnode, cnhbr, lnhbr, stat)
      INSIST(stat == 0)
      call partition_cells (nPE, 'CHACO', cnhbr, lnhbr, p)
      call blocked_partition (p, cell_bsize, cell_perm)
      
      !! Reorder cell-based arrays.
      call reorder (mesh%cnode, cell_perm)
      call reorder (mesh%cell_block, cell_perm)
      deallocate(cnhbr)
      !call reorder (cnhbr, cell_perm)
      
      !! Map the values of cell-valued arrays.
      call invert_perm (cell_perm, p)
      do j = 1, size(lnhbr,dim=2)
        lnhbr(:,j) = p(lnhbr(:,j))
      end do
      do n = 1, size(mesh%sset)
        mesh%sset(n)%elem = p(mesh%sset(n)%elem)
      end do
      !do j = 1, size(cnhbr,dim=2)
      !  cnhbr(:,j) = p(cnhbr(:,j))
      !end do
      
      !! Partition and order the cells.
      !call organize_cells (nPE, 'CHACO', mesh%cnode, mesh%lnode, cell_bsize, cell_perm)
      !call reorder (mesh%cell_block, cell_perm)  ! array is cell-based.

      !! Partition and order the nodes.
      call organize_facets (mesh%cnode, cell_bsize, node_bsize, node_perm)
      call reorder (mesh%x, node_perm)  ! array is node-based.
      
      !! Map the values of node-valued arrays.
      allocate(perm(mesh%nnode))
      call invert_perm (node_perm, perm)
      do j = 1, size(mesh%lnode,dim=2)
        mesh%lnode(:,j) = perm(mesh%lnode(:,j))
      end do
      do n = 1, size(mesh%nset)
        mesh%nset(n)%node = perm(mesh%nset(n)%node)
      end do
      deallocate(perm)

      !! Label the edges and faces of the mesh.
      !call label_cell_edges (mesh%cnode, cedge, nedge)
      call label_cell_faces (mesh%cnode, cface, nface, mesh%lnode, lface)
      cfpar = 0
      do j = 1, ncell
        n = 0
        do k = 1, 6
          if (cface(k,j) < 0) then
            n = ibset(n,k)
            cface(k,j) = -cface(k,j)
          end if
        end do
        cfpar(j) = n
      end do

      !! Partition and order the edges.
      !allocate(edge_perm(nedge))
      !call organize_facets (cedge, cell_bsize, edge_bsize, edge_perm)
      !deallocate(edge_perm)

      !! Partition and order the faces.
      allocate(face_perm(nface))
      call organize_facets (cface, cell_bsize, face_bsize, face_perm)
      call invert_perm (face_perm)
      do j = 1, mesh%nlink
        lface(:,j) = face_perm(lface(:,j))
      end do
      deallocate(face_perm)
    end if

    !! Generate the boundary face mask array.
    call allocate_collated_array (face_set_mask, nface)
    if (is_IOP) then
      call generate_face_set_mask (face_set_mask, cface, mesh%sset)
      !! We won't call interface faces boundary faces.
      !face_set_mask(lface(1,:)) = ibclr(face_set_mask(lface(1,:)), pos=0)
      !face_set_mask(lface(2,:)) = ibclr(face_set_mask(lface(2,:)), pos=0)
    end if
    
    !! Generate the node set mask array; need to tag boundary nodes later.
    call allocate_collated_array (node_set_mask, mesh%nnode)
    if (is_IOP) call generate_node_set_mask (node_set_mask, mesh%nset)

    !!
    !! At this point, the mesh is fully fleshed out on the IO processor, and
    !! everything block partitioned.  We're now ready to build the distributed
    !! mesh that will be used by the application code.  The inputs to this
    !! section of the procedure are:
    !!  o global cell node, edge, and face arrays: CNODE, CEDGE, CFACE
    !!  o partition block sizes: NODE_BSIZE, EDGE_BSIZE, FACE_BSIZE, CELL_BSIZE
    !!  o node positions X and boundary face mask array FACE_SET_MASK
    !! All are collated arrays.
    !!

    !! Identify off-process cells to include with each partition.
    if (is_IOP) then

      !! Form the sets of off-process cells.
      allocate(xcells(nPE))
      call overlapping_cells (cnode, cell_bsize, node_bsize, xcells)
      !call overlapping_cells (cedge, cell_bsize, edge_bsize, xcells)
      call overlapping_cells (cface, cell_bsize, face_bsize, xcells)
      call overlapping_cells2 (lnhbr, cell_bsize, xcells)

      !! Copy the cell sets into packed array storage.
      allocate(offP_size(nPE))
      do n = 1, nPE
        offP_size(n) = size(xcells(n))
      end do

      n = sum(offP_size)
      allocate(offP_index(n))
      offset = 0
      do n = 1, nPE
        array => to_array(xcells(n))
        offP_index(offset+1:offset+offP_size(n)) = array
        offset = offset + offP_size(n)
        deallocate(array)
      end do

      !! We're finished with the cell sets now.
      do n = 1, size(xcells)
        call clear (xcells(n))
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
    !call create (this%edge_ip, edge_bsize)
    !call localize_index_array (cedge, this%cell_ip, this%edge_ip, this%cedge, offP_index)
    !call add_offP_index (this%edge_ip, offP_index)
    !deallocate(cedge, offP_index)

    !! Create the face index partition and localize the global CFACE array,
    !! which identifies off-process faces to augment the partition with.
    call create (this%face_ip, face_bsize)
    call localize_index_array (cface, this%cell_ip, this%face_ip, this%cface, offP_index)
    call add_offP_index (this%face_ip, offP_index)
    deallocate(cface, offP_index)

    if (is_IOP) then
      allocate(offP_size(nPE), perm(mesh%nlink))
      call partition_links (lnhbr, cell_bsize, bsize, perm, offP_size, offP_index)
      call reorder (lnhbr, perm)
      call reorder (lface, perm)
      call reorder (mesh%lnode, perm)
      call reorder (mesh%link_block, perm)
      call invert_perm (perm)
      offP_index = perm(offP_index)
      deallocate(perm)
    else
      allocate(offP_size(0), offP_index(0))
    end if
    call create (this%link_ip, bsize, offP_size, offP_index)
    deallocate(offP_size, offP_index)
    call localize_index_array (lface, this%link_ip, this%face_ip, this%lface, offP_index)
    INSIST(size(offP_index) == 0)
    deallocate(lface, offP_index)
    
    !! Local mesh sizes: on-process plus off-process.
    this%nnode = local_size(this%node_ip)
    this%nedge = 0 !local_size(this%edge_ip)
    this%nface = local_size(this%face_ip)
    this%ncell = local_size(this%cell_ip)
    this%nlink = local_size(this%link_ip)

    !! On-process mesh sizes.
    this%nnode_onP = onP_size(this%node_ip)
    this%nedge_onP = 0 !onP_size(this%edge_ip)
    this%nface_onP = onP_size(this%face_ip)
    this%ncell_onP = onP_size(this%cell_ip)
    this%nlink_onP = onP_size(this%link_ip)

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
    do j = 1, this%ncell
      do n = size(this%cell_set_id), 1, -1
        if (this%cblock(j) == this%cell_set_id(n)) exit
      end do
      INSIST( n /= 0 )
      this%cell_set_mask(j) = ibset(this%cell_set_mask(j), pos=n)
    end do

    !! Link set arrays.
    n = mesh%nlblock
    call broadcast (n)
    allocate(this%link_set_id(n), this%link_set_mask(this%nlink))
    if (is_IOP) this%link_set_id = mesh%link_block_id
    call broadcast (this%link_set_id)
    call allocate_collated_array (link_set_mask, mesh%nlink)
    if (is_IOP) then
      INSIST(n <= bit_size(link_set_mask)-1)
      link_set_mask = ZERO_BITFIELD
      do j = 1, mesh%nlink
        do n = size(mesh%link_block_id), 1, -1
          if (mesh%link_block(j) == mesh%link_block_id(n)) exit
        end do
        INSIST(n /= 0)
        link_set_mask(j) = ibset(link_set_mask(j), pos=n)
      end do
    end if
    call distribute (this%link_set_mask(:this%nlink_onP), link_set_mask)
    call gather_boundary (this%link_ip, this%link_set_mask)
    deallocate(link_set_mask)
      
    allocate(this%cfpar(this%ncell))
    call distribute (this%cfpar(:this%ncell_onP), cfpar)
    call gather_boundary (this%cell_ip, this%cfpar)
    deallocate(cfpar)

    !! Distribute the boundary face bitmask array.
    allocate(this%face_set_mask(this%nface))
    call distribute (this%face_set_mask(:this%nface_onP), face_set_mask)
    call gather_boundary (this%face_ip, this%face_set_mask)
    deallocate(face_set_mask)

    !! Distribute the node set bitmask array.
    allocate(this%node_set_mask(this%nnode))
    call distribute (this%node_set_mask(:this%nnode_onP), node_set_mask)
    call gather_boundary (this%node_ip, this%node_set_mask)
    deallocate(node_set_mask)
    
    !! Mark all link faces as boundary faces.  Note that this may leave off-P
    !! values inconsistent with the corresponding on-P value.  This is deliberatte
    !this%face_set_mask(this%lface(1,:)) = ibset(this%face_set_mask(this%lface(1,:)), pos=0)
    !this%face_set_mask(this%lface(2,:)) = ibset(this%face_set_mask(this%lface(2,:)), pos=0)

    !! Distribute the cell permutation array; gives mapping to the external cell number.
    allocate(this%xcell(this%ncell))
    call distribute (this%xcell(:this%ncell_onP), cell_perm)
    call gather_boundary (this%cell_ip, this%xcell)
    deallocate(cell_perm)

    !! Distribute the node permutation array; gives mapping to the external node number.
    allocate(this%xnode(this%nnode))
    call distribute (this%xnode(:this%nnode_onP), node_perm)
    call gather_boundary (this%node_ip, this%xnode)
    deallocate(node_perm)

    !! Broadcast the side set IDs to all processors.
    if (is_IOP) n = size(mesh%sset)
    call broadcast (n)
    allocate(this%face_set_ID(n))
    if (is_IOP) this%face_set_ID = mesh%sset%ID
    call broadcast (this%face_set_ID)

    !! Distribute the global node position array.
    allocate(this%x(3,this%nnode))
    call distribute (this%x(:,:this%nnode_onP), mesh%x)
    call gather_boundary (this%node_ip, this%x)

    !! Construct the secondary indexing arrays (THIS SUBDOMAIN MESH).
    !allocate(this%fnode(4,this%nface), this%fedge(3,this%nface), this%enode(2,this%nedge))
    allocate(this%fnode(4,this%nface))
    call assemble_face_node_list (this%cface, this%cfpar, this%cnode, this%fnode)
    !call assemble_face_edge_list (this%cface, this%cedge, this%fedge)
    !call assemble_edge_node_list (this%cedge, this%cnode, this%enode)
    
    !! Tag boundary nodes in the node set mask.
    do j = 1, this%nface
      if (btest(this%face_set_mask(j),pos=0)) then
        this%node_set_mask(this%fnode(:,j)) = ibset(this%node_set_mask(this%fnode(:,j)),pos=0)
      end if
    end do
    !! NB: I don't *think* it is necessary to do a gather_boundary on node_set_mask.
    !! The issue is whether a node could be seen to belong to a boundary face on one
    !! process but not another.  This isn't possible except *perhaps* for a node on
    !! a internal interface.  In any case, if it were possible, then a gather would
    !! not necessarily fix things, and may actually break things.  NNC, Apr 2014.

    !! Edge lengths (THIS SUBDOMAIN MESH)
    !allocate(this%length(this%nedge))
    !do j = 1, this%nedge
    !  this%length(j) = edge_length(this%x(:,this%enode(:,j)))
    !end do

    !! Directed face areas and face areas (THIS SUBDOMAIN MESH)
    allocate(this%normal(3,this%nface), this%area(this%nface))
    do j = 1, this%nface
      this%normal(:,j) = quad_face_normal(this%x(:,this%fnode(:,j)))
      this%area(j) = vector_length(this%normal(:,j))
    end do

    !! Cell and corner volumes (THIS SUBDOMAIN MESH)
    allocate(this%volume(this%ncell), this%corner_volume(8,this%ncell))
    do j = 1, this%ncell
      call eval_hex_volumes (this%x(:,this%cnode(:,j)), this%volume(j), this%corner_volume(:,j))
    end do

  end subroutine create_dist_hex_mesh

  subroutine partition_cells (np, pmeth, cnhbr, lnhbr, pass)

    use GraphModule

    integer, intent(in)    :: np            ! number of partitions
    character(len=*), intent(in) :: pmeth   ! partition method
    integer, intent(inout) :: cnhbr(:,:)    ! cell neighbor array
    integer, intent(inout) :: lnhbr(:,:)    ! link neighbor array
    integer, intent(out)   :: pass(:)       ! cell partition assignment

    integer :: i, j, k, ncell, stat, n1, n2

    type(NGraphType) :: g
    integer, allocatable :: xadj(:), adjncy(:)
    real, allocatable :: ewgt(:)
    real, parameter :: LINK_WEIGHT = 1.0

    ASSERT(np > 0)
    ASSERT(size(pass) == size(cnhbr,2))
    ASSERT(size(lnhbr,1) == 2)
    
    if (np == 1) then
      pass = 1
      return
    end if

    ncell = size(cnhbr,dim=2)

    !! Create the cell adjacency graph.
    g = CreateGraph(ncell)
    do j = 1, size(cnhbr,dim=2)
      do k = 1, size(cnhbr,dim=1)
        if (cnhbr(k,j) > j) call AddEdge (g, j, cnhbr(k,j))
      end do
    end do
    
    !! Add edges between linked cells.
    do j = 1, size(lnhbr,dim=2)
      call AddEdge (g, lnhbr(1,j), lnhbr(2,j))
    end do
    
    call GetNeighborStructure (g, xadj, adjncy)
    call DeleteGraph (g)

    !! Define edge weights.  We weight the link edges heavily to discourage
    !! them from being cut.  The weights are probably best handled as part
    !! of forming the graph, but this would require additions to GraphModule
    !! that I'm not prepared to make just yet.
    allocate(ewgt(size(adjncy)))
    ewgt = 1  ! the default
    do j = 1, size(lnhbr,2)
      n1 = lnhbr(1,j)
      n2 = lnhbr(2,j)
      !! Locate the edge (n1,n2) in the graph and set the weight.
      do i = xadj(n1), xadj(n1+1) - 1
        if (adjncy(i) == n2) exit
      end do
      ASSERT(i /= xadj(n1+1))
      ewgt(i) = LINK_WEIGHT
      !! Locate the edge (n2,n1) in the graph and set the weight.
      do i = xadj(n2), xadj(n2+1) - 1
        if (adjncy(i) == n1) exit
      end do
      ASSERT(i /= xadj(n2+1))
      ewgt(i) = LINK_WEIGHT
    end do

    select case (pmeth)
    case ('CHACO')

      !! Use Chaco to assign the cells to partitions based on the cell
      !! adjacency graph.  Some Chaco control parameters are hard-wired
      !! within the wrapper.  In Chaco, graph nodes are numbered starting
      !! at one (as here), but its C-arrays use 0-based indexing, so we
      !! need to offset the values in XADJ, which point to ADJNCY, by 1.
      call chaco_f90_wrapper2 (ncell, np, xadj-1, adjncy, ewgt, pass, stat)
      pass = pass + 1   ! we want 1-based numbering of the partitions.
      ASSERT( minval(pass) == 1 .and. maxval(pass) == np )

    case default ! shouldn't be here

      ASSERT( .false. )

    end select

    deallocate(xadj, adjncy, ewgt)

  end subroutine partition_cells

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

  subroutine organize_cells (np, pmeth, cnode, lnode, bsize, perm)

    use hexahedral_mesh_support, only: cell_neighbor_info
    use GraphModule
    use permutations

    integer, intent(in)    :: np            ! number of partitions
    character(len=*), intent(in) :: pmeth   ! partition method
    integer, intent(inout) :: cnode(:,:)    ! the cell node array
    integer, intent(inout) :: lnode(:,:)    ! the link node array
    integer, intent(out)   :: bsize(:)      ! partition block size
    integer, intent(out)   :: perm(:)       ! permutation

    integer :: i, j, k, ncell, stat, n1, n2
    integer :: pass(size(cnode,2)), cnhbr(6,size(cnode,2)), rcm_perm(size(perm))
    integer :: lnhbr(2,size(lnode,2))

    type(NGraphType) :: g
    integer, allocatable :: xadj(:), adjncy(:)
    real, allocatable :: ewgt(:)
    real, parameter :: LINK_WEIGHT = 2.0

    ASSERT( np > 0 )
    ASSERT( size(bsize) == np )
    ASSERT( size(perm) == size(cnode,2) )

    ncell = size(cnode,dim=2)

    !! Generate the cell neighbor array
    !call assemble_cell_neighbor_list (cnode, cnhbr, stat)
    !call assemble_link_neighbor_list (cnode, cnhbr, lnode, lnhbr)
    call cell_neighbor_info (cnode, lnode, cnhbr, lnhbr, stat)
    INSIST( stat == 0 )

    !! Create the cell adjacency graph
    g = CreateGraph (ncell)
    do j = 1, size(cnhbr,dim=2)
      do k = 1, size(cnhbr,dim=1)
        if (cnhbr(k,j) > j) call AddEdge (g, j, cnhbr(k,j))
      end do
    end do
    
    !! Add edges between linked cells.
    do j = 1, size(lnhbr,dim=2)
      call AddEdge (g, lnhbr(1,j), lnhbr(2,j))
    end do
    
    call GetNeighborStructure (g, xadj, adjncy)
    call DeleteGraph (g)

    if (np > 1) then  ! Assign cells to partitions.
    
      !! Define edge weights.  We weight the link edges heavily to discourage
      !! them from being cut.  The weights are probably best handled as part
      !! of forming the graph, but this would require additions to GraphModule
      !! that I'm not prepared to make just yet.
      
      allocate(ewgt(size(adjncy)))
      ewgt = 1  ! the default
      do j = 1, size(lnhbr,2)
        n1 = lnhbr(1,j)
        n2 = lnhbr(2,j)
        !! Locate the edge (n1,n2) in the graph and set the weight.
        do i = xadj(n1), xadj(n1+1) - 1
          if (adjncy(i) == n2) exit
        end do
        ASSERT(i /= xadj(n1+1))
        ewgt(i) = LINK_WEIGHT
        !! Locate the edge (n2,n1) in the graph and set the weight.
        do i = xadj(n2), xadj(n2+1) - 1
          if (adjncy(i) == n1) exit
        end do
        ASSERT(i /= xadj(n2+1))
        ewgt(i) = LINK_WEIGHT
      end do

      select case (pmeth)
      case ('CHACO')

        !! Use Chaco to assign the cells to partitions based on the cell
        !! adjacency graph.  Some Chaco control parameters are hard-wired
        !! within the wrapper.  In Chaco, graph nodes are numbered starting
        !! at one (as here), but its C-arrays use 0-based indexing, so we
        !! need to offset the values in XADJ, which point to ADJNCY, by 1.

        call chaco_f90_wrapper2 (ncell, np, xadj-1, adjncy, ewgt, pass, stat)

        pass = pass + 1   ! we want 1-based numbering of the partitions.
        ASSERT( minval(pass) == 1 .and. maxval(pass) == np )

      case default ! shouldn't be here

        ASSERT( .false. )

      end select

      call blocked_partition (pass, bsize, perm)
      deallocate(ewgt)

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
 !! PARTITION_LINKS
 !!
 !! This routine partitions the interface links.  All we required is that a
 !! link is known to the processes owning the link cell neighbors (which will
 !! also own their respective link faces); whether they are on or off-process
 !! is irrelevant because no data communication between links is needed.  Here
 !! we assign a link to the process owning the first link cell neighbor, and if
 !! the second link cell neighbor belongs to a  different partition, we add the
 !! link to the set of off-process links for that partition.  The partition
 !! block sizes are returned in the BSIZE, and the associated permutation
 !! (new-to-old) in PERM. The link indices that must be included as off-process
 !! are returned in the pointer array OFFP_INDEX (allocated by the routine):
 !! the first OFFP_SIZE(1) elements for process 1, the next OFFP_SIZE(2)
 !! elements for process 2, and so on.
 !!
 !! N.B. This routine does not reorder the LNHBR array argument.
 !!

  subroutine partition_links (lnhbr, cell_bsize, bsize, perm, offP_size, offP_index)
  
    use sets
  
    integer, intent(in)  :: lnhbr(:,:)
    integer, intent(in)  :: cell_bsize(:)
    integer, intent(out) :: bsize(:)
    integer, intent(out) :: perm(:)
    integer, intent(out) :: offP_size(:)
    integer, pointer     :: offP_index(:)
    
    integer :: j, n, offset, pass(size(perm))
    type(integer_set) :: xlink(size(bsize))
    integer, pointer :: array(:)
    
    ASSERT(size(bsize) == size(cell_bsize))
    ASSERT(size(perm) == size(lnhbr,2))
    ASSERT(minval(lnhbr) >= 1 .and. maxval(lnhbr) <= sum(cell_bsize))
    ASSERT(size(offP_size) == size(bsize))
    
    !! Assign the links to processes, either as on or off-process.
    do j = 1, size(lnhbr,dim=2)
      pass(j) = cell_part(lnhbr(1,j))
      n = cell_part(lnhbr(2,j))
      if (n /= pass(j)) call add (xlink(n), j)
    end do
    
    !! Block partition permutation (new-to-old).
    call blocked_partition (pass, bsize, perm)

    !! Copy the off-process link sets into packed array storage.
    do n = 1, size(xlink)
      offP_size(n) = size(xlink(n))
    end do
    
    n = sum(offP_size)
    allocate(offP_index(n))
    offset = 0
    do n = 1, size(xlink)
      array => to_array(xlink(n))
      offP_index(offset+1:offset+offP_size(n)) = array
      offset = offset + offP_size(n)
      deallocate(array)
    end do
    
    do n = 1, size(xlink)
      call clear (xlink(n))
    end do

  contains

    integer function cell_part (n)
      integer, intent(in) :: n
      integer :: m
      m = n
      do cell_part = 1, size(cell_bsize)
        if (m <= cell_bsize(cell_part)) return
        m = m - cell_bsize(cell_part)
      end do
    end function cell_part

  end subroutine partition_links

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

    use bitfield_type
    use mesh_importer, only: side_set

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
 !! GENERATE_NODE_SET_MASK
 !!
 !! This routine generates the node-based, bitmask vector NODE_SET_MASK. Bit n,
 !! n>0, is set if the node belongs to the nth node set.  Bit 0 will be set
 !! (later) for nodes lying on the boundary of the global domain. Note that the
 !! relationship between this implicit internal node set numbering and the user-
 !! assigned ID is stored separately.
 !!
 !! N.B.  The number of node sets is limited to 31 for the typical default
 !! 4-byte integer, but this can easily be amended by going to 8-byte integers
 !! or by using a general 'bit field' derived data type.
 !!

  subroutine generate_node_set_mask (node_set_mask, nset)

    use mesh_importer, only: node_set

    integer, intent(out) :: node_set_mask(:) ! node set bitmask array
    type(node_set), intent(in) :: nset(:) ! node set array

    integer :: i, j, k

    INSIST(size(nset) <= bit_size(node_set_mask)-1)

    node_set_mask = 0
    do k = 1, size(nset)
      ASSERT(minval(nset(k)%node) >= 1)
      ASSERT(maxval(nset(k)%node) <= size(node_set_mask))
      do j = 1, nset(k)%num_node
        i = nset(k)%node(j)
        node_set_mask(i) = ibset(node_set_mask(i),pos=k)
      end do
    end do

  end subroutine generate_node_set_mask

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

    use sets

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
          if (m /= n) call add (xcells(m), j)
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

  subroutine overlapping_cells2 (lnhbr, cell_bsize, xcells)

    use sets

    integer, intent(in) :: lnhbr(:,:)     ! link cell neighbors
    integer, intent(in) :: cell_bsize(:)  ! cell partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, p1, p2

    ASSERT( all(cell_bsize >= 0) )
    ASSERT( size(xcells) == size(cell_bsize) )
    ASSERT( size(lnhbr,1) == 2 )
    ASSERT( minval(lnhbr) >= 1 .and. maxval(lnhbr) <= sum(cell_bsize) )

    do j = 1, size(lnhbr,2)
      p1 = cell_part(lnhbr(1,j))
      p2 = cell_part(lnhbr(2,j))
      if (p1 /= p2) then
        call add (xcells(p1), lnhbr(2,j))
        call add (xcells(p2), lnhbr(1,j))
      end if
    end do

  contains

    integer function cell_part (n)
      integer, intent(in) :: n
      integer :: m
      m = n
      do cell_part = 1, size(cell_bsize)
        if (m <= cell_bsize(cell_part)) return
        m = m - cell_bsize(cell_part)
      end do
    end function cell_part

  end subroutine overlapping_cells2

end module distributed_hex_mesh
