!!
!! UNSTR_MESH_FACTORY
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"

module unstr_mesh_factory

  use kinds, only: r8
  use unstr_mesh_type
  use parameter_list_type
  implicit none
  private

  public :: unstr_mesh, new_unstr_mesh

contains

  function new_unstr_mesh (params, stat, errmsg) result (mesh)

    use ext_exodus_mesh_type
    use exodus_mesh_tools
    use exodus_mesh_io, only: read_exodus_mesh
    use permutations
    use index_partitioning
    use unstr_mesh_tools
    use parallel_communication

    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg
    type(unstr_mesh), pointer :: mesh

    integer :: j, k, n, nnode, nface, ncell
    real(r8) :: csf
    integer :: cell_bsize(nPE), node_bsize(nPE), face_bsize(nPE)
    integer, allocatable :: xcnode(:), cnode(:), xcnhbr(:), cnhbr(:), pass(:)
    integer, allocatable :: xcface(:), cface(:), cfpar(:)
    integer, allocatable :: lnhbr(:,:), lface(:,:)
    integer, allocatable :: cell_perm(:)
    integer, allocatable :: node_perm(:)
    integer, allocatable :: face_perm(:)
    integer, allocatable :: offP_size(:)
    integer, pointer :: offP_index(:)
    integer, allocatable :: p(:), ssid(:), ebid(:)
    type(ext_exodus_mesh) :: exo_mesh
    character(:), allocatable :: mesh_file

    mesh => null()

    if (is_IOP) then
      call params%get ('mesh-file', mesh_file)
      call read_exodus_mesh (mesh_file, exo_mesh, stat, errmsg)
    else
      allocate(exo_mesh%coord(3,0))
    end if
    call broadcast (stat)
    if (stat /= 0) return

    !! Scale the node coordinates
    if (is_IOP) then
      call params%get ('coord-scale-factor', csf, default=1.0d0)
      if (csf /= 1.0_r8) exo_mesh%coord = csf * exo_mesh%coord
    end if

    !! Define an empty mesh link structure
    if (is_IOP) call exo_mesh%set_no_links

    !! Create internal interfaces
    if (is_IOP) then
      if (params%is_parameter('gap-element-block-ids')) then
        call params%get ('gap-element-block-ids', ebid)
        call convert_cells_to_links (exo_mesh, ebid, stat, errmsg)
      end if
    end if
    call broadcast (stat)
    if (stat /= 0) then
      if (is_IOP) n = len(errmsg)
      call broadcast (n)
      if (.not.is_IOP) allocate(character(len=n)::errmsg)
      call broadcast (errmsg)
      return
    end if
    if (is_IOP) then
      if (params%is_parameter('interface-side-set-ids')) then
        call params%get ('interface-side-set-ids', ssid)
        call create_internal_interfaces (exo_mesh, ssid, stat, errmsg)
      end if
    end if
    call broadcast (stat)
    if (stat /= 0) then
      if (is_IOP) n = len(errmsg)
      call broadcast (n)
      if (.not.is_IOP) allocate(character(len=n)::errmsg)
      call broadcast (errmsg)
      return
    end if

    call exo_mesh%get_concat_elem_conn (xcnode, cnode)

    allocate(cell_perm(exo_mesh%num_elem))
    allocate(node_perm(exo_mesh%num_node))

    if (is_IOP) call get_cell_neighbor_array (xcnode, cnode, exo_mesh%xlnode, exo_mesh%lnode, xcnhbr, cnhbr, lnhbr, stat)
    call broadcast (stat)
    if (stat /= 0) return

    if (is_IOP) then
      ncell = exo_mesh%num_elem
      nnode = exo_mesh%num_node

      ! cell adjacency graph from cnhbr and lnhbr
      ! call graph partitioner to get part assignment
      allocate(pass(exo_mesh%num_elem))
      call partition_cells (xcnhbr, cnhbr, lnhbr, nPE, pass)
      deallocate(xcnhbr, cnhbr)

      ! generate permutation that makes it a block partition (block sizes too)
      !allocate(cell_perm(exo_mesh%num_elem))
      call blocked_partition (pass, cell_bsize, cell_perm)

      !! Reorder cell-based arrays.
      call reorder (xcnode, cnode, cell_perm)

      !! Map the values of cell-valued arrays.
      allocate(p(size(cell_perm)))
      call invert_perm (cell_perm, p)
      do j = 1, size(lnhbr,dim=2)
        lnhbr(:,j) = p(lnhbr(:,j))
      end do
      do n = 1, size(exo_mesh%sset)
        exo_mesh%sset(n)%elem = p(exo_mesh%sset(n)%elem)
      end do
      deallocate(p)

      ! partition and order the nodes
      call organize_facets (xcnode, cnode, cell_bsize, node_bsize, node_perm)

      !! Reorder node-based arrays.
      call reorder (exo_mesh%coord, node_perm)

      !! Map the values of node-valued arrays.
      allocate(p(size(node_perm)))
      call invert_perm (node_perm, p)
      do j = 1, size(exo_mesh%lnode)
        exo_mesh%lnode(j) = p(exo_mesh%lnode(j))
      end do
      do n = 1, size(exo_mesh%nset)
        exo_mesh%nset(n)%node = p(exo_mesh%nset(n)%node)
      end do
      deallocate(p)

      !! input: (xcnode,cnode), (xlnode,lnode)
      !! output: (xcface,cface), cfpar, lface, face_bsize
      !! label mesh faces
      call label_mesh_faces (xcnode, cnode, exo_mesh%xlnode, exo_mesh%lnode, nface, xcface, cface, lface)
      allocate(cfpar(exo_mesh%num_elem))
      cfpar = 0
      do j = 1, exo_mesh%num_elem
        associate (list => cface(xcface(j):xcface(j+1)-1))
          n = 0
          do k = 1, size(list)
            if (list(k) < 0) then
              n = ibset(n,k)
              list(k) = -list(k)
            end if
          end do
          cfpar(j) = n
        end associate
      end do

      !! Partition and order the faces.
      allocate(face_perm(nface))
      call organize_facets (xcface, cface, cell_bsize, face_bsize, face_perm)

      !! Map the values of face-valued arrays
      allocate(p(size(face_perm)))
      call invert_perm (face_perm, p)
      do j = 1, exo_mesh%nlink
        lface(:,j) = p(lface(:,j))
      end do
      deallocate(face_perm, p)

    else
      allocate(xcface(1), cface(0), lnhbr(2,0), lface(2,0))
      xcface(1) = 1
    end if

    allocate(mesh)

    !! Identify off-process ghost cells to include with each partition.
    call select_ghost_cells (cell_bsize, xcnode, cnode, node_bsize, &
                             xcface, cface, face_bsize, lnhbr, offP_size, offP_index)

    !! Create the cell index partition; include the off-process cells from above.
    call create (mesh%cell_ip, cell_bsize, offP_size, offP_index)
    deallocate(offP_size, offP_index)

    mesh%ncell = mesh%cell_ip%local_size()
    mesh%ncell_onP = mesh%cell_ip%onP_size()

    !! Distribute the cell permutation array; gives mapping to the external cell number.
    allocate(mesh%xcell(mesh%ncell))
    call distribute (mesh%xcell(:mesh%ncell_onP), cell_perm)
    call gather_boundary (mesh%cell_ip, mesh%xcell)
    deallocate(cell_perm)

    !! Create the node index partition and localize the global CNODE array,
    !! which identifies off-process nodes to augment the partition with.
    call init_cell_node_data (mesh, node_bsize, xcnode, cnode)

    !! Distribute the node permutation array; gives mapping to the external node number.
    allocate(mesh%xnode(mesh%nnode))
    call distribute (mesh%xnode(:mesh%nnode_onP), node_perm)
    call gather_boundary (mesh%node_ip, mesh%xnode)
    deallocate(node_perm)

    !! Create the face index partition and localize the global CFACE array,
    !! which identifies off-process faces to augment the partition with.
    if (.not.is_IOP) allocate(cfpar(0))
    call init_cell_face_data (mesh, face_bsize, xcface, cface, cfpar)
    deallocate(cfpar)

    !! Initialize the interface link data components: %NLINK, %NLINK_ONP, %LFACE, %LINK_IP, %LINK_SET_ID, %LINK_SET_MASK
    call init_link_data (mesh, exo_mesh, lnhbr, lface)
    deallocate(lnhbr, lface)

    !! Initialize the face node data: %XFNODE, %FNODE
    call init_face_node_data (mesh)

    !! %FACE_SET_ID, %FACE_SET_MASK
    call init_face_set_data (mesh, exo_mesh, xcface, cface)
    !! %NODE_SET_ID, %NODE_SET_MASK
    call init_node_set_data (mesh, exo_mesh)
    !! %CELL_SET_ID, %CELL_SET_MASK
    call init_cell_set_data (mesh, exo_mesh)

    !! MESH GEOMETRY

    allocate(mesh%x(3,mesh%nnode))
    call distribute (mesh%x(:,:mesh%nnode_onP), exo_mesh%coord)
    call gather_boundary (mesh%node_ip, mesh%x)

    allocate(mesh%volume(mesh%ncell), mesh%normal(3,mesh%nface), mesh%area(mesh%nface))
    call mesh%compute_geometry

  end function new_unstr_mesh

!!!! AUXILIARY PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partition_cells (xcnhbr, cnhbr, lnhbr, npart, part)

    use graph_type
    use graph_partitioner_factory
    use parameter_list_type

    integer, intent(in)  :: xcnhbr(:), cnhbr(:) ! cell neighbor array
    integer, intent(in)  :: lnhbr(:,:)  ! link neighbor array
    integer, intent(in)  :: npart   ! number of partitions
    integer, intent(out) :: part(:) ! cell partition assignment

    integer :: i, j, k, ncell, stat, n1, n2

    type(graph), allocatable :: g
    integer, allocatable :: xadj(:), adjncy(:)
    real, allocatable :: ewgt(:)
    real, parameter :: LINK_WEIGHT = 1.0
    class(graph_partitioner), allocatable :: gpart
    type(parameter_list) :: params

    ASSERT(npart > 0)
    ASSERT(size(part) == size(xcnhbr)-1)

    if (npart == 1) then
      part = 1
      return
    end if

    ncell = size(xcnhbr)-1

    !! Create the cell adjacency graph.
    allocate(g)
    call g%init(ncell)
    do j = 1, ncell
      associate(list => cnhbr(xcnhbr(j):xcnhbr(j+1)-1))
        do k = 1, size(list)
          if (list(k) > 0) call g%add_edge (j, list(k))
        end do
      end associate
    end do

    !! Add edges between linked cells.
    do j = 1, size(lnhbr,dim=2)
      call g%add_edge (lnhbr(1,j), lnhbr(2,j))
    end do

    call g%get_adjacency (xadj, adjncy)
    deallocate(g)

    !! Define edge weights.  We weight the link edges heavily to discourage
    !! them from being cut.  The weights are probably best handled as part
    !! of forming the graph, but this would require additions to GRAPH_TYPE
    !! that I'm not prepared to make just yet.
    allocate(ewgt(size(adjncy)))
    ewgt = 1.0  ! the default
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

    call params%set ('partitioner', 'chaco')
    call alloc_graph_partitioner (gpart, params)
    call gpart%compute (ncell, xadj, adjncy, ewgt, npart, part, stat)

  end subroutine partition_cells

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

    ASSERT(size(pass) == size(perm))
    ASSERT(minval(pass) >= 1 .and. maxval(pass) <= size(bsize))

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

    ASSERT(is_perm(perm))

  end subroutine blocked_partition

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

  subroutine organize_facets (xfacet, facet, cell_bsize, bsize, perm)

    use permutations

    integer, intent(inout) :: xfacet(:), facet(:)    ! cell facets
    integer, intent(in)    :: cell_bsize(:) ! cell partition block sizes
    integer, intent(out)   :: bsize(:)      ! facet partition block sizes
    integer, intent(out)   :: perm(:)       ! permutation

    integer :: j, n, offset, pass(size(perm)), perm1(size(perm))

    ASSERT(size(bsize) == size(cell_bsize))
    ASSERT(sum(cell_bsize) == size(xfacet)-1)
    ASSERT(size(facet) == xfacet(size(xfacet))-1)
    ASSERT(minval(facet) == 1 .and. maxval(facet) == size(perm))

    !! Generate a good ordering of the facets: Number the facets consecutively
    !! as they are encountered in FACET.  PERM1 is old-to-good numbering.
    n = 0
    perm1 = 0
    do j = 1, size(facet)
      if (perm1(facet(j)) /= 0) cycle ! already numbered
      n = n + 1
      perm1(facet(j)) = n
    end do
    ASSERT(is_perm(perm1))

    call invert_perm (perm1)  ! PERM1 is now good-to-old numbering

    !! Assign facets to partitions: a simple greedy algorithm.
    pass = 0
    offset = 1
    do n = 1, size(cell_bsize)
      associate (list => facet(xfacet(offset):xfacet(offset+cell_bsize(n))-1))
        do j = 1, size(list)
          if (pass(list(j)) == 0) pass(list(j)) = n
        end do
      end associate
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
    do j = 1, size(facet)
      facet(j) = perm1(facet(j))
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

    use integer_set_type

    integer, intent(in)  :: lnhbr(:,:)
    integer, intent(in)  :: cell_bsize(:)
    integer, intent(out) :: bsize(:)
    integer, intent(out) :: perm(:)
    integer, intent(out) :: offP_size(:)
    integer, pointer     :: offP_index(:)

    integer :: j, n, offset, pass(size(perm))
    type(integer_set) :: xlink(size(bsize))

    ASSERT(size(bsize) == size(cell_bsize))
    ASSERT(size(perm) == size(lnhbr,2))
    ASSERT(minval(lnhbr) >= 1 .and. maxval(lnhbr) <= sum(cell_bsize))
    ASSERT(size(offP_size) == size(bsize))

    !! Assign the links to processes, either as on or off-process.
    do j = 1, size(lnhbr,dim=2)
      pass(j) = cell_part(lnhbr(1,j))
      n = cell_part(lnhbr(2,j))
      if (n /= pass(j)) call xlink(n)%add (j)
    end do

    !! Block partition permutation (new-to-old).
    call blocked_partition (pass, bsize, perm)

    !! Copy the off-process link sets into packed array storage.
    do n = 1, size(xlink)
      offP_size(n) = xlink(n)%size()
    end do

    n = sum(offP_size)
    allocate(offP_index(n))
    offset = 0
    do n = 1, size(xlink)
      call xlink(n)%copy_to_array (offP_index(offset+1:))
      offset = offset + offP_size(n)
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

  subroutine select_ghost_cells (cell_bsize, xcnode, cnode, node_bsize, &
      xcface, cface, face_bsize, lnhbr, offP_size, offP_index)

    use integer_set_type
    use parallel_communication, only: is_IOP, nPE

    integer, intent(in) :: cell_bsize(:), node_bsize(:), face_bsize(:)
    integer, intent(in) :: xcnode(:), cnode(:)
    integer, intent(in) :: xcface(:), cface(:)
    integer, intent(in) :: lnhbr(:,:)
    integer, allocatable, intent(out) :: offP_size(:)
    integer, pointer, intent(out) :: offP_index(:)  ! concession to the caller

    integer :: n, offset
    type(integer_set), allocatable :: ghosts(:)

    if (is_IOP) then
      allocate(ghosts(nPE))
      call overlapping_cells (xcnode, cnode, cell_bsize, node_bsize, ghosts)
      call overlapping_cells (xcface, cface, cell_bsize, face_bsize, ghosts)
      call overlapping_cells2 (lnhbr, cell_bsize, ghosts)
      !! Copy the sets into packed array storage
      allocate(offP_size(nPE))
      do n = 1, nPE
        offP_size(n) = ghosts(n)%size()
      end do
      n = sum(offP_size)
      allocate(offP_index(n))
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

  subroutine overlapping_cells (xfacet, facet, cell_bsize, bsize, xcells)

    use integer_set_type

    integer, intent(in) :: xfacet(:), facet(:)     ! cell facets
    integer, intent(in) :: cell_bsize(:)  ! cell partition block sizes
    integer, intent(in) :: bsize(:)       ! facet partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, k, n, m, offset

    ASSERT(size(bsize) == size(cell_bsize))
    ASSERT(all(bsize >= 0))
    ASSERT(all(cell_bsize >= 0))
    ASSERT(size(xcells) == size(bsize))
    ASSERT(sum(cell_bsize) == size(xfacet)-1)
    ASSERT(size(facet) == xfacet(size(xfacet))-1)
    ASSERT(minval(facet) >= 1 .and. maxval(facet) <= sum(bsize))

    offset = 0
    do n = 1, size(cell_bsize)  ! loop over partitions
      do j = offset+1, offset+cell_bsize(n) ! loop over cells in partition N.
        associate (list => facet(xfacet(j):xfacet(j+1)-1))
          do k = 1, size(list) ! loop over the facets of cell J
            m = fpart(list(k))
            if (m /= n) call xcells(m)%add (j)
          end do
        end associate
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

    use integer_set_type

    integer, intent(in) :: lnhbr(:,:)     ! link cell neighbors
    integer, intent(in) :: cell_bsize(:)  ! cell partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, p1, p2

    ASSERT(all(cell_bsize >= 0))
    ASSERT(size(xcells) == size(cell_bsize))
    ASSERT(size(lnhbr,1) == 2)
    ASSERT(minval(lnhbr) >= 1 .and. maxval(lnhbr) <= sum(cell_bsize))

    do j = 1, size(lnhbr,2)
      p1 = cell_part(lnhbr(1,j))
      p2 = cell_part(lnhbr(2,j))
      if (p1 /= p2) then
        call xcells(p1)%add (lnhbr(2,j))
        call xcells(p2)%add (lnhbr(1,j))
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

  !! This subroutine initializes the node partition and cell node data
  !! components: NODE_IP, NNODE, NNODE_ONP, XCNODE, and CNODE.

  subroutine init_cell_node_data (this, bsize, xcnode, cnode)

    use parallel_communication, only: is_IOP
    use index_partitioning, only: create, localize_index_struct, add_offP_index

    type(unstr_mesh), intent(inout) :: this
    integer, intent(in) :: bsize(:)
    integer, intent(in) :: xcnode(:), cnode(:)

    integer :: j
    integer, allocatable :: count_g(:), count_l(:)
    integer, pointer :: offP_index(:)

    call create (this%node_ip, bsize)

    !! Translate the global indexing array into global row sizes.
    if (is_IOP) then
      count_g = xcnode(2:) - xcnode(:size(xcnode)-1)
    else
      allocate(count_g(0))
    end if

    call localize_index_struct (count_g, cnode, this%cell_ip, this%node_ip, count_l, this%cnode, offP_index)
    call add_offP_index (this%node_ip, offP_index)
    deallocate(count_g, offP_index)

    !! Translate the local row sizes into the local indexing array.
    allocate(this%xcnode(1+size(count_l)))
    this%xcnode(1) = 1
    do j = 1, size(count_l)
      this%xcnode(j+1) = this%xcnode(j) + count_l(j)
    end do
    deallocate(count_l)

    this%nnode = this%node_ip%local_size()
    this%nnode_onP = this%node_ip%onP_size()

  end subroutine init_cell_node_data

  !! This subroutine initializes the face partition and cell face data
  !! components: FACE_IP, NFACE, NFACE_ONP, XCFACE, CFACE, and CFPAR.

  subroutine init_cell_face_data (this, bsize, xcface, cface, cfpar)

    use parallel_communication, only: is_IOP, distribute
    use index_partitioning, only: create, localize_index_struct, add_offP_index, gather_boundary

    type(unstr_mesh), intent(inout) :: this
    integer, intent(in) :: bsize(:)
    integer, intent(in) :: xcface(:), cface(:)
    integer, intent(in) :: cfpar(:)

    integer :: j
    integer, allocatable :: count_g(:), count_l(:)
    integer, pointer :: offP_index(:)

    call create (this%face_ip, bsize)

    !! Translate the global indexing array into global row sizes.
    if (is_IOP) then
      count_g = xcface(2:) - xcface(:size(xcface)-1)
    else
      allocate(count_g(0))
    end if

    call localize_index_struct (count_g, cface, this%cell_ip, this%face_ip, count_l, this%cface, offP_index)
    call add_offP_index (this%face_ip, offP_index)
    deallocate(count_g, offP_index)

    !! Translate the local row sizes into the local indexing array.
    allocate(this%xcface(1+size(count_l)))
    this%xcface(1) = 1
    do j = 1, size(count_l)
      this%xcface(j+1) = this%xcface(j) + count_l(j)
    end do
    deallocate(count_l)

    !! Distribute the CFPAR mask array.
    allocate(this%cfpar(this%ncell))
    call distribute (this%cfpar(:this%ncell_onP), cfpar)
    call gather_boundary (this%cell_ip, this%cfpar)

    this%nface = this%face_ip%local_size()
    this%nface_onP = this%face_ip%onP_size()

  end subroutine init_cell_face_data

  !! This subroutine initializes the face-node connectivity data components,
  !! stored in the pair of arrays XFNODE and FNODE.  The face-node lists are
  !! normalized (the smallest node index begins the list) and oriented
  !! consistently with the CFPAR data.

  subroutine init_face_node_data (this)

    use cell_topology, only: get_cell_face_sizes, get_face_nodes

    type(unstr_mesh), intent(inout) :: this

    integer :: j, k, n
    integer, allocatable :: fsize(:), list(:)

    ASSERT(size(this%xcnode) == size(this%xcface))
    ASSERT(size(this%cfpar) == size(this%xcface)-1)
    ASSERT(size(this%cnode) == this%xcnode(size(this%xcnode))-1)
    ASSERT(size(this%cface) == this%xcface(size(this%xcface))-1)
    ASSERT(minval(this%cnode) > 0)
    ASSERT(minval(this%cface) > 0)

    !! Determine storage requirements; XFNODE(k) to store num nodes on face k.
    allocate(fsize(this%nface))
    fsize = 0
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1), &
                 cell_faces => this%cface(this%xcface(j):this%xcface(j+1)-1))
        call get_cell_face_sizes (cell_nodes, list)
        where (fsize(cell_faces) == 0) fsize(cell_faces) = list
      end associate
    end do

    !! Generate the XFNODE indexing array; exclusive prefix sum of FSIZE.
    allocate(this%xfnode(this%nface+1))
    this%xfnode(1) = 1
    do j = 1, this%nface
      this%xfnode(j+1) = this%xfnode(j) + fsize(j)
    end do
    deallocate(fsize)

    !! Fill the FNODE array.
    allocate(this%fnode(this%xfnode(this%nface+1)-1))
    this%fnode = 0
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1), &
                 cell_faces => this%cface(this%xcface(j):this%xcface(j+1)-1))
        do k = 1, size(cell_faces)
          n = cell_faces(k)
          associate (face_nodes => this%fnode(this%xfnode(n):this%xfnode(n+1)-1))
            if (face_nodes(1) == 0) then
              call get_face_nodes (cell_nodes, k, list, normalize=.true., &
                                   reverse=btest(this%cfpar(j),pos=k))
              face_nodes = list
            end if
          end associate
        end do
      end associate
    end do
    ASSERT(minval(this%fnode) > 0)

  end subroutine init_face_node_data

  !! This subroutine initializes the interface link data.

  subroutine init_link_data (this, exo_mesh, lnhbr, lface)

    use parallel_communication, only: is_IOP, nPE, collate, broadcast, distribute
    use permutations, only: reorder, invert_perm
    use index_partitioning, only: create, localize_index_array, gather_boundary
    use bitfield_type
    use ext_exodus_mesh_type

    type(unstr_mesh), intent(inout) :: this
    type(ext_exodus_mesh), intent(inout) :: exo_mesh
    integer, intent(inout) :: lnhbr(:,:), lface(:,:)

    integer :: j, n
    integer, allocatable :: offP_size(:), cell_bsize(:), bsize(:), perm(:)
    integer, pointer :: offP_index(:) => null()
    type(bitfield), allocatable :: link_set_mask(:)

    !! Partition the links
    allocate(cell_bsize(merge(nPE,0,is_IOP)))
    call collate (cell_bsize, this%cell_ip%onP_size())
    if (is_IOP) then
      allocate(bsize(nPE), perm(exo_mesh%nlink), offP_size(nPE))
      call partition_links (lnhbr, cell_bsize, bsize, perm, offP_size, offP_index)
      call reorder (lface, perm)
      call reorder (exo_mesh%link_block, perm)
      call invert_perm (perm)
      offP_index = perm(offP_index)
      deallocate(perm)
    else
      allocate(bsize(0), offP_size(0), offP_index(0))
    end if

    !! THIS%LINK_IP: create the link index partition; include off-process links from above.
    call create (this%link_ip, bsize, offP_size, offP_index)
    deallocate(bsize, offP_size, offP_index)

    !! THIS%LFACE: distribute and localize the link face indexing array.
    call localize_index_array (lface, this%link_ip, this%face_ip, this%lface, offP_index)
    INSIST(size(offP_index) == 0)
    deallocate(offP_index)

    this%nlink = this%link_ip%local_size()
    this%nlink_onP = this%link_ip%onP_size()

    !! THIS%LINK_SET_ID
    n = exo_mesh%nlblock
    call broadcast (n)
    allocate(this%link_set_id(n))
    if (is_IOP) this%link_set_id = exo_mesh%link_block_id
    call broadcast (this%link_set_id)

    !! Convert the array of link block IDs to a bit mask.
    allocate(link_set_mask(merge(exo_mesh%nlink,0,is_IOP)))
    if (is_IOP) then
      INSIST(n <= bit_size(link_set_mask)-1)
      link_set_mask = ZERO_BITFIELD
      do j = 1, exo_mesh%nlink
        do n = size(exo_mesh%link_block_id), 1, -1
          if (exo_mesh%link_block(j) == exo_mesh%link_block_id(n)) exit
        end do
        INSIST(n /= 0)
        link_set_mask(j) = ibset(link_set_mask(j), pos=n)
      end do
    end if

    !! THIS%LINK_SET_MASK
    allocate(this%link_set_mask(this%nlink))
    call distribute (this%link_set_mask(:this%nlink_onP), link_set_mask)
    call gather_boundary (this%link_ip, this%link_set_mask)

  end subroutine init_link_data

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

    type(unstr_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: exo_mesh

    integer :: i, j, n, nnode_tot
    integer, allocatable :: node_set_mask(:)
    logical, allocatable :: bnode(:)

    ASSERT(allocated(this%xfnode))
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
      if (btest(this%face_set_mask(j),pos=0)) then  ! boundary face
        associate (face_nodes => this%fnode(this%xfnode(j):this%xfnode(j+1)-1))
          bnode(face_nodes) = .true.
        end associate
      end if
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

    type(unstr_mesh), intent(inout) :: this
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

  subroutine init_face_set_data (this, exo_mesh, xcface, cface)

    use bitfield_type
    use exodus_mesh_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary

    type(unstr_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: exo_mesh
    integer, intent(in) :: xcface(:), cface(:)

    integer :: i, j, n, nface_tot
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
        ASSERT(minval(exo_mesh%sset(n)%elem) >= 1)
        ASSERT(maxval(exo_mesh%sset(n)%elem) <= size(xcface)-1)
        ASSERT(minval(exo_mesh%sset(n)%face) >= 1)
        do i = 1, exo_mesh%sset(n)%num_side
          j = xcface(exo_mesh%sset(n)%elem(i)) + exo_mesh%sset(n)%face(i) - 1 ! index of side in CFACE
          ASSERT(j < xcface(exo_mesh%sset(n)%elem(i)+1))
          face_set_mask(cface(j)) = ibset(face_set_mask(cface(j)), pos=n)
        end do
      end do

      !! Count the references to each face.  A count of 1 indicates a boundary
      !! face, 2 an interior face.  Anything else exposes a bad mesh topology.
      allocate(tag(size(face_set_mask)))
      tag = 0
      do j = 1, size(cface)
        tag(cface(j)) = tag(cface(j)) + 1
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

end module unstr_mesh_factory
