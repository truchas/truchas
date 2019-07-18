!!
!! UNSTR_MESH_FACTORY
!!
!! Provides a procedure for instantiating a new UNSTR_MESH object that stores
!! a distributed unstructured mixed-element mesh.  Information about the mesh
!! is passed using a PARAMETER_LIST object.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  NEW_UNSTR_MESH(PARAMS, STAT, ERRMSG) returns a pointer to a newly allocated
!!    UNSTR_MESH object that has been initialized according to the information
!!    specified by the parameter list PARAMS.  If an error is encountered, the
!!    integer STAT is assigned a non-zero value and the alloctable deferred-
!!    length character string ERRMSG assigned an explanatory message.  The
!!    following parameter list parameters are recognized:
!!
!!    'mesh-file' -- path to an ExodusII mesh file (required)
!!    'coord-scale-factor' -- a multiplicative scaling factor applied to the
!!        node coordinates (optional)
!!    'gap-element-block-ids' -- a list of Exodus element block IDs whose
!!        elements should be converted into interface links. The elements in the
!!        specified element blocks must satisfy certain constraints. (optional)
!!    'interface-side-set-ids' -- a list of Exodus side set IDs defining an
!!        internal mesh interface along which the mesh is sliced open and
!!        stitched together with interface links. (optional)
!!    'exodus-block-modulus' -- Immediately after the Exodus mesh is read the
!!        element block IDs are overwritten with their values modulo the value
!!        of this parameter.  This is optional; 0 is treated as if the parameter
!!        was not specified.
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

  function new_unstr_mesh (params, stat, errmsg) result (this)

    use ext_exodus_mesh_type
    use exodus_mesh_tools
    use exodus_mesh_io, only: read_exodus_mesh
    use permutations
    use simple_partitioning_methods, only: get_block_partition, read_partition
    use index_partitioning
    use unstr_mesh_tools
    use parallel_communication
    use truchas_logging_services
    use string_utilities, only: i_to_c

    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg
    type(unstr_mesh), pointer :: this

    integer :: j, k, n, nnode, nface, ncell, new_id, exodus_block_modulus, pfirst
    integer :: cell_psize(nPE), node_psize(nPE), face_psize(nPE)
    integer, allocatable :: xcnode(:), cnode(:), xcnhbr(:), cnhbr(:), part(:)
    integer, allocatable :: xndcell(:), ndcell(:)
    integer, allocatable :: xcface(:), cface(:), cfpar(:), lnhbr(:,:), lface(:,:)
    integer, allocatable :: cell_perm(:), node_perm(:), offP_size(:), offP_index(:)
    integer, allocatable :: perm(:), ssid(:), ebid(:)
    type(ext_exodus_mesh) :: mesh
    character(:), allocatable :: mesh_file, msg, string
    real(r8) :: csf
    logical :: have_parent_node ! TEMPORARY

    this => null()

    !! Read the Exodus mesh file into MESH.
    if (is_IOP) then
      call params%get ('mesh-file', mesh_file)
      call TLS_info ('  reading ExodusII mesh file "' // mesh_file // '"')
      call read_exodus_mesh (mesh_file, mesh, stat, errmsg)
      if (stat /= 0) errmsg = 'error reading mesh file: ' // errmsg
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    !! Overwrite the Exodus block IDs with their congruent values.
    if (is_IOP) then
      call params%get ('exodus-block-modulus', exodus_block_modulus, default=0)
      if (exodus_block_modulus > 0) then
        do n = 1, mesh%num_eblk
          associate (id => mesh%eblk(n)%id)
            new_id = modulo(id, exodus_block_modulus)
            if (new_id /= id) then
              msg = '  element block ' // i_to_c(id) // ' merged with block ' // i_to_c(new_id)
              call TLS_info (msg, TLS_VERB_NORMAL)
              id = new_id
            end if
          end associate
        end do
      end if
    end if

    !! Define an empty mesh link structure
    if (is_IOP) call mesh%set_no_links

    !! Create internal interfaces from gap element blocks.
    if (is_IOP) then
      if (params%is_parameter('gap-element-block-ids')) then
        call TLS_info ('  processing gap element blocks', TLS_VERB_NORMAL)
        call params%get ('gap-element-block-ids', ebid)
        call convert_cells_to_links (mesh, ebid, stat, errmsg)
        if (stat /= 0) errmsg = 'error processing gap element blocks: ' // errmsg
      else
        stat = 0
      end if
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    !! Create internal interfaces from side sets.
    if (is_IOP) then
      if (params%is_parameter('interface-side-set-ids')) then
        call TLS_info ('  processing interface side sets', TLS_VERB_NORMAL)
        call params%get ('interface-side-set-ids', ssid)
        call create_internal_interfaces (mesh, ssid, stat, errmsg)
        if (stat /= 0) errmsg = 'error processing interface side sets: ' // errmsg
      else
        stat = 0
      end if
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    ncell = mesh%num_elem
    nnode = mesh%num_node

    !! Flatten the Exodus element block structure.
    call mesh%get_concat_elem_conn (xcnode, cnode)

    !! Generate the cell and link neighbor arrays.
    if (is_IOP) then
      call TLS_info ('  finding cell neighbors', TLS_VERB_NORMAL)
      call get_cell_neighbor_array (xcnode, cnode, mesh%xlnode, mesh%lnode, xcnhbr, cnhbr, lnhbr, stat)
      if (stat /= 0) errmsg = 'get_cell_neighbor_array: invalid mesh topology detected'
    else
      allocate(xcnhbr(1), cnhbr(0), lnhbr(2,0))
      xcnhbr(1) = 1
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    !! Partition and order the cells.
    allocate(cell_perm(ncell))
    if (is_IOP) then
      call TLS_info ('  partitioning the mesh cells', TLS_VERB_NORMAL)
      !! Partition the cell neighbor graph.
      allocate(part(mesh%num_elem))
      call params%get ('partitioner', string, default='chaco')
      if (nPE == 1) then
        part = 1
        stat = 0
      else if (string == 'block') then
        call get_block_partition (nPE, part)
        stat = 0
      else if (string == 'file') then
        call params%get ('partition-file', string)
        call params%get ('first-partition', pfirst, default=0)
        call read_partition (string, pfirst, nPE, part, stat, errmsg)
        if (stat /= 0) errmsg = 'error reading cell partition: ' // errmsg
      else
        call partition_cells (params, xcnhbr, cnhbr, lnhbr, nPE, part, stat, errmsg)
        if (stat /= 0) errmsg = 'error computing cell partition: ' // errmsg
      end if
    end if
    call broadcast_status (stat, errmsg)
    if (stat /= 0) return

    if (is_IOP) then
      !! Compute the partition sizes and the permutation making this a block partition.
      call blocked_partition (part, cell_psize, cell_perm)
      !! Reorder cell-based arrays.
      call reorder (xcnode, cnode, cell_perm)
      call reorder (xcnhbr, cnhbr, cell_perm)
      !! Map the values of cell-valued arrays.
      allocate(perm(size(cell_perm)))
      call invert_perm (cell_perm, perm)
      do j = 1, size(cnhbr)
        if (cnhbr(j) > 0) cnhbr(j) = perm(cnhbr(j))
      end do
      do j = 1, size(lnhbr,dim=2)
        lnhbr(:,j) = perm(lnhbr(:,j))
      end do
      do n = 1, size(mesh%sset)
        mesh%sset(n)%elem = perm(mesh%sset(n)%elem)
      end do
      deallocate(perm)
    end if

    !! Partition and order the nodes.
    allocate(node_perm(nnode))
    if (is_IOP) then
      call TLS_info ('  partitioning the mesh nodes', TLS_VERB_NORMAL)
      call partition_facets (xcnode, cnode, cell_psize, node_psize, node_perm)
      !! Reorder node-based arrays.
      call reorder (mesh%coord, node_perm)
      if (allocated(mesh%parent_node)) call reorder(mesh%parent_node, node_perm) ! TEMPORARY
      !! Map the values of node-valued arrays.
      allocate(perm(size(node_perm)))
      call invert_perm (node_perm, perm)
      do j = 1, size(cnode)
        cnode(j) = perm(cnode(j))
      end do
      do j = 1, size(mesh%lnode)
        mesh%lnode(j) = perm(mesh%lnode(j))
      end do
      do n = 1, size(mesh%nset)
        mesh%nset(n)%node = perm(mesh%nset(n)%node)
      end do
      if (allocated(mesh%parent_node)) then ! TEMPORARY
        do j = 1, size(mesh%parent_node)
          mesh%parent_node(j) = perm(mesh%parent_node(j))
        end do
      end if
      deallocate(perm)
    end if

    !! Enumerate and partition the mesh faces.
    allocate(cfpar(ncell))
    if (is_IOP) then
      call TLS_info ('  numbering the mesh faces', TLS_VERB_NORMAL)
      call label_mesh_faces (xcnode, cnode, mesh%xlnode, mesh%lnode, nface, xcface, cface, lface)
      !! Extract the relative face orientation info.
      cfpar = 0
      do j = 1, mesh%num_elem
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
      call TLS_info ('  partitioning the mesh faces', TLS_VERB_NORMAL)
      allocate(perm(nface))
      call partition_facets (xcface, cface, cell_psize, face_psize, perm)
      call invert_perm (perm)
      do j = 1, size(cface)
        cface(j) = perm(cface(j))
      end do
      do j = 1, mesh%nlink
        lface(:,j) = perm(lface(:,j))
      end do
      deallocate(perm)
    else
      allocate(xcface(1), cface(0), lface(2,0))
      xcface(1) = 1
    end if

    !! Identify off-process ghost cells to include with each partition.
    call TLS_info ('  identifying off-process ghost cells', TLS_VERB_NORMAL)
    call select_ghost_cells (cell_psize, xcnhbr, cnhbr, xcnode, cnode, node_psize, &
                             xcface, cface, face_psize, lnhbr, offP_size, offP_index)
    deallocate(xcnhbr, cnhbr)

    !! Begin initializing the UNSTR_MESH result object.
    call TLS_info ('  generating parallel mesh structure')
    allocate(this)

    !! Create the cell index partition; include the off-process cells from above.
    call this%cell_ip%init (cell_psize, offP_size, offP_index)
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
    call init_cell_node_data (this, node_psize, xcnode, cnode)

    !! Distribute the node permutation array; gives mapping to the external node number.
    allocate(this%xnode(this%nnode))
    call distribute (this%xnode(:this%nnode_onP), node_perm)
    call gather_boundary (this%node_ip, this%xnode)
    deallocate(node_perm)

    !! Distribute the node parent array (TEMPORARY)
    have_parent_node = allocated(mesh%parent_node)
    call broadcast (have_parent_node)
    if (have_parent_node) then
      allocate(this%parent_node(this%nnode))
      if (.not.is_IOP) allocate(mesh%parent_node(0))
      call distribute (this%parent_node(:this%nnode_onP), mesh%parent_node)
      call gather_boundary (this%node_ip, this%parent_node)
    end if

    !! Create the face index partition and localize the global CFACE array,
    !! which identifies off-process faces to augment the partition with.
    call init_cell_face_data (this, face_psize, xcface, cface, cfpar)
    deallocate(cfpar)

    !! Initialize the interface link data components.
    call init_link_data (this, mesh, lnhbr, lface)
    deallocate(lnhbr, lface)

    !! Initialize the secondary face-node indexing array.
    call init_face_node_data (this)
    call init_face_cell_data (this)

    !! Generate the cell neighbor data for each subdomain.
    call get_cell_neighbor_array (this%xcnode, this%cnode, this%xcnhbr, this%cnhbr, stat)
    INSIST(stat == 0)

    !! Initialize the node, face, and cell set data.
    call init_face_set_data (this, mesh, xcface, cface)
    call init_node_set_data (this, mesh)
    call init_cell_set_data (this, mesh)
    deallocate(xcface, cface)

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

    ! Generate node->cell connectivity in staggered array (FOR LOCAL DOMAIN 1:_onP ONLY)
    call init_node_to_cell_data(this)
    call init_cell_to_cell_through_node_data(this)

    !! Initialize the mesh geometry data components.
    allocate(this%volume(this%ncell), this%normal(3,this%nface), this%area(this%nface))
    call this%compute_geometry

  end function new_unstr_mesh

!!!! AUXILIARY PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partition_cells (params, xcnhbr, cnhbr, lnhbr, npart, part, stat, errmsg)

    use graph_type
    use graph_partitioner_factory

    type(parameter_list) :: params
    integer, intent(in)  :: xcnhbr(:), cnhbr(:) ! cell neighbor array
    integer, intent(in)  :: lnhbr(:,:)  ! link neighbor array
    integer, intent(in)  :: npart   ! number of partitions
    integer, intent(out) :: part(:) ! cell partition assignment
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, k, ncell, n1, n2

    type(graph), allocatable :: g
    integer, allocatable :: xadj(:), adjncy(:)
    real, allocatable :: ewgt(:)
    real, parameter :: LINK_WEIGHT = 1.0
    class(graph_partitioner), allocatable :: gpart

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

    call alloc_graph_partitioner (gpart, params)
    call gpart%compute (ncell, xadj, adjncy, ewgt, npart, part, stat, errmsg)

  end subroutine partition_cells

  !! This auxiliary subroutine partitions the mesh facets of one type (nodes,
  !! edges, or faces) and generates a new numbering of them for which the
  !! partition becomes a block partition and the facets within partitions are
  !! well-ordered.  This is a serial procedure that operates on facet data for
  !! the entire mesh.  FACET(XFACET(j):XFACET(j+1)-1) are the facets of one
  !! type belonging to cell j.  It is assumed that the cells have already been
  !! partitioned and renumbered, and that XFACET and FACET reflects that new
  !! numbering.  PSIZE returns the number of facets in each partition, and
  !! the new facet numbering is returned in the permutation array PERM, which
  !! maps new facet numbers to old.
  !!
  !! NB: The caller is responsible for applying the renumbering.  That means
  !! permuting any facet-based arrays, and mapping the values of all facet-
  !! valued arrays, including FACET itself.
  !!
  !! The facet partition is based on the cell partition given by the cell
  !! partition sizes CELL_PSIZE.  A facet is assigned to the partition of one
  !! of the cells it is adjacent to.  We use a simple greedy algorithm: for
  !! n=1,2,... assign any unassigned facet adjacent to a cell in partition n
  !! to partition n.  NB: While this may lead to a poorly balanced partition
  !! for facets, it is not clear that this significantly impacts performance.
  !!
  !! To order the facets in a partition we rely on the well-orderedness of the
  !! cells.  We merely number them as they are encountered in the FACET array.

  subroutine partition_facets (xfacet, facet, cell_psize, psize, perm)

    use permutations

    integer, intent(in)  :: xfacet(:), facet(:) ! cell facets
    integer, intent(in)  :: cell_psize(:)       ! cell partition block sizes
    integer, intent(out) :: psize(:)            ! facet partition block sizes
    integer, intent(out) :: perm(:)             ! permutation

    integer :: j, n, offset, part(size(perm)), perm1(size(perm))

    ASSERT(size(psize) == size(cell_psize))
    ASSERT(sum(cell_psize) == size(xfacet)-1)
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
    part = 0
    offset = 1
    do n = 1, size(cell_psize)
      associate (list => facet(xfacet(offset):xfacet(offset+cell_psize(n))-1))
        do j = 1, size(list)
          if (part(list(j)) == 0) part(list(j)) = n
        end do
      end associate
      offset = offset + cell_psize(n)
    end do

    !! Partition assignment relative to the good ordering.
    call reorder (part, perm1)

    !! Block partition permutation (new-to-good).
    call blocked_partition (part, psize, perm)

    !! Total facet permutation (new-to-old).
    do j = 1, size(perm)
      perm(j) = perm1(perm(j))
    end do

  end subroutine partition_facets

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
 !! block sizes are returned in the PSIZE, and the associated permutation
 !! (new-to-old) in PERM. The link indices that must be included as off-process
 !! are returned in the pointer array OFFP_INDEX (allocated by the routine):
 !! the first OFFP_SIZE(1) elements for process 1, the next OFFP_SIZE(2)
 !! elements for process 2, and so on.
 !!
 !! N.B. This routine does not reorder the LNHBR array argument.
 !!

  subroutine partition_links (lnhbr, cell_psize, xlnode, lnode, node_psize, psize, perm, offP_size, offP_index)

    use integer_set_type

    integer, intent(in)  :: lnhbr(:,:)
    integer, intent(in)  :: cell_psize(:)
    integer, intent(in)  :: xlnode(:), lnode(:)
    integer, intent(in)  :: node_psize(:)
    integer, intent(out) :: psize(:)
    integer, intent(out) :: perm(:)
    integer, intent(out) :: offP_size(:)
    integer, allocatable, intent(out) :: offP_index(:)

    integer :: j, k, n, offset, pass(size(perm))
    type(integer_set) :: xlink(size(psize))

    ASSERT(size(psize) == size(cell_psize))
    ASSERT(size(perm) == size(lnhbr,2))
    ASSERT(minval(lnhbr) >= 1 .and. maxval(lnhbr) <= sum(cell_psize))
    ASSERT(size(offP_size) == size(psize))

    !! Assign the links to processes, either as on or off-process.
    do j = 1, size(lnhbr,dim=2)
      pass(j) = cell_part(lnhbr(1,j))
      n = cell_part(lnhbr(2,j))
      if (n /= pass(j)) call xlink(n)%add (j)
!      do k = xlnode(j), xlnode(j+1)-1
!        n = node_part(lnode(k))
!        if (n /= pass(j)) call xlink(n)%add (j)
!      end do
    end do

    !! Block partition permutation (new-to-old).
    call blocked_partition (pass, psize, perm)

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
      do cell_part = 1, size(cell_psize)
        if (m <= cell_psize(cell_part)) return
        m = m - cell_psize(cell_part)
      end do
    end function cell_part

    integer function node_part (n)
      integer, intent(in) :: n
      integer :: m
      m = n
      do node_part = 1, size(node_psize)
        if (m <= node_psize(node_part)) return
        m = m - node_psize(node_part)
      end do
    end function node_part

  end subroutine partition_links

  !! Given a partition assignment array PART, this auxiliary subroutine
  !! computes the permutation array PERM that makes the partition a block
  !! partition.  The partition sizes are returned in PSIZE.  PERM maps the
  !! new numbering to the original, and preserves the relative order of
  !! elements in the same partition.  Specifically, 1) PART(PERM(:)) is
  !! sorted (non-decreasing), and 2) j < k whenever PERM(j) < PERM(k) and
  !! PART(PERM(j)) = PART(PERM(k)).

  subroutine blocked_partition (part, psize, perm)

    use permutations

    integer, intent(in)  :: part(:)   ! partition assignment
    integer, intent(out) :: psize(:)  ! partition block size
    integer, intent(out) :: perm(:)   ! permutation

    integer :: j, n, next(size(psize))

    ASSERT(size(part) == size(perm))
    ASSERT(minval(part) >= 1 .and. maxval(part) <= size(psize))

    !! Compute the block size of each partition.
    psize = 0
    do j = 1, size(part)
      psize(part(j)) = psize(part(j)) + 1
    end do

    !! NEXT(j) is the next free cell number for partition j.
    next(1) = 1
    do n = 2, size(psize)
      next(n) = next(n-1) + psize(n-1)
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
  !! contains a face or node that belongs to a subdomain and the cell itself
  !! does not belong to the subdomain, then it is added to the subdomain as
  !! a ghost cell.  This ensures that every face and node belonging to the
  !! subdomain will have complete cell support.  The motivation for this
  !! stems from the finite element context in which the equation for a DoF
  !! located at a facet depends on calculations over all the cells that
  !! contain the facet.  An additional criterion is if a cell adjacent to
  !! an interface belongs to a partition then the cell it is linked to across
  !! the interface is added as a ghost if it does not also belong to the
  !! partition.  This will guarantee that the associated pair of interface
  !! faces is present in the subdomain so that an interface condition can be
  !! fully formed.  This is a serial procedure.

  subroutine select_ghost_cells (cell_psize, xcnhbr, cnhbr, xcnode, cnode, node_psize, &
      xcface, cface, face_psize, lnhbr, offP_size, offP_index)

    use integer_set_type
    use parallel_communication, only: is_IOP, nPE

    integer, intent(in) :: cell_psize(:), node_psize(:), face_psize(:)
    integer, intent(in) :: xcnhbr(:), cnhbr(:)
    integer, intent(in) :: xcnode(:), cnode(:)
    integer, intent(in) :: xcface(:), cface(:)
    integer, intent(in) :: lnhbr(:,:)
    integer, allocatable, intent(out) :: offP_size(:), offP_index(:)

    integer :: n, offset
    type(integer_set), allocatable :: ghosts(:)

    if (is_IOP) then
      allocate(ghosts(nPE))
      call all_cell_neighbors (xcnode, cnode, cell_psize, ghosts)
      !NNC, Feb 2016.  The above call produces the same ghosts below, and more.
      !call overlapping_cells (xcnhbr, cnhbr, cell_psize, cell_psize, ghosts)
      !call overlapping_cells (xcnode, cnode, cell_psize, node_psize, ghosts)
      !call overlapping_cells (xcface, cface, cell_psize, face_psize, ghosts)
      call overlapping_cells2 (lnhbr, cell_psize, ghosts)
      !! Copy the sets into packed array storage
      offP_size = ghosts%size()
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

  !! This auxiliary subroutine identifies, for each partition, the entire
  !! layer of first neighbor cells from other partitions.  These are cells
  !! that belong to other partitions but which share a node with a cell in
  !! the partition.  This is driven by the least squares operator used by
  !! fluid flow, which uses a cell-based stencil that consists of a cell
  !! and all cells it shares a node with.  The solid mechanics node-node
  !! connectivity also requires this complete layer of ghost cells.  This
  !! subsumes the more refined criteria effected by OVERLAPPING_CELLS.

  subroutine all_cell_neighbors (xcnode, cnode, cell_psize, xcells)

    use integer_set_type

    integer, intent(in) :: xcnode(:), cnode(:)
    integer, intent(in) :: cell_psize(:)
    type(integer_set), intent(inout) :: xcells(:)

    integer :: i, j, k, n, offset, jlower, jupper
    integer, allocatable :: nhbr(:)
    type(integer_set) :: nsupp(maxval(cnode))

    ASSERT(all(cell_psize >= 0))
    ASSERT(size(xcells) == size(cell_psize))
    ASSERT(sum(cell_psize) == size(xcnode)-1)
    ASSERT(size(cnode) == xcnode(size(xcnode))-1)
    ASSERT(minval(cnode) >= 0)

    !! For each node, generate the set of cells that contain it.
    do j = 1, size(xcnode)-1
      associate(jnode => cnode(xcnode(j):xcnode(j+1)-1))
        do k = 1, size(jnode)
          call nsupp(jnode(k))%add(j)
        end do
      end associate
    end do

    !! For each cell, scan all the cells that contain one of its nodes, and add
    !! those that belong to another partition to the ghost set of the partition.
    offset = 0
    do n = 1, size(cell_psize)
      jlower = offset + 1
      jupper = offset + cell_psize(n)
      do j = jlower, jupper ! loop over cells in partition N
        associate (jnode => cnode(xcnode(j):xcnode(j+1)-1))
          do k = 1, size(jnode) ! for each node of cell J
            nhbr = nsupp(jnode(k))
            do i = 1, size(nhbr) ! loop over cells containing the node
              if (nhbr(i) < jlower .or. nhbr(i) > jupper) call xcells(n)%add(nhbr(i))
            end do
          end do
        end associate
      end do
      offset = offset + cell_psize(n)
    end do

  end subroutine all_cell_neighbors

  !! This auxiliary subroutine identifies for each partition those cells that
  !! contain a facet in the partition but that themselves belong to a different
  !! partition.  For a given partition the identitified cell set will be the
  !! minimal one that if added to the cells already belonging to the partition
  !! will ensure that every facet belonging to the partition has complete
  !! cell-support.  The motivation for this stems from the finite element
  !! context in which the equation for a DoF located at a facet depends on
  !! calculations over all the cells that contain that facet.
  !!
  !! FACET(XFACET(j):XFACET(j+1)-1) are the indices of the facets of one type
  !! (node, edge, or face) contained in cell j.  The cell partitioning is
  !! described by the cell partition block sizes CELL_PSIZE, and the facet
  !! partitioning by their block sizes PSIZE.  The cells identified for
  !! partition n are added to the set XCELLS(n).  Thus this routine can be
  !! called for different FACET arrays, each call adding to XCELLS the cells
  !! for that type of facet.  This is a serial procedure.

  subroutine overlapping_cells (xfacet, facet, cell_psize, psize, xcells)

    use integer_set_type

    integer, intent(in) :: xfacet(:), facet(:)  ! cell facets
    integer, intent(in) :: cell_psize(:)        ! cell partition block sizes
    integer, intent(in) :: psize(:)             ! facet partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, k, n, m, offset

    ASSERT(size(psize) == size(cell_psize))
    ASSERT(all(psize >= 0))
    ASSERT(all(cell_psize >= 0))
    ASSERT(size(xcells) == size(psize))
    ASSERT(sum(cell_psize) == size(xfacet)-1)
    ASSERT(size(facet) == xfacet(size(xfacet))-1)
    ASSERT(minval(facet) >= 0 .and. maxval(facet) <= sum(psize))

    offset = 0
    do n = 1, size(cell_psize)  ! loop over partitions
      do j = offset+1, offset+cell_psize(n) ! loop over cells in partition N.
        associate (list => facet(xfacet(j):xfacet(j+1)-1))
          do k = 1, size(list) ! loop over the facets of cell J
            if (list(k) == 0) cycle
            m = fpart(list(k))
            if (m /= n) call xcells(m)%add (j)
          end do
        end associate
      end do
      offset = offset + cell_psize(n)
    end do

  contains

    integer function fpart (n)
      integer, intent(in) :: n
      integer :: m
      m = n
      do fpart = 1, size(psize)
        if (m <= psize(fpart)) return
        m = m - psize(fpart)
      end do
    end function fpart

  end subroutine overlapping_cells

  !! This auxiliary subroutine identifies for each partition those cells that
  !! that are linked across an interface to a cell in the partition but than
  !! themselves belong to a different partition.  For a given partition the
  !! identitified cell set will be the minimal one that if added to the cells
  !! already belonging to the partition will ensure that every face linked to
  !! an interface face owned by the partition is also present on the partition.
  !! This permits the full evaluation of an interface condition at a pair of
  !! interface faces without communication.  Note that the cell containing
  !! such interface faces are necessarily also present in the partition.
  !!
  !! The LNHBR array is a link-based array: LNHBR(:,j) are the indices of the
  !! two cells linked by link j.  The cell partitioning is described by the
  !! cell partition block sizes CELL_PSIZE.  The cells identified for partition
  !! n are added to the set XCELLS(n).  This is a serial procedure.

  subroutine overlapping_cells2 (lnhbr, cell_psize, xcells)

    use integer_set_type

    integer, intent(in) :: lnhbr(:,:)     ! link cell neighbors
    integer, intent(in) :: cell_psize(:)  ! cell partition block sizes
    type(integer_set), intent(inout) :: xcells(:)

    integer :: j, p1, p2

    ASSERT(all(cell_psize >= 0))
    ASSERT(size(xcells) == size(cell_psize))
    ASSERT(size(lnhbr,1) == 2)
    ASSERT(minval(lnhbr) >= 1 .and. maxval(lnhbr) <= sum(cell_psize))

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
      do cell_part = 1, size(cell_psize)
        if (m <= cell_psize(cell_part)) return
        m = m - cell_psize(cell_part)
      end do
    end function cell_part

  end subroutine overlapping_cells2

  !! This subroutine initializes the node partition and cell node data
  !! components: NODE_IP, NNODE, NNODE_ONP, XCNODE, and CNODE.

  subroutine init_cell_node_data (this, psize, xcnode, cnode)

    use parallel_communication, only: is_IOP
    use index_partitioning, only: localize_index_struct

    type(unstr_mesh), intent(inout) :: this
    integer, intent(in) :: psize(:)
    integer, intent(in) :: xcnode(:), cnode(:)

    integer :: j
    integer, allocatable :: count_g(:), count_l(:), offP_index(:)

    call this%node_ip%init (psize)

    !! Translate the global indexing array into global row sizes.
    if (is_IOP) then
      count_g = xcnode(2:) - xcnode(:size(xcnode)-1)
    else
      allocate(count_g(0))
    end if

    call localize_index_struct (count_g, cnode, this%cell_ip, this%node_ip, count_l, this%cnode, offP_index)
    call this%node_ip%add_offP_index (offP_index)
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

  subroutine init_cell_face_data (this, psize, xcface, cface, cfpar)

    use parallel_communication, only: is_IOP, distribute
    use index_partitioning, only: localize_index_struct, gather_boundary

    type(unstr_mesh), intent(inout) :: this
    integer, intent(in) :: psize(:)
    integer, intent(in) :: xcface(:), cface(:)
    integer, intent(in) :: cfpar(:)

    integer :: j
    integer, allocatable :: count_g(:), count_l(:), offP_index(:)

    call this%face_ip%init (psize)

    !! Translate the global indexing array into global row sizes.
    if (is_IOP) then
      count_g = xcface(2:) - xcface(:size(xcface)-1)
    else
      allocate(count_g(0))
    end if

    call localize_index_struct (count_g, cface, this%cell_ip, this%face_ip, count_l, this%cface, offP_index)
    call this%face_ip%add_offP_index (offP_index)
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

    use cell_topology, only: cell_face_sizes, get_face_nodes

    type(unstr_mesh), intent(inout) :: this

    integer :: j, k, n
    integer, allocatable :: fsize(:), fnodes(:)
    integer, pointer :: list(:)

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
        list => cell_face_sizes(cell_nodes)
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
              call get_face_nodes (cell_nodes, k, fnodes, normalize=.true., &
                                   reverse=btest(this%cfpar(j),pos=k))
              face_nodes = fnodes
            end if
          end associate
        end do
      end associate
    end do
    ASSERT(minval(this%fnode) > 0)

  end subroutine init_face_node_data

  !! This subroutine initializes the face-cell connectivity component %FCELL.
  !! FCELL(:,j) are the IDs of the two cells containing the oriented face j:
  !! its normal points out of cell FCELL(1,j) into cell FCELL(2,j). If the
  !! face is on the subdomain boundary, then one of the two values will be 0.
  !! However, true domain boundary faces are always on-process and oriented
  !! such that FCELL(1,j) is non-zero and FCELL(2,j) is 0.

  subroutine init_face_cell_data(this)
    type(unstr_mesh), intent(inout) :: this
    integer :: j, k, n
    allocate(this%fcell(2,this%nface))
    this%fcell = 0
    do j = 1, this%ncell
      associate (cell_faces => this%cface(this%xcface(j):this%xcface(j+1)-1))
        do k = 1, size(cell_faces)
          n = cell_faces(k)
          if (btest(this%cfpar(j),pos=k)) then ! face oriented inward wrto cell
            this%fcell(2,n) = j
          else
            this%fcell(1,n) = j
          end if
        end do
      end associate
    end do
  end subroutine init_face_cell_data

  !! This subroutine generates a staggered array that stores
  !! the cell id of all cells using a particular node. This
  !! will become invalid near the boundary, where the neighboring cell
  !! may not exist.
  !! The cells are not ordered in any particular fashion.

  subroutine init_node_to_cell_data(this)

    type(unstr_mesh), intent(inout) :: this    


    integer :: j, n, i
    integer, allocatable ::nconnect(:)   

    allocate(nconnect(this%nnode))
    nconnect = 0
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1))
        do n = 1, size(cell_nodes)
          nconnect(cell_nodes(n)) = nconnect(cell_nodes(n)) + 1
        end do
      end associate
    end do

    allocate(this%xndcell(this%nnode+1))

    this%xndcell(1) = 1
    do n = 1, this%nnode
      this%xndcell(n+1) = this%xndcell(n) + nconnect(n)
    end do
    deallocate(nconnect)
    
    allocate(this%ndcell(this%xndcell(this%nnode+1)-1))    
    this%ndcell = -1
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1))
        do n = 1, size(cell_nodes)
          i = this%xndcell(cell_nodes(n))
          do while(this%ndcell(i) /= -1)
            i = i + 1
          end do
          this%ndcell(i) = j         
        end do
      end associate
    end do
    
  end subroutine init_node_to_cell_data

  !! This subroutine generates a staggered array that stores
  !! the cell id of all cells that share a node with a cell. This
  !! is only valid for the local mesh (1: mesh%ncell_onP)
  !! and is not ordered in any particular fashion

  subroutine init_cell_to_cell_through_node_data(this)

    use integer_vector_type
    
    type(unstr_mesh), intent(inout) :: this    


    integer :: j, n, i, k
    type(integer_vector), allocatable :: connected_cells(:)

    allocate(connected_cells(this%ncell_onP))
    do n = 1, this%nnode
      associate(cn => this%ndcell(this%xndcell(n):this%xndcell(n+1)-1))
        do j = 1, size(cn)
          if(cn(j) <= this%ncell_onP) then
             other_cell:do i = 1, size(cn)
               ! Don't connect to the cell we are adding to
               if(cn(i) == cn(j)) cycle
               
               ! Add cell if not yet added
               do k = 1, connected_cells(cn(j))%size()
                 ! Check if already added
                    if(connected_cells(cn(j))%at(k) == cn(i)) then
                    cycle other_cell
                 end if 
               end do
               ! Made it this far, new cell, add
               call connected_cells(cn(j))%push_back(cn(i))                               
             end do other_cell
          end if
        end do

      end associate
    end do
    
    allocate(this%xcnc(this%ncell_onP+1))
    this%xcnc(1) = 1
    do j = 1, this%ncell_onP
      this%xcnc(j+1) = this%xcnc(j) + connected_cells(j)%size()
    end do

    allocate(this%cnc(this%xcnc(this%ncell_onP+1)-1))
    k = 0
    do j = 1, this%ncell_onP
      do i = 1, connected_cells(j)%size()
        k = k + 1
        this%cnc(k) = connected_cells(j)%at(i)
      end do
    end do
    deallocate(connected_cells)
    
  end subroutine init_cell_to_cell_through_node_data
  
  !! This subroutine initializes the interface link data.

  subroutine init_link_data (this, mesh, lnhbr, lface)

    use parallel_communication, only: is_IOP, nPE, collate, broadcast, distribute
    use permutations, only: reorder, invert_perm
    use index_partitioning, only: localize_index_array, gather_boundary, localize_index_struct
    use bitfield_type
    use ext_exodus_mesh_type

    type(unstr_mesh), intent(inout) :: this
    type(ext_exodus_mesh), intent(inout) :: mesh
    integer, intent(inout) :: lnhbr(:,:), lface(:,:)

    integer :: j, n
    integer, allocatable :: offP_size(:), offP_index(:), cell_psize(:), psize(:), perm(:)
    integer, allocatable :: count_g(:), count_l(:), node_psize(:)
    type(bitfield), allocatable :: link_set_mask(:)

    !! Partition the links
    allocate(cell_psize(merge(nPE,0,is_IOP)), node_psize(merge(nPE,0,is_IOP)))
    call collate (cell_psize, this%cell_ip%onP_size())
    call collate (node_psize, this%node_ip%onP_size())
    if (is_IOP) then
      allocate(psize(nPE), perm(mesh%nlink), offP_size(nPE))
      call partition_links (lnhbr, cell_psize, mesh%xlnode, mesh%lnode, node_psize, psize, perm, offP_size, offP_index)
      call reorder (lface, perm)
      call reorder (lnhbr, perm)
      call reorder (mesh%xlnode, mesh%lnode, perm)
      call reorder (mesh%link_block, perm)
      call reorder (mesh%link_cell_id, perm)
      call invert_perm (perm)
      offP_index = perm(offP_index)
      deallocate(perm)
    else
      allocate(psize(0), offP_size(0), offP_index(0))
    end if

    !! THIS%LINK_IP: create the link index partition; include off-process links from above.
    call this%link_ip%init (psize, offP_size, offP_index)
    deallocate(psize, offP_size, offP_index)

    this%nlink = this%link_ip%local_size()
    this%nlink_onP = this%link_ip%onP_size()

    !! THIS%LNODE: distribute and localize the link node indexing array.
    call localize_index_array (lface, this%link_ip, this%face_ip, this%lface, offP_index)
    INSIST(size(offP_index) == 0)
    deallocate(offP_index)

    !! THIS%LNHBR: distribute and localize the link cell neighbor array.
    call localize_index_array (lnhbr, this%link_ip, this%cell_ip, this%lnhbr, offP_index)
    INSIST(size(offP_index) == 0)
    deallocate(offP_index)

    !! THIS%XLNODE, THIS%LNODE: distribute and localize the link node arrays.
    if (is_IOP) then  ! convert xlnode into equivalent count array
      count_g = mesh%xlnode(2:) - mesh%xlnode(:size(mesh%xlnode)-1)
    else
      allocate(count_g(0), mesh%lnode(0))
    end if
    call localize_index_struct (count_g, mesh%lnode, this%link_ip, this%node_ip, &
                                count_l, this%lnode, offP_index)
    INSIST(size(offP_index) == 0)
    deallocate(offP_index)
    allocate(this%xlnode(1+size(count_l)))
    this%xlnode(1) = 1
    do j = 1, size(count_l) ! convert count_l array into equivalent xlnode array
      this%xlnode(j+1) = this%xlnode(j) + count_l(j)
    end do
    deallocate(count_l)

    !! THIS%LINK_SET_ID
    n = mesh%nlblock
    call broadcast (n)
    allocate(this%link_set_id(n))
    if (is_IOP) this%link_set_id = mesh%link_block_id
    call broadcast (this%link_set_id)

    !! Convert the array of link block IDs to a bit mask.
    allocate(link_set_mask(merge(mesh%nlink,0,is_IOP)))
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

    !! THIS%LINK_SET_MASK
    allocate(this%link_set_mask(this%nlink))
    call distribute (this%link_set_mask(:this%nlink_onP), link_set_mask)
    call gather_boundary (this%link_ip, this%link_set_mask)

    !! THIS%LINK_CELL_ID
    allocate(this%link_cell_id(this%nlink))
    if (.not.is_IOP) allocate(mesh%link_cell_id(0))
    call distribute (this%link_cell_id(:this%nlink_onP), mesh%link_cell_id)
    call gather_boundary (this%link_ip, this%link_cell_id)

  end subroutine init_link_data

  !! This subroutine initializes the node set data components.  NODE_SET_MASK
  !! is a node-based integer mask array, with bit n > 0 set if the node belongs
  !! to the nth node set.  Bit 0 is set if the node belongs to a boundary face.
  !! NODE_SET_ID stores the user-assigned integer IDs for the side sets, and is
  !! replicated on each process.

  subroutine init_node_set_data (this, mesh)

    use exodus_mesh_type
    use bitfield_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary, scatter_boundary_or

    type(unstr_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: mesh

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
      INSIST(size(mesh%nset)+1 <= bit_size(node_set_mask))
      node_set_mask = 0
      do n = 1, size(mesh%nset)
        ASSERT(minval(mesh%nset(n)%node) >= 1)
        ASSERT(maxval(mesh%nset(n)%node) <= size(node_set_mask))
        do i = 1, mesh%nset(n)%num_node
          j = mesh%nset(n)%node(i)
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
    if (is_IOP) n = size(mesh%nset)
    call broadcast (n)
    allocate(this%node_set_id(n))
    if (is_IOP) this%node_set_id = mesh%nset%id
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
  !!
  !! NB: This implementation does not assume the element blocks have unique IDs.
  !! Duplicates are coalesced and the returned CELL_SET_ID array is the list of
  !! the unique IDs (important).  A consequence is that the order of CELL_SET_ID
  !! no longer corresponds to the order of the element blocks in the Exodus mesh
  !! but this was never a guaranteed property, and should not have been assumed.

  subroutine init_cell_set_data (this, mesh)

    use exodus_mesh_type
    use integer_set_type
    use permutations, only: reorder
    use parallel_communication, only: is_IOP, distribute, broadcast, collate
    use index_partitioning, only: gather_boundary

    type(unstr_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: mesh

    integer :: i, j, n, offset, ncell_tot
    integer, allocatable :: cell_set_mask(:), cell_perm(:)
    type(integer_set) :: id_set

    !! Initialize the list of cell set IDs (%CELL_SET_ID), eliminating duplicates.
    if (is_IOP) then
      do n = 1, size(mesh%eblk)
        call id_set%add (mesh%eblk(n)%id)
      end do
      n = id_set%size()
    end if
    call broadcast (n)
    allocate(this%cell_set_id(n))
    if (is_IOP) this%cell_set_id = id_set
    call broadcast (this%cell_set_id)

    ncell_tot = this%cell_ip%global_size()

    !! Generate the global cell_set mask array (original cell ordering)
    allocate(cell_set_mask(merge(ncell_tot,0,is_IOP)))
    if (is_IOP) then
      INSIST(size(this%cell_set_id)+1 <= bit_size(cell_set_mask))
      offset = 0
      do n = 1, size(mesh%eblk)
        do j = size(this%cell_set_id), 1, -1
          if (this%cell_set_id(j) == mesh%eblk(n)%id) exit
        end do
        INSIST(j > 0)
        do i = 1, mesh%eblk(n)%num_elem
          cell_set_mask(offset+i) = ibset(0, pos=j)
        end do
        offset = offset + mesh%eblk(n)%num_elem
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

  end subroutine init_cell_set_data

  !! This subroutine initializes the face set data components.  FACE_SET_MASK
  !! is a face-based integer mask array, with bit n > 0 set if the face belongs
  !! to the nth face set.  CELL_SET_ID stores the user-assigned integer IDs for
  !! the cell sets, and is replicated on each process.  Currently, each Exodus
  !! side set is mapped to a face set whose ID is the side set ID.  An exodus
  !! "side" is a pair of indices identifying a cell and one of its sides.  This
  !! naturally identifies a unique mesh face, but any orientation information
  !! implicit with the side-of-a-cell description is lost.

  subroutine init_face_set_data (this, mesh, xcface, cface)

    use bitfield_type
    use exodus_mesh_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary

    type(unstr_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: mesh
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
      INSIST(size(mesh%sset)+1 <= bit_size(face_set_mask))
      face_set_mask = ZERO_BITFIELD
      do n = 1, size(mesh%sset)
        ASSERT(minval(mesh%sset(n)%elem) >= 1)
        ASSERT(maxval(mesh%sset(n)%elem) <= size(xcface)-1)
        ASSERT(minval(mesh%sset(n)%face) >= 1)
        do i = 1, mesh%sset(n)%num_side
          j = xcface(mesh%sset(n)%elem(i)) + mesh%sset(n)%face(i) - 1 ! index of side in CFACE
          ASSERT(j < xcface(mesh%sset(n)%elem(i)+1))
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
    if (is_IOP) n = size(mesh%sset)
    call broadcast (n)
    allocate(this%face_set_id(n))
    if (is_IOP) this%face_set_ID = mesh%sset%id
    call broadcast (this%face_set_id)

  end subroutine init_face_set_data

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

end module unstr_mesh_factory
