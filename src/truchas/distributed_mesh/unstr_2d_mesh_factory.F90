#include "f90_assert.fpp"

module unstr_2d_mesh_factory

  use unstr_2d_mesh_type
  use truchas_logging_services
  implicit none
  private

  public :: new_unstr_2d_mesh, unstr_2d_mesh

  interface new_unstr_2d_mesh
    procedure new_unstr_2d_mesh_regular
  end interface

contains

  function new_unstr_2d_mesh_regular(xmin, xmax, nx, eps) result(this)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use ext_exodus_mesh_type
    use parallel_communication, only: is_IOP
    use parameter_list_type

    real(r8), intent(in) :: xmin(:), xmax(:)
    integer,  intent(in) :: nx(:)
    real(r8), intent(in), optional :: eps
    type(unstr_2d_mesh), pointer :: this

    integer :: i, j, n
    real(r8) :: x(2)
    type(ext_exodus_mesh) :: mesh ! temporary serial base mesh
    type(parameter_list) :: params
    character(:), allocatable :: errmsg

    ASSERT(size(xmin) == 2)
    ASSERT(size(xmax) == 2)
    ASSERT(size(nx) == 2)
    ASSERT(all(xmin < xmax))
    ASSERT(all(nx > 0))

    if (is_IOP) then

      mesh%num_dim = 2
      mesh%num_elem = product(nx)
      mesh%num_node = product(nx+1)

      !! Generate the node coordinates.
      allocate(mesh%coord(2,mesh%num_node))
      do j = 1, nx(2)+1
        x(2) = ((nx(2)-j-1)/real(nx(2),r8))*xmin(2) + ((j-1)/real(nx(2),r8))*xmax(2)
        do i = 1, nx(1)+1
          x(1) = ((nx(1)-i-1)/real(nx(1),r8))*xmin(1) + ((i-1)/real(nx(1),r8))*xmax(1)
          mesh%coord(:,node_index(i,j)) = x
        end do
      end do
      if (present(eps)) call randomize_coord(eps)

      !! All cells go into a single element block (ID=1)
      mesh%num_eblk = 1
      allocate(mesh%eblk(mesh%num_eblk))
      mesh%eblk(1)%id = 1
      mesh%eblk(1)%num_elem = mesh%num_elem
      mesh%eblk(1)%num_nodes_per_elem = 4
      mesh%eblk(1)%elem_type = 'QUAD'
      allocate(mesh%eblk(1)%connect(4,mesh%eblk(1)%num_elem))
      !! Generate the regular hex connectivity
      do concurrent (i = 1:nx(1), j = 1:nx(2))
        n = cell_index(i,j)
        mesh%eblk(1)%connect(1,n) = node_index(i,j)
        mesh%eblk(1)%connect(2,n) = node_index(i+1,j)
        mesh%eblk(1)%connect(3,n) = node_index(i+1,j+1)
        mesh%eblk(1)%connect(4,n) = node_index(i,j+1)
      end do

      !! Generate side sets (left/right/bottom/top)
      mesh%num_sset = 4
      allocate(mesh%sset(mesh%num_sset))

      !! Side set 1 (x = xmin)
      mesh%sset(1)%id = 1
      mesh%sset(1)%num_side = nx(2)
      allocate(mesh%sset(1)%elem(nx(2)), mesh%sset(1)%face(nx(2)))
      do concurrent (j = 1:nx(2))
        mesh%sset(1)%elem(j) = cell_index(1,j)
      end do
      mesh%sset(1)%face = 4

      !! Side set 2 (x = xmax)
      mesh%sset(2)%id = 2
      mesh%sset(2)%num_side = nx(2)
      allocate(mesh%sset(2)%elem(nx(2)), mesh%sset(2)%face(nx(2)))
      do concurrent (j = 1:nx(2))
        mesh%sset(2)%elem(j) = cell_index(nx(1),j)
      end do
      mesh%sset(2)%face = 2

      !! Side set 3 (y = ymin)
      mesh%sset(3)%id = 3
      mesh%sset(3)%num_side = nx(1)
      allocate(mesh%sset(3)%elem(nx(1)), mesh%sset(3)%face(nx(1)))
      do concurrent (i = 1:nx(1))
        mesh%sset(3)%elem(i) = cell_index(i,1)
      end do
      mesh%sset(3)%face = 1

      !! Side set 4 (y = ymax)
      mesh%sset(4)%id = 4
      mesh%sset(4)%num_side = nx(1)
      allocate(mesh%sset(4)%elem(nx(1)), mesh%sset(4)%face(nx(1)))
      do concurrent (i = 1:nx(1))
        mesh%sset(4)%elem(i) = cell_index(i,nx(2))
      end do
      mesh%sset(4)%face = 3

      !! No node sets, no internal interfaces
      mesh%num_nset = 0
      allocate(mesh%nset(0))
      call mesh%set_no_links

    end if

    this => new_unstr_2d_mesh_aux(mesh, params, errmsg)
    INSIST(associated(this))

  contains

    pure integer function cell_index(i, j)
      integer, intent(in) :: i, j
      cell_index = i + (j-1)*nx(1)
    end function

    pure integer function node_index(i, j)
      integer, intent(in) :: i, j
      node_index = i + (j-1)*(nx(1)+1)
    end function

    ! Apply a random perturbation to the coordinates (preserves the domain)
    subroutine randomize_coord(eps)
      real(r8), intent(in) :: eps
      integer :: i, j, n
      logical :: mask(2)
      real(r8) :: dx(2)
      ASSERT(eps >= 0)
      if (eps == 0) return
      do j = 1, nx(2)
        mask(2) = (j > 1 .and. j < nx(2))
        do i = 1, nx(1)
          mask(1) = (i > 1 .and. i < nx(1))
          call random_number(dx)  ! in [0,1)
          dx = eps*(2*dx - 1)     ! in [-eps, eps)
          n = node_index(i,j)
          mesh%coord(:,n) = mesh%coord(:,n) + merge(dx, 0.0_r8, mask)
        end do
      end do
    end subroutine randomize_coord

  end function new_unstr_2d_mesh_regular

  function new_unstr_2d_mesh_aux(mesh, params, errmsg) result(this)

    use unstr_2d_mesh_type
    use ext_exodus_mesh_type
    use permutations
    use simple_partitioning_methods, only: get_block_partition, read_partition
    use index_partitioning
    use unstr_2d_mesh_tools
    use parallel_communication
    use parameter_list_type

    type(ext_exodus_mesh), intent(inout) :: mesh
    type(parameter_list),  intent(inout) :: params
    character(:), allocatable, intent(out) :: errmsg
    type(unstr_2d_mesh), pointer :: this

    integer :: j, k, n, nnode, nface, ncell, stat, pfirst
    integer :: cell_psize(nPE), node_psize(nPE), face_psize(nPE)
    integer, allocatable :: xcnode(:), cnode(:), xcnhbr(:), cnhbr(:), part(:)
    integer, allocatable :: xcface(:), cface(:), cfpar(:)
    integer, allocatable :: cell_perm(:), node_perm(:), offP_size(:), offP_index(:)
    integer, allocatable :: perm(:)
    character(:), allocatable :: string

    this => null()
    ncell = mesh%num_elem
    nnode = mesh%num_node

    !! Flatten the Exodus element block structure.
    call mesh%get_concat_elem_conn(xcnode, cnode)

    !! Generate the cell neighbor array.
    if (is_IOP) then
      call TLS_info('  finding cell neighbors', TLS_VERB_NORMAL)
      call get_cell_neighbor_array(xcnode, cnode, xcnhbr, cnhbr, stat)
      if (stat /= 0) errmsg = 'get_cell_neighbor_array: invalid mesh topology detected'
    else
      allocate(xcnhbr(1), cnhbr(0))
      xcnhbr(1) = 1
    end if
    call broadcast_status(stat, errmsg)
    if (stat /= 0) return

    !! Partition and order the cells.
    allocate(cell_perm(ncell))
    if (is_IOP) then
      call TLS_info('  partitioning the mesh cells', TLS_VERB_NORMAL)
      !! Partition the cell neighbor graph.
      allocate(part(mesh%num_elem))
      call params%get('partitioner', string, default='chaco')
      if (nPE == 1) then
        part = 1
        stat = 0
      else if (string == 'block') then
        call get_block_partition(nPE, part)
        stat = 0
      else if (string == 'file') then
        call params%get('partition-file', string)
        call params%get('first-partition', pfirst, default=0)
        call read_partition(string, pfirst, nPE, part, stat, errmsg)
        if (stat /= 0) errmsg = 'error reading cell partition: ' // errmsg
      else
        call partition_cells(params, xcnhbr, cnhbr, nPE, part, stat, errmsg)
        if (stat /= 0) errmsg = 'error computing cell partition: ' // errmsg
      end if
    end if
    call broadcast_status(stat, errmsg)
    if (stat /= 0) return

    if (is_IOP) then
      !! Compute the partition sizes and the permutation making this a block partition.
      call blocked_partition(part, cell_psize, cell_perm)
      !! Reorder cell-based arrays.
      call reorder(xcnode, cnode, cell_perm)
      call reorder(xcnhbr, cnhbr, cell_perm)
      !! Map the values of cell-valued arrays.
      allocate(perm(size(cell_perm)))
      call invert_perm(cell_perm, perm)
      do j = 1, size(cnhbr)
        if (cnhbr(j) > 0) cnhbr(j) = perm(cnhbr(j))
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
      call partition_facets(xcnode, cnode, cell_psize, node_psize, node_perm)
      !! Reorder node-based arrays.
      call reorder(mesh%coord, node_perm)
      !! Map the values of node-valued arrays.
      allocate(perm(size(node_perm)))
      call invert_perm(node_perm, perm)
      do j = 1, size(cnode)
        cnode(j) = perm(cnode(j))
      end do
      do n = 1, size(mesh%nset)
        mesh%nset(n)%node = perm(mesh%nset(n)%node)
      end do
      deallocate(perm)
    end if

    !! Enumerate and partition the mesh faces.
    allocate(cfpar(ncell))
    if (is_IOP) then
      call TLS_info('  numbering the mesh faces', TLS_VERB_NORMAL)
      call label_mesh_faces(xcnode, cnode, nface, xcface, cface)
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
      call TLS_info('  partitioning the mesh faces', TLS_VERB_NORMAL)
      allocate(perm(nface))
      call partition_facets(xcface, cface, cell_psize, face_psize, perm)
      call invert_perm (perm)
      do j = 1, size(cface)
        cface(j) = perm(cface(j))
      end do
      deallocate(perm)
    else
      allocate(xcface(1), cface(0))
      xcface(1) = 1
    end if

    !! Identify off-process ghost cells to include with each partition.
    call TLS_info('  identifying off-process ghost cells', TLS_VERB_NORMAL)
    call select_ghost_cells(cell_psize, xcnhbr, cnhbr, xcnode, cnode, node_psize, &
                             xcface, cface, face_psize, offP_size, offP_index)
    deallocate(xcnhbr, cnhbr)

    !! Begin initializing the unstr_2d_mesh result object.
    call TLS_info('  generating parallel mesh structure')
    allocate(this)

    !! Create the cell index partition; include the off-process cells from above.
    call this%cell_ip%init(cell_psize, offP_size, offP_index)
    deallocate(offP_size, offP_index)

    this%ncell = this%cell_ip%local_size()
    this%ncell_onP = this%cell_ip%onP_size()

    !! Distribute the cell permutation array; gives mapping to the external cell number.
    allocate(this%xcell(this%ncell))
    call distribute(this%xcell(:this%ncell_onP), cell_perm)
    call gather_boundary(this%cell_ip, this%xcell)
    deallocate(cell_perm)

    !! Create the node index partition and localize the global CNODE array,
    !! which identifies off-process nodes to augment the partition with.
    call init_cell_node_data(this, node_psize, xcnode, cnode)

    !! Distribute the node permutation array; gives mapping to the external node number.
    allocate(this%xnode(this%nnode))
    call distribute(this%xnode(:this%nnode_onP), node_perm)
    call gather_boundary(this%node_ip, this%xnode)
    deallocate(node_perm)

    !! Create the face index partition and localize the global CFACE array,
    !! which identifies off-process faces to augment the partition with.
    call init_cell_face_data(this, face_psize, xcface, cface, cfpar)
    deallocate(cfpar)

    !! Initialize the secondary face-node indexing array.
    call init_face_node_data(this)
    call init_face_cell_data(this)

    !! Generate the cell neighbor data for each subdomain.
    call get_cell_neighbor_array(this%xcnode, this%cnode, this%xcnhbr, this%cnhbr, stat)
    INSIST(stat == 0)

    !! Initialize the node, face, and cell set data.
    call init_face_set_data (this, mesh, xcface, cface)
    call init_node_set_data (this, mesh)
    call init_cell_set_data (this, mesh)
    deallocate(xcface, cface)

    !! Scale the node coordinates and distribute.
    if (.not.is_IOP) then
      if (.not.allocated(mesh%coord)) allocate(mesh%coord(2,0))
    end if
    allocate(this%x(2,this%nnode))
    call distribute(this%x(:,:this%nnode_onP), mesh%coord)
    call gather_boundary(this%node_ip, this%x)

    !! Initialize the mesh geometry data components.
    allocate(this%volume(this%ncell), this%normal(2,this%nface), this%area(this%nface))
    call this%compute_geometry

  end function new_unstr_2d_mesh_aux

!!!! AUXILIARY PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! When a procedure executed only on the IO processor returns a status
  !! and possible error message, these need to be communicated to the other
  !! processes.  This auxiliary procedure performs that communication.
  !! The IO process value of STAT is broadcast, and if not 0, the IO process
  !! value of ERRMSG, which must be allocated, is broadcast also.

  subroutine broadcast_status(stat, errmsg)
    use parallel_communication, only: is_IOP, broadcast
    integer, intent(inout) :: stat
    character(:), allocatable, intent(inout) :: errmsg
    integer :: n
    call broadcast(stat)
    if (stat /= 0) then
      if (is_IOP) n = len(errmsg)
      call broadcast(n)
      if (.not.is_IOP) allocate(character(len=n)::errmsg)
      call broadcast(errmsg)
    end if
  end subroutine broadcast_status

  subroutine partition_cells(params, xcnhbr, cnhbr, npart, part, stat, errmsg)

    use graph_type
    use graph_partitioner_factory
    use parameter_list_type

    type(parameter_list) :: params
    integer, intent(in)  :: xcnhbr(:), cnhbr(:) ! cell neighbor array
    integer, intent(in)  :: npart   ! number of partitions
    integer, intent(out) :: part(:) ! cell partition assignment
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, k, ncell

    type(graph), allocatable :: g
    integer, allocatable :: xadj(:), adjncy(:)
    real, allocatable :: ewgt(:)
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
          if (list(k) > 0) call g%add_edge(j, list(k))
        end do
      end associate
    end do

    call g%get_adjacency(xadj, adjncy)
    deallocate(g)

    !! Define edge weights.
    allocate(ewgt(size(adjncy)))
    ewgt = 1.0  ! the default

    call alloc_graph_partitioner(gpart, params)
    call gpart%compute(ncell, xadj, adjncy, ewgt, npart, part, stat, errmsg)

  end subroutine partition_cells

  !! Given a partition assignment array PART, this auxiliary subroutine
  !! computes the permutation array PERM that makes the partition a block
  !! partition.  The partition sizes are returned in PSIZE.  PERM maps the
  !! new numbering to the original, and preserves the relative order of
  !! elements in the same partition.  Specifically, 1) PART(PERM(:)) is
  !! sorted (non-decreasing), and 2) j < k whenever PERM(j) < PERM(k) and
  !! PART(PERM(j)) = PART(PERM(k)).

  subroutine blocked_partition(part, psize, perm)

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

  !! This auxiliary subroutine partitions the mesh facets of one type (nodes
  !! or faces) and generates a new numbering of them for which the
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

  subroutine partition_facets(xfacet, facet, cell_psize, psize, perm)

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

    call invert_perm(perm1)  ! PERM1 is now good-to-old numbering

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
    call reorder(part, perm1)

    !! Block partition permutation (new-to-good).
    call blocked_partition(part, psize, perm)

    !! Total facet permutation (new-to-old).
    do j = 1, size(perm)
      perm(j) = perm1(perm(j))
    end do

  end subroutine partition_facets

  !! This auxiliary subroutine identifies the off-process ghost cells that
  !! should be added to each subdomain.  The criterion is that if a cell
  !! contains a face or node that belongs to a subdomain and the cell itself
  !! does not belong to the subdomain, then it is added to the subdomain as
  !! a ghost cell.  This ensures that every face and node belonging to the
  !! subdomain will have complete cell support.  The motivation for this
  !! stems from the finite element context in which the equation for a DoF
  !! located at a facet depends on calculations over all the cells that
  !! contain the facet. This is a serial procedure.

  subroutine select_ghost_cells (cell_psize, xcnhbr, cnhbr, xcnode, cnode, node_psize, &
      xcface, cface, face_psize, offP_size, offP_index)

    use integer_set_type
    use parallel_communication, only: is_IOP, nPE

    integer, intent(in) :: cell_psize(:), node_psize(:), face_psize(:)
    integer, intent(in) :: xcnhbr(:), cnhbr(:)
    integer, intent(in) :: xcnode(:), cnode(:)
    integer, intent(in) :: xcface(:), cface(:)
    integer, allocatable, intent(out) :: offP_size(:), offP_index(:)

    integer :: n, offset
    type(integer_set), allocatable :: ghosts(:)

    if (is_IOP) then
      allocate(ghosts(nPE))
      call all_cell_neighbors(xcnode, cnode, cell_psize, ghosts)
      !! Copy the sets into packed array storage
      offP_size = ghosts%size()
      n = sum(offP_size)
      allocate(offP_index(n))
      offset = 0
      do n = 1, nPE
        call ghosts(n)%copy_to_array(offP_index(offset+1:))
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

  subroutine all_cell_neighbors(xcnode, cnode, cell_psize, xcells)

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

  !! This subroutine initializes the node partition and cell node data
  !! components: NODE_IP, NNODE, NNODE_ONP, XCNODE, and CNODE.

  subroutine init_cell_node_data (this, psize, xcnode, cnode)

    use parallel_communication, only: is_IOP
    use index_partitioning, only: localize_index_struct

    type(unstr_2d_mesh), intent(inout) :: this
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

    call localize_index_struct(count_g, cnode, this%cell_ip, this%node_ip, count_l, this%cnode, offP_index)
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

    type(unstr_2d_mesh), intent(inout) :: this
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

    call localize_index_struct(count_g, cface, this%cell_ip, this%face_ip, count_l, this%cface, offP_index)
    call this%face_ip%add_offP_index(offP_index)
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
    call distribute(this%cfpar(:this%ncell_onP), cfpar)
    call gather_boundary(this%cell_ip, this%cfpar)

    this%nface = this%face_ip%local_size()
    this%nface_onP = this%face_ip%onP_size()

  end subroutine init_cell_face_data

  !! This subroutine initializes the face-node connectivity data components,
  !! stored in the array FNODE. The face-node lists are oriented consistently
  !! with the CFPAR data.

  subroutine init_face_node_data(this)

    use cell_topology_2d, only: get_face_nodes

    type(unstr_2d_mesh), intent(inout) :: this

    integer :: j, k, n
    integer, allocatable :: fnodes(:)

    ASSERT(size(this%xcnode) == size(this%xcface))
    ASSERT(size(this%cfpar) == size(this%xcface)-1)
    ASSERT(size(this%cnode) == this%xcnode(size(this%xcnode))-1)
    ASSERT(size(this%cface) == this%xcface(size(this%xcface))-1)
    ASSERT(minval(this%cnode) > 0)
    ASSERT(minval(this%cface) > 0)

    !! Fill the FNODE array.
    allocate(this%fnode(2,this%nface))
    this%fnode = 0
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1), &
                 cell_faces => this%cface(this%xcface(j):this%xcface(j+1)-1))
        do k = 1, size(cell_faces)
          n = cell_faces(k)
          if (this%fnode(1,n) == 0) then
            call get_face_nodes(cell_nodes, k, fnodes, reverse=btest(this%cfpar(j),pos=k))
            this%fnode(:,n) = fnodes
          end if
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
    type(unstr_2d_mesh), intent(inout) :: this
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

  !! This subroutine initializes the face set data components.  FACE_SET_MASK
  !! is a face-based integer mask array, with bit n > 0 set if the face belongs
  !! to the nth face set.  CELL_SET_ID stores the user-assigned integer IDs for
  !! the cell sets, and is replicated on each process.  Currently, each Exodus
  !! side set is mapped to a face set whose ID is the side set ID.  An exodus
  !! "side" is a pair of indices identifying a cell and one of its sides.  This
  !! naturally identifies a unique mesh face, but any orientation information
  !! implicit with the side-of-a-cell description is lost.

  subroutine init_face_set_data(this, mesh, xcface, cface)

    use bitfield_type
    use exodus_mesh_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary

    type(unstr_2d_mesh), intent(inout) :: this
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
    call distribute(this%face_set_mask(:this%nface_onP), face_set_mask)
    call gather_boundary(this%face_ip, this%face_set_mask)
    deallocate(face_set_mask)

    !! Initialize the list of cell set IDs (%FACE_SET_ID)
    if (is_IOP) n = size(mesh%sset)
    call broadcast (n)
    allocate(this%face_set_id(n))
    if (is_IOP) this%face_set_ID = mesh%sset%id
    call broadcast(this%face_set_id)

  end subroutine init_face_set_data

  !! This subroutine initializes the node set data components.  NODE_SET_MASK
  !! is a node-based integer mask array, with bit n > 0 set if the node belongs
  !! to the nth node set.  Bit 0 is set if the node belongs to a boundary face.
  !! NODE_SET_ID stores the user-assigned integer IDs for the side sets, and is
  !! replicated on each process.

  subroutine init_node_set_data(this, mesh)

    use exodus_mesh_type
    use bitfield_type
    use parallel_communication, only: is_IOP, distribute, broadcast
    use index_partitioning, only: gather_boundary, scatter_boundary_or

    type(unstr_2d_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: mesh

    integer :: i, j, n, nnode_tot
    integer, allocatable :: node_set_mask(:)
    logical, allocatable :: bnode(:)

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
    call distribute(this%node_set_mask(:this%nnode_onP), node_set_mask)
    call gather_boundary(this%node_ip, this%node_set_mask)
    deallocate(node_set_mask)

    !! Initialize the list of node set IDs (%NODE_SET_ID)
    if (is_IOP) n = size(mesh%nset)
    call broadcast(n)
    allocate(this%node_set_id(n))
    if (is_IOP) this%node_set_id = mesh%nset%id
    call broadcast(this%node_set_id)

    !! Tag boundary nodes in the node set mask (bit 0).
    allocate(bnode(this%nnode))
    bnode = .false.
    do j = 1, this%nface
      if (btest(this%face_set_mask(j),pos=0)) then  ! boundary face
        bnode(this%fnode(:,j)) = .true.
      end if
    end do
    call scatter_boundary_or(this%node_ip, bnode)
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

  subroutine init_cell_set_data(this, mesh)

    use exodus_mesh_type
    use integer_set_type
    use permutations, only: reorder
    use parallel_communication, only: is_IOP, distribute, broadcast, collate
    use index_partitioning, only: gather_boundary

    type(unstr_2d_mesh), intent(inout) :: this
    class(exodus_mesh), intent(in) :: mesh

    integer :: i, j, n, offset, ncell_tot
    integer, allocatable :: cell_set_mask(:), cell_perm(:)
    type(integer_set) :: id_set

    !! Initialize the list of cell set IDs (%CELL_SET_ID), eliminating duplicates.
    if (is_IOP) then
      do n = 1, size(mesh%eblk)
        call id_set%add(mesh%eblk(n)%id)
      end do
      n = id_set%size()
    end if
    call broadcast(n)
    allocate(this%cell_set_id(n))
    if (is_IOP) this%cell_set_id = id_set
    call broadcast(this%cell_set_id)

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
    call collate(cell_perm, this%xcell(:this%ncell_onP))
    if (is_IOP) call reorder(cell_set_mask, cell_perm)
    deallocate(cell_perm)

    !! Initialize the distributed cell set mask (%CELL_SET_MASK)
    allocate(this%cell_set_mask(this%ncell))
    call distribute(this%cell_set_mask(:this%ncell_onP), cell_set_mask)
    call gather_boundary(this%cell_ip, this%cell_set_mask)
    deallocate(cell_set_mask)

  end subroutine init_cell_set_data

end module unstr_2d_mesh_factory
