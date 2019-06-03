!!
!! UNSTR_MESH_TYPE
!!
!! This module provides a derived type that encapsulates the data describing a
!! distributed unstructured mixed-element mesh.  Supported element types are
!! hexes, tets, pyramids, and wedges/prisms.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Revised May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type UNSTR_MESH that encapsulates the data
!!  describing a distributed unstructured mixed-element mesh.  On each process
!!  the object describes a complete mesh of some subdomain that references only
!!  local entities (cells, faces, nodes), and can rightly be considered as a
!!  serial mesh for that subdomain.  Globally, the mesh-conforming subdomains
!!  will overlap, perhaps only along a boundary but more generally on some
!!  collection of cells.  Additional mesh data describes this overlap and
!!  provides for communication between overlapping entities.
!!
!!  Objects of this type are intended to be used by trusted code, and so its
!!  data components are public.  However, the components must be treated as
!!  read-only because an object may be shared amongst multiple clients.  The
!!  following data components are accessible:
!!
!!    nnode, nface, ncell - the number of nodes, faces, and cells in the mesh.
!!
!!    xcnode, cnode - pair of rank-1 integer arrays storing the cell-node data:
!!        cnode(xcnode(j):xcnode(j+1)-1) is the ordered list of node indices
!!        defining cell j.  The shape of xcnode is [ncell+1] and the shape of
!!        cnode is [xcnode(ncell+1)-1].
!!
!!    xcface, cface - pair of rank-1 integer arrays storing the cell-face data:
!!        cface(xcface(j):xcface(j+1)-1) is the ordered list of face indices
!!        belonging to cell j.  The shape of xcface is [ncell+1] and the shape
!!        of cface is [xcface(ncell+1)-1].
!!
!!    xfnode, fnode - pair of rank-1 integer arrays storing the face-node data:
!!        fnode(xfnode(j):xfnode(j+1)-1) is the ordered list of node indices
!!        defining the oriented face j.  The shape of xfnode is [nface+1] and
!!        the shape of fnode is [xfnode(nface+1)-1].
!!
!!    cfpar - an integer bit mask array storing the relative cell face
!!        orientations: btest(cfpar(j),k) is true when face k of cell j is
!!        inward oriented with respect to cell j, and false when it is
!!        outward oriented.  The shape of cfpar is [ncell].
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
!!    nnode_onP, nface_onP, ncell_onP - the number of local nodes, faces, and
!!        cells that that are uniquely owned (on-process).
!!
!!    node_ip, face_ip, cell_ip - derived types that describe the partitioning
!!        and overlap of nodes, edges, faces, and cells, including information
!!        necessary to communicate off-process data between processes.
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
!!  GEOMETRY DATA:
!!
!!    x - the rank-2 real array of node coordinates; x(:,j) is the position in
!!        R^3 of node j.  Its shape is [3,nnode].
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
!!  INTERFACE LINKS:
!!
!!    nlink, nlink_onP - the number of interface links and uniquely owned
!!        (on-process) interface links in the (subdomain) mesh.
!!
!!    lface - the rank-1 integer link-face array: lface(:,j) are the indices
!!        of the two opposing mesh faces for link j.
!!
!!    link_set_id - a rank-1 integer array storing the unique link set IDs.
!!        This data is replicated on each process.
!!
!!    link_set_mask - a rank-1 bitmask array: btest(link_set_mask(j),k)
!!        returns true if link j belongs to the link set with ID link_set_id(k).
!!
!!    link_ip - derived type that describes the partitioning and overlap of
!!        links, including information necessary to comminicate off-process
!!        data between processes.
!!

#include "f90_assert.fpp"

module unstr_mesh_type

  use kinds, only: r8
  use base_mesh_class
  use index_partitioning
  use parallel_communication
  use bitfield_type
  use cell_topology
  use f08_intrinsics
  implicit none
  private

  type, extends(base_mesh), public :: unstr_mesh
    integer, allocatable :: xcnode(:), cnode(:) ! cell nodes
    integer, allocatable :: xcface(:), cface(:) ! cell faces
    integer, allocatable :: xfnode(:), fnode(:) ! face nodes
    integer, allocatable :: xcnhbr(:), cnhbr(:) ! cell neighbors
    integer, allocatable :: cfpar(:)  ! relative cell face orientation (bit mask)
    integer, allocatable :: fcell(:,:)  ! face cell neighbors
    real(r8), allocatable :: normal(:,:)
    real(r8), allocatable :: cell_centroid(:,:)
    real(r8), allocatable :: face_centroid(:,:)
    real(r8), allocatable :: face_normal_dist(:)  ! minimum distance from face along edges
    !! Mesh interface links.
    integer :: nlink = 0, nlink_onP = 0
    integer, allocatable :: lface(:,:)        ! pointer due to localize_index_array
    integer, allocatable :: link_set_id(:)    ! user-assigned ID for each link block
    type(bitfield), allocatable :: link_set_mask(:)  ! link block index
    type(ip_desc) :: link_ip
    !! Additional link data aiding transition from old mesh.
    integer, allocatable :: link_cell_id(:)   ! external cell ID the link was derived from (or 0)
    integer, allocatable :: lnhbr(:,:)        ! link cell neighbors (2)
    integer, allocatable :: xlnode(:), lnode(:) ! link nodes
    integer, allocatable :: parent_node(:)    ! node parents (global ID)
  contains
    procedure :: get_global_cnode_array
    procedure :: get_global_cface_array
    procedure :: compute_geometry
    procedure :: write_profile
    procedure :: check_bndry_face_set
    procedure :: get_link_set_ids
    procedure :: init_cell_centroid
    procedure :: init_face_centroid
    procedure :: init_face_normal_dist
    procedure :: nearest_node
    procedure :: nearest_cell
  end type unstr_mesh

contains

  !! Compute the geometric data components from the node coordinates.
  subroutine compute_geometry (this)
    use cell_geometry, only: cell_volume, face_normal, vector_length
    class(unstr_mesh), intent(inout) :: this
    integer :: j
    ASSERT(allocated(this%volume))
    ASSERT(allocated(this%normal))
    ASSERT(allocated(this%area))
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1))
        this%volume(j) = cell_volume(this%x(:,cell_nodes))
      end associate
    end do
    do j = 1, this%nface
      associate(face_nodes => this%fnode(this%xfnode(j):this%xfnode(j+1)-1))
        this%normal(:,j) = face_normal(this%x(:,face_nodes))
        this%area(j) = vector_length(this%normal(:,j))
      end associate
    end do
  end subroutine compute_geometry

  subroutine init_cell_centroid(this)
    class(unstr_mesh), intent(inout) :: this
    integer :: j
    if (allocated(this%cell_centroid)) return
    allocate(this%cell_centroid(3,this%ncell))
    do j = 1, this%ncell
      associate(cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1))
        this%cell_centroid(:,j) = sum(this%x(:,cell_nodes),dim=2)/size(cell_nodes)
      end associate
    end do
  end subroutine init_cell_centroid

  subroutine init_face_centroid(this)
    class(unstr_mesh), intent(inout) :: this
    integer :: j
    if (allocated(this%face_centroid)) return
    allocate(this%face_centroid(3,this%nface))
    do j = 1, this%nface
      associate(face_nodes => this%fnode(this%xfnode(j):this%xfnode(j+1)-1))
        this%face_centroid(:,j) = sum(this%x(:,face_nodes),dim=2)/size(face_nodes)
      end associate
    end do
  end subroutine init_face_centroid

  subroutine init_face_normal_dist(this)
    class(unstr_mesh), intent(inout) :: this
    integer :: i,j,k

    if (allocated(this%face_normal_dist)) return
    allocate(this%face_normal_dist(this%nface))

    this%face_normal_dist(:) = huge(1.0_r8)
    do i = 1, this%ncell
      associate(faces => this%cface(this%xcface(i):this%xcface(i+1)-1), &
          cell_nodes => this%cnode(this%xcnode(i):this%xcnode(i+1)-1))

        do k = 1, size(faces)

          select case(size(cell_nodes))
          case (4)
            associate (face_list => TET4_FACES(TET4_XFACE(k):TET4_XFACE(k+1)-1))
              do j = 1, size(cell_nodes)
                if (findloc(face_list, j) == 0) then
                  this%face_normal_dist(faces(k)) = min(this%face_normal_dist(faces(k)), &
                      abs(dot_product(this%x(:,cell_nodes(j))-this%face_centroid(:,faces(k)), &
                      this%normal(:,faces(k))/this%area(faces(k)))))
                end if
              end do
            end associate
          case (5)
            associate (face_list => PYR5_FACES(PYR5_XFACE(k):PYR5_XFACE(k+1)-1))
              do j = 1, size(cell_nodes)
                if (findloc(face_list, j) == 0) then
                  this%face_normal_dist(faces(k)) = min(this%face_normal_dist(faces(k)), &
                      abs(dot_product(this%x(:,cell_nodes(j))-this%face_centroid(:,faces(k)), &
                      this%normal(:,faces(k))/this%area(faces(k)))))
                end if
              end do
            end associate
          case (6)
            associate (face_list => WED6_FACES(WED6_XFACE(k):WED6_XFACE(k+1)-1))
              do j = 1, size(cell_nodes)
                if (findloc(face_list, j) == 0) then
                  this%face_normal_dist(faces(k)) = min(this%face_normal_dist(faces(k)), &
                      abs(dot_product(this%x(:,cell_nodes(j))-this%face_centroid(:,faces(k)), &
                      this%normal(:,faces(k))/this%area(faces(k)))))
                end if
              end do
            end associate
          case (8)
            associate (face_list => HEX8_FACES(HEX8_XFACE(k):HEX8_XFACE(k+1)-1))
              do j = 1, size(cell_nodes)
                if (findloc(face_list, j) == 0) then
                  this%face_normal_dist(faces(k)) = min(this%face_normal_dist(faces(k)), &
                      abs(dot_product(this%x(:,cell_nodes(j))-this%face_centroid(:,faces(k)), &
                      this%normal(:,faces(k))/this%area(faces(k)))))
                end if
              end do
            end associate
          end select
        end do
      end associate
    end do
  end subroutine init_face_normal_dist

  !! Creates the global ragged CNODE array on the IO process, 0-sized on others.
  subroutine get_global_cnode_array (this, xcnode, cnode)
    class(unstr_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: xcnode(:), cnode(:)
    associate (xcnode_onP => this%xcnode(:this%ncell_onP+1), &
                cnode_onP => this%cnode(:this%xcnode(this%ncell_onP+1)-1))
      call get_global_ragged_array (xcnode_onP, this%node_ip%global_index(cnode_onP), xcnode, cnode)
    end associate
  end subroutine get_global_cnode_array

  !! Creates the global ragged CFACE array on the IO process, 0-sized on others.
  subroutine get_global_cface_array (this, xcface, cface)
    class(unstr_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: xcface(:), cface(:)
    associate (xcface_onP => this%xcface(:this%ncell_onP+1), &
                cface_onP => this%cface(:this%xcface(this%ncell_onP+1)-1))
      call get_global_ragged_array (xcface_onP, this%face_ip%global_index(cface_onP), xcface, cface)
    end associate
  end subroutine get_global_cface_array

  !! Auxiliary subroutine creates a global ragged array on the IO process,
  !! 0-sized on others, given a distributed ragged array.
  subroutine get_global_ragged_array (xarray_l, array_l, xarray, array)
    use parallel_communication, only: nPE, is_IOP, global_sum, collate, distribute
    integer, intent(in) :: xarray_l(:), array_l(:)
    integer, allocatable, intent(out) :: xarray(:), array(:)
    integer :: offset
    ASSERT(size(xarray_l) >= 1)
    ASSERT(xarray_l(1) == 1)
    ASSERT(size(array_l) == xarray_l(size(xarray_l))-1)
    ASSERT(all(xarray_l(2:) - xarray_l(:size(xarray_l)-1) >= 0))
#ifdef NAG_COMPILER_WORKAROUND
    !! NAGFOR 6.0(1052)/GCC 4.8.3 produce bad code under optimization (-O2)
    !! with the one-line allocation in parallel, getting the wrong value for
    !! the global_sum expression, apparently it gets the value of its argument
    !! instead.  Looking at the intermediate C code suggests that tha problem
    !! may not be here but in global_sum and perhaps a race condition?
    offset = global_sum(size(xarray_l)-1)
    allocate(xarray(1+merge(offset,0,is_IOP)))
#else
    allocate(xarray(1+merge(global_sum(size(xarray_l)-1),0,is_IOP)))
#endif
    offset = excl_prefix_sum(size(array_l))
    xarray(1) = 1
    call collate (xarray(2:), xarray_l(2:)+offset)
#ifdef NAG_COMPILER_WORKAROUND
    !! Same comments as above.
    offset = global_sum(size(array_l))
    allocate(array(merge(offset,0,is_IOP)))
#else
    allocate(array(merge(global_sum(size(array_l)),0,is_IOP)))
#endif
    call collate (array, array_l)
    if (is_IOP) then
      ASSERT(size(xarray) >= 1)
      ASSERT(xarray(1) == 1)
      ASSERT(all(xarray(2:) - xarray(:size(xarray)-1) >= 0))
      ASSERT(size(array) == xarray(size(xarray))-1)
    end if
  contains
    integer function excl_prefix_sum (n) result (psum)
      integer, intent(in) :: n
      integer :: j
      integer, allocatable :: array(:)
      allocate(array(merge(nPE,0,is_IOP)))
      call collate (array, n)
      if (is_IOP) then
        do j = 2, nPE
          array(j) = array(j) + array(j-1)
        end do
      end if
      call distribute (psum, array)
      psum = psum - n
    end function
  end subroutine get_global_ragged_array

 !! Writes to the tty and output file a profile of the distributed mesh:
 !! numbers of nodes, faces, and cells assigned to each processor; numbers
 !! of on-process and off-process objects.

  subroutine write_profile (this)

    use parallel_communication, only: nPE, broadcast, collate
    use truchas_logging_services

    class(unstr_mesh), intent(in) :: this

    integer :: n
    character(80) :: line
    integer, dimension(nPE) :: nnode_vec, nface_vec, ncell_vec
    integer, dimension(2,nPE) :: nvec, fvec, cvec

    call collate (nnode_vec, this%nnode)
    call collate (nface_vec, this%nface)
    call collate (ncell_vec, this%ncell)

    call broadcast (nnode_vec)
    call broadcast (nface_vec)
    call broadcast (ncell_vec)

    call TLS_info ('  UNSTR_MESH Profile:')
    write(line,fmt='(4x,a3,a,4a9)') 'PE', '|', 'nnode', 'nface', 'ncell'
    call TLS_info (line)
    call TLS_info ('    ---+'//repeat('-',27))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,3i9)') n, '|', nnode_vec(n),  nface_vec(n), ncell_vec(n)
      call TLS_info (line)
    end do

    call collate (nvec(1,:), this%node_ip%offP_size())
    call collate (nvec(2,:), this%node_ip%onP_size())
    call broadcast (nvec)

    call collate (fvec(1,:), this%face_ip%offP_size())
    call collate (fvec(2,:), this%face_ip%onP_size())
    call broadcast (fvec)

    call collate (cvec(1,:), this%cell_ip%offP_size())
    call collate (cvec(2,:), this%cell_ip%onP_size())
    call broadcast (cvec)

    call TLS_info ('  Mesh Communication Profile:')
    write(line,fmt='(4x,3x,1x,a11,2a16)')  'Nodes', 'Faces', 'Cells'
    call TLS_info (line)
    write(line,fmt='(4x,a3,a,3a16)') 'PE', '|', ('off-PE   on-PE', n=1,3)
    call TLS_info (line)
    call TLS_info ('    ---+'//repeat('-',48))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,3(i7,i9))') n, '|', nvec(:,n), fvec(:,n), cvec(:,n)
      call TLS_info (line)
    end do

  end subroutine write_profile

  !! This looks for boundary faces that do not belong to any face set, nor any
  !! link. If any are found, a warning message is written.  In such cases it
  !! may not be possible to specify a complete set of boundary conditions.

  subroutine check_bndry_face_set(this)

    use string_utilities, only: i_to_c
    use truchas_logging_services

    class(unstr_mesh), intent(in) :: this

    integer, parameter :: MAX_PRINT = 10
    integer :: j, n, nqf, array(nPE)
    logical, allocatable :: mask(:)
    real(r8), allocatable :: xc(:,:), xc_l(:,:)
    character(70), allocatable :: msg(:)
    type(bitfield) :: bitmask

    !! Boundary faces that do not belong to any face set.
    bitmask = ibset(ZERO_BITFIELD, pos=0)
    mask = (this%face_set_mask == bitmask)

    !! These boundary faces are associated with links and are okay.
    mask(this%lface(1,:)) = .false.
    mask(this%lface(2,:)) = .false.

    !! Count the number of questionable boundary faces (on-process).
    nqf = count(mask(:this%nface_onP))
    call collate(array, nqf)
    if (is_IOP) nqf = sum(array)
    call broadcast(nqf)

    if (nqf == 0) return  ! nothing to see here

    !! Distribute the number of faces each process will supply data for.
    if (is_IOP) then
      n = MAX_PRINT
      do j = 1, nPE
        array(j) = min(array(j), n)
        n = n - array(j)
      end do
      allocate(xc(3,sum(array)))
    else
      allocate(xc(3,0))
    end if
    call distribute(n, array)
    allocate(xc_l(3,n))

    !! Collect the face centroid data.
    n = 0
    do j = 1, this%nface_onP
      if (.not.mask(j)) cycle
      n = n + 1
      if (n > size(xc_l,dim=2)) exit
      associate (fnode => this%fnode(this%xfnode(j):this%xfnode(j+1)-1))
        xc_l(:,n) = sum(this%x(:,fnode),dim=2) / size(fnode)
      end associate
    end do
    call collate(xc, xc_l)

    !! Write the warning message.
    allocate(msg(3+size(xc,dim=2)))
    msg(1) = i_to_c(nqf) // ' boundary faces do not belong to any face set or interface.'
    msg(2) = 'This may make it impossible to specify a complete set of BC.'
    msg(3) = 'Some of the face centroids are'
    do j = 1, size(xc,dim=2)
      write(msg(3+j),'("(",es13.5,2(",",es13.5)," )")') xc(:,j)
    end do
    call TLS_info(repeat('*',79))
    call TLS_warn(msg)
    call TLS_info(repeat('*',79))

  end subroutine check_bndry_face_set

  subroutine get_link_set_ids(this, mask, setids)

    class(unstr_mesh), intent(in) :: this
    logical, intent(in) :: mask(:)
    integer, allocatable, intent(out) :: setids(:)

    integer :: j
    type(bitfield) :: bitmask

    ASSERT(size(mask) == this%nface)

    bitmask = ZERO_BITFIELD
    do j = 1, this%nlink_onP
      if (any(mask(this%lface(:,j)))) bitmask = ior(bitmask, this%link_set_mask(j))
    end do
    bitmask = global_ior(bitmask)

    setids = pack(this%link_set_id, mask=btest(bitmask, pos=[(j,j=1,size(this%link_set_id))]))

  end subroutine get_link_set_ids

  !! This function identifies the global cell that is nearest the given point,
  !! and returns its local index on the process that owns the cell, and 0 on
  !! all other processes. Here nearest means the cell whose centroid is nearest
  !! the given point. Least cell index breaks ties. Note that this means the
  !! cell may not actually contain the point at all. The immediate purpose is
  !! in solution probe initialization, and it may not be suitable for anything
  !! else.

  function nearest_cell(this, point)

    class(unstr_mesh), intent(in) :: this
    real(r8), intent(in) :: point(:)
    integer :: nearest_cell

    integer :: j, min_cell, min_PE
    real(r8) :: min_dist, d, array(nPE), centroid(size(this%x,dim=1))

    ASSERT(size(point) == size(this%x,dim=1))

    !! Compute the minimum distance and cell index for the local mesh subdomain
    min_dist = huge(min_dist)
    do j = 1, this%ncell_onP
      associate(cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1))
        centroid = sum(this%x(:,cell_nodes),dim=2)/size(cell_nodes)
      end associate
      d = norm2(centroid-point)
      if (d < min_dist) then
        min_dist = d
        min_cell = j
      end if
    end do

    !! Determine the nearest cell and its owner globally
    call collate(array, min_dist)
    if (is_IOP) min_PE = minloc(array,dim=1)
    call broadcast(min_PE)
    nearest_cell = merge(min_cell, 0, (this_PE == min_PE))

  end function nearest_cell

  !! This function identifies the global node that is nearest the given point,
  !! and returns its local index on the process that owns the node, and 0 on
  !! all other processes. Least node index breaks ties.  The immediate purpose
  !! is in solution probe initialization, and it may not be suitable for
  !! anything else.

  function nearest_node(this, point)

    class(unstr_mesh), intent(in) :: this
    real(r8), intent(in) :: point(:)
    integer :: nearest_node

    integer :: j, min_node, min_PE
    real(r8) :: min_dist, d, array(nPE)

    ASSERT(size(point) == size(this%x,dim=1))

    !! Compute the minimum distance and node index for the local mesh subdomain.
    min_dist = huge(min_dist)
    do j = 1, this%nnode_onP
      d = norm2(this%x(:,j)-point)
      if (d < min_dist) then
        min_dist = d
        min_node = j
      end if
    end do

    !! Determine the nearest node and its owner globally.
    call collate(array, min_dist)
    if (is_IOP) min_PE = minloc(array,dim=1)
    call broadcast(min_PE)
    nearest_node = merge(min_node, 0, (this_PE == min_PE))

  end function nearest_node

end module unstr_mesh_type
