!!
!! UNSTR_2D_MESH_TYPE
!!
!! This module provides a derived type that encapsulates the data describing
!! a distributed unstructured 2D mixed-element mesh. It formally supports
!! arbitrary polygons, but likely triangles and quadrilaterals in practice.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type UNSTR_2D_MESH that encapsulates the
!!  data describing a distributed unstructured 2D mixed-element mesh. On each
!!  process the object describes a complete mesh of some subdomain that only
!!  references local entities (cells, faces, nodes), and can rightly be
!!  considered as a serial mesh for that subdomain. Globally, the subdomains
!!  will overlap, perhaps only along a boundary but more generally on some
!!  collection of cells.  Additional mesh data describes this overlap and
!!  provides for communication between overlapping entities.
!!
!!  Objects of this type are intended to be used by trusted code, and so its
!!  data components are public.  However, the components must be treated as
!!  read-only because an object may be shared amongst multiple clients. The
!!  following data components are accessible:
!!
!!    nnode, nface, ncell - the number of nodes, faces, and cells in the
!!        subdomain mesh.
!!
!!    cstart -- rank-1 integer array giving the array location where the cell
!!        connectivity data for each cell starts.
!!
!!    cnode -- rank-1 integer array storing the cell node connectivity data:
!!        cnode(cstart(j):cstart(j+1)-1) is the list of node indices that
!!        define cell j.
!!
!!    cface -- rank-1 integer array storing the cell face connectivity data:
!!        cface(cstart(j):cstart(j+1)-1) is the list of faces belonging to
!!        cell j in canonical order.
!!
!!    cnhbr -- rank-1 integer array storing the cell neighbor connectivity
!!        data: cnhbr(cstart(j):cstart(j+1)-1) is the list of face neighbors
!!        of cell j, and corresponds to its list of faces.
!!
!!    cfpar - an integer bit mask array storing the relative cell face
!!        orientations: btest(cfpar(j),k) is true when face k of cell j is
!!        inward oriented with respect to cell j, and false when it is
!!        outward oriented.
!!
!!    fnode - rank-2 integer array storing the two nodes defining each face:
!!        fnode(:,j) is the ordered pair of node indices defining the oriented
!!        face j.
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
!!        R^2 of node j.  Its shape is [2,nnode].
!!
!!    area - the rank-1 real array of face areas; area(j) is the area of face j.
!!        Its shape is [nface].
!!
!!    volume - the rank-1 real array of cell volumes; volume(j) is the volume
!!        of cell j.  Its shape is [ncell].
!!
!!    normal - the rank-2 real array of oriented face areas; normal(:,j) is the
!!        oriented area of face j.  Its shape is [2,nface].
!!
!!    unit-normal - the rank-2 real array of face unit-normals; unit_normal(:,j) is
!!        the unit-normal of face j.  Its shape is [2,nface].
!!        The unit-normal is required to be stored separately due to the fact that the
!!        axis of symmetry has a zero area. This would reduce its area-weighted normal
!!        to zero. The normal to the axis of symmetry is needed for axisymmetric
!!        VOF calculations.
!!

#include "f90_assert.fpp"

module unstr_2d_mesh_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_base_mesh_class
  use index_partitioning
  use parallel_communication
  use bitfield_type
  use cell_topology
  use f08_intrinsics
  implicit none
  private

  type, extends(unstr_base_mesh), public :: unstr_2d_mesh
    integer, allocatable :: cnode(:) ! cell nodes connectivity
    integer, allocatable :: cface(:) ! cell faces connectivity
    integer, allocatable :: cnhbr(:) ! cell neighbors connectivity
    integer, allocatable :: cstart(:) ! start array index for cell connectivity
    integer, allocatable :: fnode(:,:) ! face nodes connectivity
    integer, allocatable :: cfpar(:)  ! relative cell face orientation (bit mask)
    integer, allocatable :: fcell(:,:)  ! face cell neighbors
    real(r8), allocatable :: normal(:,:)
    real(r8), allocatable :: unit_normal(:,:)
    real(r8), allocatable :: cell_centroid(:,:)
    real(r8), allocatable :: face_centroid(:,:)
  contains
    procedure :: get_global_cnode_array
    procedure :: get_global_cface_array
    procedure :: compute_geometry
    procedure :: write_profile
    procedure :: check_bndry_face_set
    procedure :: get_link_set_bitmask
    procedure :: get_link_set_ids
    procedure :: init_cell_centroid
    procedure :: init_face_centroid
    procedure :: nearest_node
    procedure :: nearest_cell
    procedure :: cell_node_list_view
    procedure :: cell_face_list_view
    procedure :: face_node_list_view
  end type unstr_2d_mesh

contains

  !! Compute the geometric data components from the node coordinates.
  subroutine compute_geometry (this)
    use cell_geometry, only: cell_volume, vector_length
    class(unstr_2d_mesh), intent(inout) :: this
    integer :: j
    ASSERT(allocated(this%volume))
    ASSERT(allocated(this%normal))
    ASSERT(allocated(this%area))
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%cstart(j):this%cstart(j+1)-1))
        this%volume(j) = cell_volume(this%x(:,cell_nodes))
      end associate
    end do
    do j = 1, this%nface
      this%normal(1,j) = this%x(2,this%fnode(2,j)) - this%x(2,this%fnode(1,j))
      this%normal(2,j) = this%x(1,this%fnode(1,j)) - this%x(1,this%fnode(2,j))
      this%area(j) = vector_length(this%normal(:,j))
      this%unit_normal(1,j) = this%normal(1,j) / this%area(j)
      this%unit_normal(2,j) = this%normal(2,j) / this%area(j)
    end do
  end subroutine compute_geometry

  subroutine init_cell_centroid(this)
    use cell_geometry, only: cell_centroid_2d
    class(unstr_2d_mesh), intent(inout) :: this
    integer :: j
    if (allocated(this%cell_centroid)) return
    allocate(this%cell_centroid(2,this%ncell))
    do j = 1, this%ncell
      associate(cell_nodes => this%cnode(this%cstart(j):this%cstart(j+1)-1))
        this%cell_centroid(:,j) = cell_centroid_2d(this%x(:,cell_nodes))
      end associate
    end do
  end subroutine init_cell_centroid

  subroutine init_face_centroid(this)
    class(unstr_2d_mesh), intent(inout) :: this
    integer :: j
    if (allocated(this%face_centroid)) return
    allocate(this%face_centroid(2,this%nface))
    do j = 1, this%nface
      this%face_centroid(:,j) = 0.5_r8 * sum(this%x(:,this%fnode(:,j)), dim=2)
    end do
  end subroutine init_face_centroid

  !! Creates the global ragged CNODE array on the IO process, 0-sized on others.
  subroutine get_global_cnode_array (this, cstart, cnode)
    class(unstr_2d_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cstart(:), cnode(:)
    associate (cstart_onP => this%cstart(:this%ncell_onP+1), &
                cnode_onP => this%cnode(:this%cstart(this%ncell_onP+1)-1))
      call get_global_ragged_array (cstart_onP, this%node_ip%global_index(cnode_onP), cstart, cnode)
    end associate
  end subroutine get_global_cnode_array

  !! Creates the global ragged CFACE array on the IO process, 0-sized on others.
  subroutine get_global_cface_array (this, cstart, cface)
    class(unstr_2d_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: cstart(:), cface(:)
    associate (cstart_onP => this%cstart(:this%ncell_onP+1), &
                cface_onP => this%cface(:this%cstart(this%ncell_onP+1)-1))
      call get_global_ragged_array (cstart_onP, this%face_ip%global_index(cface_onP), cstart, cface)
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

    class(unstr_2d_mesh), intent(in) :: this

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

    call TLS_info ('  unstr_2d_mesh Profile:')
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

    class(unstr_2d_mesh), intent(in) :: this

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
      xc_l(:,n) = sum(this%x(:,this%fnode(:,j))) / 2
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

  !! Returns a scalar bit mask for use in bit operations with the link_set_mask
  !! array component.  The corresponding bit is set for each link set ID given
  !! in the array SETIDS.  STAT returns a non-zero value if an unknown link set
  !! ID is specified, and the optional allocatable deferred-length character
  !! ERRMSG is assigned an explanatory message if present.

  subroutine get_link_set_bitmask (this, setids, bitmask, stat, errmsg)
    use string_utilities, only: i_to_c
    class(unstr_2d_mesh), intent(in) :: this
    integer, intent(in) :: setids(:)
    type(bitfield), intent(out) :: bitmask
    integer, intent(out) :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: i, j
    bitmask = ZERO_BITFIELD
    do i = 1, size(setids)
      do j = size(this%link_set_id), 1, -1
        if (setids(i) == this%link_set_id(j)) exit
      end do
      if (j == 0) then
        stat = 1
        if (present(errmsg)) errmsg = 'unknown link set ID: ' // i_to_c(setids(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do
    stat = 0
  end subroutine get_link_set_bitmask

  subroutine get_link_set_ids(this, mask, setids)

    class(unstr_2d_mesh), intent(in) :: this
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

    class(unstr_2d_mesh), intent(in) :: this
    real(r8), intent(in) :: point(:)
    integer :: nearest_cell

    integer :: j, min_cell, min_PE
    real(r8) :: min_dist, d, array(nPE), centroid(size(this%x,dim=1))

    ASSERT(size(point) == size(this%x,dim=1))

    !! Compute the minimum distance and cell index for the local mesh subdomain
    min_dist = huge(min_dist)
    do j = 1, this%ncell_onP
      associate(cell_nodes => this%cnode(this%cstart(j):this%cstart(j+1)-1))
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

    class(unstr_2d_mesh), intent(in) :: this
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

  !! Returns the nodes of the given cell
  function cell_node_list_view(this, n) result(view)
    class(unstr_2d_mesh), intent(in), target :: this
    integer, intent(in) :: n
    integer, pointer, contiguous :: view(:)
    view => this%cnode(this%cstart(n):this%cstart(n+1)-1)
  end function

  !! Returns the faces of the given cell
  function cell_face_list_view(this, n) result(view)
    class(unstr_2d_mesh), intent(in), target :: this
    integer, intent(in) :: n
    integer, pointer, contiguous :: view(:)
    view => this%cface(this%cstart(n):this%cstart(n+1)-1)
  end function

  !! Returns the nodes of the given face
  function face_node_list_view(this, n) result(view)
    class(unstr_2d_mesh), intent(in), target :: this
    integer, intent(in) :: n
    integer, pointer, contiguous :: view(:)
    view => this%fnode(:,n)
  end function

end module unstr_2d_mesh_type
