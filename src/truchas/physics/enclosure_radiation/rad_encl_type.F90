!!
!! RAD_ENCL_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for Fortran 2008, May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include <f90_assert.fpp>

module rad_encl_type

  use kinds
  use index_map_type
  implicit none
  private

  type, public :: rad_encl
    !! Local surface mesh.
    integer :: nnode = 0, nface = 0
    integer, allocatable :: xface(:), fnode(:)
    real(r8), allocatable :: coord(:,:)
    !! Mappings to external numberings.
    integer, allocatable :: node_map(:), face_map(:)
    !! Face block data.
    integer, allocatable :: face_block_id(:), face_block(:)
    !! Partitioning and inter-process communication data.
    integer :: nnode_onP = 0, nface_onP = 0
    type(index_map) :: node_ip, face_ip
  contains
    procedure :: init
  end type rad_encl

contains

  subroutine init (this, file, fmap, csf)

    use rad_encl_file_type
    use permutations
    use parallel_communication

    class(rad_encl), intent(out) :: this
    type(rad_encl_file), intent(in) :: file
    integer, intent(in) :: fmap(:)
    real(r8), intent(in), optional :: csf  ! Coordinate scaling factor

    integer :: j, n, nnode, nface
    integer, allocatable :: fsize(:), fnode(:), group_ids(:), gnum(:)
    integer, allocatable :: node_map(:), face_map(:), fsize_l(:)
    real(r8), allocatable :: coord(:,:)
    integer :: face_bsize(nPE), node_bsize(nPE)

    !! Read the enclosure data from the file
    if (is_IOP) then

      !! Read the surface mesh.
      call file%get_surface(fsize, fnode, coord)
      nface = size(fsize)
      nnode = size(coord,dim=2)

      !! Read the face block data.
      call file%get_group_info(gnum, group_ids)

    else  ! allocate dummy arrays
      nface = 0
      nnode = 0
      allocate (fsize(0), fnode(0), coord(3,0), gnum(0))
    end if

    !! Reorder and block partition the global surface mesh.
    allocate(node_map(nnode), face_map(nface))
    call collate(face_map, fmap)
    call collate(face_bsize, size(fmap))
    if (is_IOP) then
      !! Reorder the faces; reorder the face block index array.
      call organize_faces (face_map, fsize, fnode)
      call reorder (gnum, face_map)
      !! Partition and reorder the nodes; reorder the node coordinate array.
      call organize_nodes (fsize, fnode, face_bsize, node_bsize, node_map)
      call reorder (coord, node_map)
    end if

    !! Create the face index partition; no off-process (ghost) faces.
    call this%face_ip%init (face_bsize)

    !! Create the node index partition; no off-process (ghost) nodes yet.
    call this%node_ip%init (node_bsize)

    !! Localize the global face connectivity structure arrays.  This identifies
    !! off-process nodes that need to be added as ghost nodes on this process.
    call this%face_ip%localize_index_array(fsize, fnode, this%node_ip, fsize_l, this%fnode)
    deallocate(fsize, fnode)

    !! Generate the local face indexing array from the local face sizes.
    allocate(this%xface(1+size(fsize_l)))
    this%xface(1) = 1
    do j = 1, size(fsize_l)
      this%xface(j+1) = this%xface(j) + fsize_l(j)
    end do
    deallocate(fsize_l)

    !! Local surface mesh sizes: on-process plus off-process.
    this%nnode = this%node_ip%local_size
    this%nface = this%face_ip%local_size

    !! On-process surface mesh sizes
    this%nnode_onP = this%node_ip%onp_size
    this%nface_onP = this%face_ip%onp_size

    !! Distribute the node coordinate array.
    allocate(this%coord(3,this%nnode))
    call this%node_ip%distribute(coord, this%coord)
    call this%node_ip%gather_offp(this%coord)
    deallocate(coord)

    !! Distribute the node map array.
    allocate(this%node_map(this%nnode))
    call this%node_ip%distribute(node_map, this%node_map)
    call this%node_ip%gather_offp(this%node_map)
    deallocate(node_map)

    !! Distribute the face map array.
    allocate(this%face_map(this%nface))
    call this%face_ip%distribute(face_map, this%face_map)
    !call this%face_ip%gather_offp(this%face_map)
    deallocate(face_map)

    !! Replicate the face block ID array.
    if (is_IOP) n = size(group_ids)
    call broadcast (n)
    allocate(this%face_block_id(n))
    if (is_IOP) this%face_block_id = group_ids
    call broadcast (this%face_block_id)
    if (is_IOP) deallocate(group_ids)

    !! Distribute the face block index array.
    allocate(this%face_block(this%nface))
    call this%face_ip%distribute(gnum, this%face_block)
    !call this%face_ip%gather_offp(this%face_block)
    deallocate(gnum)

    !! Scale coordinates
    if (present(csf) .and. csf /= 1.0_r8) this%coord = csf * this%coord

  end subroutine init

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ORGANIZE_FACES
 !!
 !! This auxillary routine organizes the surface mesh faces in preparation
 !! for creating the distributed enclosure data structure.  This involves
 !! generating a new face numbering in which the specified coloring is a
 !! blocked coloring and the faces within the same-colored blocks are
 !! well-ordered.  The renumbering is returned in FACE_MAP (new-to-old) and
 !! the block sizes in FACE_BSIZE.
 !!
 !! The face connectivity structure arrays FSIZE and FNODE are reordered to
 !! reflect the new numbering (FCOLOR is not), but the caller must ensure that
 !! all other face-based arrays are reordered using the returned FACE_MAP, and
 !! the caller must also remap any face-valued data using its inverse.
 !!
 !! For faces belonging to a same-colored block, well-ordered usually means
 !! an order that gives good memory locality: adjacent faces, which share nodes,
 !! should be numbered close together.  For our uses of an enclosure this may
 !! not be particularly relevant.  We do nothing here, but this is where it
 !! should be done.
 !!

  subroutine organize_faces (face_map, fsize, fnode)

    use permutations

    integer, intent(in)    :: face_map(:)         ! new-to-old face map
    integer, intent(inout) :: fsize(:), fnode(:)  ! face connectivity

    integer :: j, k, first, xface(1+size(fsize)), fnode_new(size(fnode))

    ASSERT(all(fsize >= 0))
    ASSERT(sum(fsize) == size(fnode))
    ASSERT(size(fsize) == size(face_map))

    !! Generate face indexing into FNODE from the face sizes:
    !! FNODE(XFACE(j):XFACE(j+1)-1) are the nodes of face j.
    xface(1) = 1
    do j = 1, size(fsize)
      xface(j+1) = xface(j) + fsize(j)
    end do

    !! Reorder face sizes array.
    call reorder (fsize, face_map)

    !! Reorder the face node array.
    first = 1
    do j = 1, size(fsize)
      do k = 0, fsize(j)-1
        fnode_new(k+first) = fnode(k+xface(face_map(j)))
      end do
      first = first + fsize(j)
    end do
    fnode = fnode_new

  end subroutine organize_faces

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ORGANIZE_NODES
 !!
 !! This auxillary routine organizes the surface mesh nodes in preparation for
 !! creating the distributed enclosure data structure.  This involves coloring
 !! the nodes and generating a new node numbering in which the coloring is a
 !! blocked coloring and the nodes within same-colored blocks are well-ordered.
 !! The computed renumbering is returned in NODE_MAP (new-to-old) and the block
 !! sizes in NODE_BSIZE.
 !!
 !! The values of the face connectivity structure array FNODE are mapped to the
 !! new numbering, but the caller must ensure that all node-based arrays are
 !! reordered using the returned NODE_MAP, and the caller must also map any
 !! other node-valued data using its inverse.
 !!
 !! It is assumed that the faces of the surface mesh, described by the face
 !! connectivity arrays FSIZE and FNODE, have already been organized by a call
 !! to ORGANIZE_FACES, and that all nodes are referenced by the connectivity
 !! data.
 !!
 !! The coloring of the nodes is based on the specifed blocked coloring of the
 !! faces given by FACE_BSIZE.  A node is colored like one of the faces that
 !! contain it.  We use a simple greedy algorithm: for colors n = 1, 2, ...,
 !! color all uncolored nodes that belong to an n-colored cell with color n.
 !! This may lead to a poorly balanced coloring, but this can be addressed
 !! later if it becomes a significant issue.
 !!
 !! To order the nodes in a same-colored block we rely on the well-orderedness
 !! of the faces.  We merely number them as they are encountered in the FNODE
 !! array.
 !!

  subroutine organize_nodes (fsize, fnode, face_bsize, node_bsize, node_map)

    use permutations

    integer, intent(in)    :: fsize(:)      ! face sizes
    integer, intent(inout) :: fnode(:)      ! packed face connectivities
    integer, intent(in)    :: face_bsize(:) ! face partition block sizes
    integer, intent(out)   :: node_bsize(:) ! node partition block sizes
    integer, intent(out)   :: node_map(:)   ! new-to-old node map

    integer :: i, j, n, ioff, joff, pass(size(node_map)), map(size(node_map))

    ASSERT(all(fsize >= 0))
    ASSERT(sum(fsize) == size(fnode))
    ASSERT(size(face_bsize) == size(node_bsize))
    ASSERT(minval(fnode) == 1 .and. maxval(fnode) == size(node_map))

    !! Generate a good ordering of the nodes: number the nodes consecutively
    !! as they are encountered in FNODE.  NODE_MAP is old-to-good numbering.
    n = 0
    node_map = 0
    do j = 1, size(fnode)
      if (node_map(fnode(j)) /= 0) cycle ! numbered this one already
      n = n + 1
      node_map(fnode(j)) = n
    end do
    ASSERT(is_perm(node_map))

    call invert_perm (node_map)  ! NODE_MAP is now good-to-old numbering

    !! Partition the nodes: a simple greedy algorithm.
    pass = 0  ! partition assignment
    ioff = 0  ! offset into the FSIZE array
    joff = 0  ! offset into the FNODE array
    do n = 1, size(face_bsize)
      do i = ioff+1, ioff+face_bsize(n)
        do j = joff+1, joff+fsize(i)
          if (pass(fnode(j)) == 0) pass(fnode(j)) = n
        end do
        joff = joff + fsize(i)
      end do
      ioff = ioff + face_bsize(n)
    end do
    ASSERT(all(pass > 0))

    !! Partition assignment relative to the good ordering.
    call reorder (pass, node_map)

    !! Get the node mapping that makes the partition blocked (new-to-good).
    call blocked_coloring_map (pass, node_bsize, map)

    !! Total node mapping (new-to-old).
    call reorder (node_map, map)
    ASSERT(is_perm(node_map))

    !! Map the values of the FNODE array to new values.
    call invert_perm (node_map, map)  ! MAP is old-to-new
    do j = 1, size(fnode)
      fnode(j) = map(fnode(j))
    end do

  end subroutine organize_nodes

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BLOCKED_COLORING_MAP
 !!
 !! Given a coloring vector COLOR, this auxillary routine computes a mapping
 !! (or renumbering) that makes the coloring a blocked coloring, and the sizes
 !! of the resulting blocks.  That is, the array COLOR(MAP(:)) is such that
 !! the first BSIZE(1) elements are color 1, the next BSIZE(2) elements are
 !! color 2, and so on.  The mapping preserves the relative order of elements
 !! having the same color.
 !!

  subroutine blocked_coloring_map (color, bsize, map)

    use permutations

    integer, intent(in)  :: color(:)  ! coloring
    integer, intent(out) :: bsize(:)  ! block size
    integer, intent(out) :: map(:)    ! mapping permutation

    integer :: j, n, next(size(bsize))

    ASSERT(size(color) == size(map))
    ASSERT(minval(color) >= 1 .and. maxval(color) <= size(bsize))

    !! Compute the block size of each color.
    bsize = 0
    do j = 1, size(color)
      bsize(color(j)) = bsize(color(j)) + 1
    end do

    !! NEXT(j) is the next free index for block j.
    next(1) = 1
    do n = 2, size(bsize)
      next(n) = next(n-1) + bsize(n-1)
    end do

    !! Generate the mapping (new-to-old)
    do j = 1, size(color)
      map(next(color(j))) = j
      next(color(j)) = next(color(j)) + 1
    end do
    ASSERT(is_perm(map))

  end subroutine blocked_coloring_map

end module rad_encl_type
