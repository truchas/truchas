!!
!! BNDRY_NODE_GROUP_BUILDER_TYPE
!!
!! This module defines an auxiliary object that constructs a grouped list of
!! mesh nodes specified incrementally using face set IDs. Its principal use is
!! in the instantiation of boundary condition objects.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! March 2025
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bndry_node_group_builder_type

  use simpl_mesh_type
  use bndry_face_group_builder_type
  implicit none
  private

  type, public :: bndry_node_group_builder
    private
    type(simpl_mesh), pointer :: mesh => null() ! reference only
    type(bndry_face_group_builder) :: builder
    logical :: no_overlap = .true.
  contains
    procedure :: init
    procedure :: add_face_group
    procedure :: get_face_groups
    procedure :: get_node_groups
  end type

contains

  subroutine init(this, mesh, no_overlap)
    class(bndry_node_group_builder), intent(out) :: this
    type(simpl_mesh), target :: mesh
    logical, intent(in), optional :: no_overlap
    this%mesh => mesh
    if (present(no_overlap)) this%no_overlap = no_overlap
    call this%builder%init(mesh, bndry_only=.true., omit_offp=.false., no_overlap=this%no_overlap)
  end subroutine

  subroutine add_face_group(this, setids, stat, errmsg)
    class(bndry_node_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:) ! NB: face set IDs
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
  end subroutine

  subroutine get_face_groups(this, ngroup, xgroup, faces, fnode, nodes)

    use integer_map_type

    class(bndry_node_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), faces(:), fnode(:,:), nodes(:)

    integer :: i, j, n
    type(integer_map) :: node_map

    call this%builder%get_face_groups(ngroup, xgroup, faces)

    !! Tag nodes contained in the specified boundary faces.
    do j = 1, size(faces)
      do i = 1, size(this%mesh%fnode,dim=1)
        call node_map%set(this%mesh%fnode(i,faces(j)), 0)
      end do
    end do

    n = node_map%size()
    allocate(nodes(n))

    !! Generate the list of tagged nodes and the inverse mapping of
    !! nodes to tagged nodes.
    n = 0
    do j = 1, this%mesh%nnode
      if (node_map%contains(j)) then
        n = n + 1
        nodes(n) = j
        call node_map%set(j, n)
      end if
    end do

    !! Extract the section of the mesh%fnode array for the specified boundary
    !! faces and re-index the values to point to the tagged nodes.
    allocate(fnode(size(this%mesh%fnode,dim=1),size(faces)))
    do j = 1, size(faces)
      do i = 1, size(this%mesh%fnode,dim=1)
        fnode(i,j) = node_map%val(this%mesh%fnode(i,faces(j)))
      end do
    end do

  end subroutine get_face_groups

  subroutine get_node_groups(this, ngroup, xgroup, index, stat)

    use parallel_communication, only: global_any

    class(bndry_node_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:)
    integer, intent(out) :: stat

    integer :: i, j, n, gsize
    logical, allocatable :: gmask(:), mask(:)

    type :: array_box
      integer, allocatable :: array(:)
    end type
    type(array_box), allocatable :: glist(:)

    stat = 0

    call this%builder%get_face_groups(ngroup, xgroup, index)

    allocate(gmask(this%mesh%nnode), glist(ngroup))
    do i = 1, ngroup
      !! Tag all nodes that belong to a face in the group
      gmask = .false.
      associate (gface => index(xgroup(i):xgroup(i+1)-1))
        do j = 1, size(gface)
          associate (fnode => this%mesh%fnode(:,gface(j)))
            gmask(fnode) = .true.
          end associate
        end do
      end associate
      call this%mesh%node_imap%scatter_offp_or(gmask)

      !! Generate the list of node indices
      gsize = count(gmask)
      allocate(glist(i)%array(gsize))
      n = 0
      do j = 1, size(gmask)
        if (gmask(j)) then
          n = n + 1
          glist(i)%array(n) = j
        end if
      end do

      !! Check for overlapping groups if requested.
      if (this%no_overlap) then
        if (i == 1) then
          if (i < ngroup) mask = gmask ! save for next pass
        else
          if (global_any(mask(glist(i)%array))) stat = 1
          mask(glist(i)%array) = .true.
        end if
      end if
    end do

    n = 0
    do i = 1, ngroup
      n = n + size(glist(i)%array)
    end do

    deallocate(index)
    allocate(index(n))
    xgroup(1) = 1
    do i = 1, ngroup
      xgroup(i+1) = xgroup(i) + size(glist(i)%array)
      index(xgroup(i):xgroup(i+1)-1) = glist(i)%array
    end do

  end subroutine get_node_groups

end module bndry_node_group_builder_type
