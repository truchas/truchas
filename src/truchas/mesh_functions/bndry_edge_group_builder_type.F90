!!
!! BNDRY_EDGE_GROUP_BUILDER_TYPE
!!
!! This module defines a builder that incrementally converts boundary face-set
!! selections into grouped edge-index representations. It is primarily used to
!! construct edge-based electromagnetic boundary functions.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, January 2024, modified July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The builder retains an unowned reference to the mesh. After initialization,
!! zero or more boundary-face groups are added using face-set IDs. Groups are
!! numbered in insertion order. All result arrays are allocated by the getter
!! procedures.
!!
!! Face overlap and derived-edge overlap are independent policies, both rejected
!! by default. Distinct, nonoverlapping face groups may legitimately share
!! edges. When edge-overlap checking is enabled, overlap is reported through a
!! status value but does not prevent construction of the edge groups.
!!
!! Two edge representations are available. The direct representation
!! concatenates the groups and may repeat edges across groups. The mapped
!! representation instead returns one ordered, unique edge array and a grouped
!! map into it. The original face groups are also available directly.
!!

#include "f90_assert.fpp"

module bndry_edge_group_builder_type

  use simpl_mesh_type
  use bndry_face_group_builder_type
  implicit none
  private

  type :: array_box
    integer, allocatable :: array(:)
  end type

  type, public :: bndry_edge_group_builder
    private
    type(simpl_mesh), pointer :: mesh => null() ! reference only -- not owned
    type(bndry_face_group_builder) :: builder
    logical :: no_face_overlap = .true., no_edge_overlap = .true.
  contains
    procedure :: init
    procedure :: add_face_group
    procedure :: get_face_groups
    procedure :: get_edge_groups
    procedure :: get_unique_edge_groups
    procedure, private :: build_edge_groups
  end type

contains

  !! Initialize an empty builder for MESH. Overlap between face groups and
  !! overlap between their derived edge groups are rejected by default. The
  !! optional arguments independently override those policies. Added faces are
  !! required to be boundary faces, and off-process faces are retained.
  subroutine init(this, mesh, no_face_overlap, no_edge_overlap)
    class(bndry_edge_group_builder), intent(out) :: this
    type(simpl_mesh), target :: mesh
    logical, intent(in), optional :: no_face_overlap, no_edge_overlap
    this%mesh => mesh
    if (present(no_face_overlap)) this%no_face_overlap = no_face_overlap
    if (present(no_edge_overlap)) this%no_edge_overlap = no_edge_overlap
    call this%builder%init(mesh, bndry_only=.true., omit_offp=.false., &
        no_overlap=this%no_face_overlap)
  end subroutine

  !! Add one boundary-face group containing the union of the face sets in
  !! SETIDS. Face-set lookup, boundary validation, and optional face-overlap
  !! errors are returned through STAT and ERRMSG.
  subroutine add_face_group(this, setids, stat, errmsg)
    class(bndry_edge_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:) ! NB: face set IDs
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
  end subroutine

  !! Return the direct face-group representation maintained by the contained
  !! face builder. The faces in group N are
  !! FACE(XGROUP(N):XGROUP(N+1)-1).
  subroutine get_face_groups(this, ngroup, xgroup, face)
    class(bndry_edge_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), face(:)
    call this%builder%get_face_groups(ngroup, xgroup, face)
  end subroutine

  !! Return the direct edge-group representation. The edges in group N are
  !! EDGE(XGROUP(N):XGROUP(N+1)-1). Edges are ordered and unique within a group,
  !! but may occur in multiple groups. Optional listed on-process edges are
  !! omitted. STAT is 1 if prohibited edge overlap is found, and 0 otherwise.
  subroutine get_edge_groups(this, ngroup, xgroup, edge, stat, omit_edge_list)

    class(bndry_edge_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), edge(:)
    integer, intent(out) :: stat
    integer, intent(in), optional :: omit_edge_list(:)

    type(array_box), allocatable :: glist(:)
    integer :: i, n

    call this%build_edge_groups(glist, stat, omit_edge_list)

    ngroup = size(glist)
    allocate(xgroup(ngroup+1))
    xgroup(1) = 1
    do i = 1, ngroup
      xgroup(i+1) = xgroup(i) + size(glist(i)%array)
    end do

    allocate(edge(xgroup(ngroup+1)-1))
    do n = 1, ngroup
      edge(xgroup(n):xgroup(n+1)-1) = glist(n)%array
    end do

  end subroutine get_edge_groups

  !! Return the mapped edge-group representation. EDGE is ordered and unique
  !! across all groups, and the edges in group N are
  !! EDGE(MAP(XGROUP(N):XGROUP(N+1)-1)). MAP may contain repeated entries.
  !! Optional omission and STAT have the same semantics as GET_EDGE_GROUPS.
  subroutine get_unique_edge_groups(this, ngroup, xgroup, map, edge, stat, omit_edge_list)

    class(bndry_edge_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), map(:), edge(:)
    integer, intent(out) :: stat
    integer, intent(in), optional :: omit_edge_list(:)

    type(array_box), allocatable :: glist(:)
    integer, allocatable :: slot(:)
    integer :: entity, i, j, n

    call this%build_edge_groups(glist, stat, omit_edge_list)

    ngroup = size(glist)
    allocate(xgroup(ngroup+1))
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + size(glist(n)%array)
    end do

    allocate(slot(this%mesh%nedge), source=0)
    do n = 1, ngroup
      do j = 1, size(glist(n)%array)
        entity = glist(n)%array(j)
        ASSERT(entity >= 1 .and. entity <= size(slot))
        slot(entity) = -1
      end do
    end do

    allocate(edge(count(slot /= 0)))
    i = 0
    do entity = 1, size(slot)
      if (slot(entity) /= 0) then
        i = i + 1
        edge(i) = entity
        slot(entity) = i
      end if
    end do

    allocate(map(xgroup(ngroup+1)-1))
    i = 0
    do n = 1, ngroup
      do j = 1, size(glist(n)%array)
        i = i + 1
        map(i) = slot(glist(n)%array(j))
        ASSERT(edge(map(i)) == glist(n)%array(j))
      end do
    end do

    ASSERT(size(map) == xgroup(ngroup+1)-1)
    ASSERT(all(map >= 1 .and. map <= size(edge)))
    if (size(edge) > 1) then
      ASSERT(all(edge(2:) > edge(:size(edge)-1)))
    end if

  end subroutine get_unique_edge_groups

  !! Derive one ordered edge list from each stored face group, synchronize edge
  !! membership across processes, apply optional on-process omissions, and
  !! report prohibited overlap between the resulting groups.
  subroutine build_edge_groups(this, glist, stat, omit_edge_list)

    use parallel_communication, only: global_any

    class(bndry_edge_group_builder), intent(in) :: this
    type(array_box), allocatable, intent(out) :: glist(:)
    integer, intent(out) :: stat
    integer, intent(in), optional :: omit_edge_list(:)

    integer :: i, j, n, ngroup, gsize
    integer, allocatable :: xgroup(:), face(:)
    logical, allocatable :: gmask(:), mask(:)

    stat = 0

    call this%builder%get_face_groups(ngroup, xgroup, face)

    allocate(gmask(this%mesh%nedge), glist(ngroup))
    do i = 1, ngroup
      !! Tag all edges that belong to a face in the group
      gmask = .false.
      associate (gface => face(xgroup(i):xgroup(i+1)-1))
        do j = 1, size(gface)
          associate (fedge => this%mesh%fedge(:,gface(j)))
            gmask(fedge) = .true.
          end associate
        end do
      end associate
      call this%mesh%edge_imap%scatter_offp_or(gmask)
      if (present(omit_edge_list)) then
        do j = 1, size(omit_edge_list) !NB: want to be tolerant of repeated edge indices
          gmask(omit_edge_list(j)) = .false.
        end do
      end if
      !NB: next call overwrites any off-process edges in omit_edge_list
      call this%mesh%edge_imap%gather_offp(gmask)

      !! Generate the list of edge indices
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
      if (this%no_edge_overlap) then
        if (i == 1) then
          if (i < ngroup) mask = gmask ! save for next pass
        else
          if (global_any(mask(glist(i)%array))) stat = 1
          mask(glist(i)%array) = .true.
        end if
      end if
    end do

  end subroutine build_edge_groups

end module bndry_edge_group_builder_type
