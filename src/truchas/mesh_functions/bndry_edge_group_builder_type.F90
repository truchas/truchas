!!
!! BNDRY_EDGE_GROUP_BUILDER_TYPE
!!
!! This module defines an auxiliary object that constructs a grouped list of
!! mesh edges specified incrementally using face set IDs. Its principal use is
!! in the instantiation of boundary condition objects.
!!
!! Face overlap and edge overlap are independent policies. Distinct face groups
!! may legitimately share edges even when faces may belong to only one group.
!! GET_EDGE_GROUPS returns the direct grouped representation, where XGROUP
!! partitions EDGE and EDGE may contain duplicates across groups.
!! GET_UNIQUE_EDGE_GROUPS instead returns an ordered unique EDGE array and a
!! grouped MAP, where EDGE(MAP(XGROUP(N):XGROUP(N+1)-1)) gives group N.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! January 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    type(simpl_mesh), pointer :: mesh => null() ! reference only
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

  subroutine add_face_group(this, setids, stat, errmsg)
    class(bndry_edge_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:) ! NB: face set IDs
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
  end subroutine

  subroutine get_face_groups(this, ngroup, xgroup, face)
    class(bndry_edge_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), face(:)
    call this%builder%get_face_groups(ngroup, xgroup, face)
  end subroutine

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
