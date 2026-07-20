!!
!! BNDRY_FACE_GROUP_BUILDER_TYPE
!!
!! This module defines a builder that incrementally converts mesh face-set
!! selections into grouped face-index representations. It is primarily used to
!! construct sparse boundary functions and related mesh-data objects.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, December 2017, modified July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The builder retains an unowned reference to the mesh. After initialization,
!! zero or more face groups are added using face-set IDs. Groups are numbered in
!! insertion order and may optionally be restricted to boundary faces, omit
!! off-process faces, or reject overlap with preceding groups.
!!
!! Two final representations are available. The direct representation
!! concatenates the groups and may repeat faces across groups when overlap is
!! allowed. The mapped representation instead returns one ordered, unique face
!! array and a grouped map into it. All result arrays are allocated by the
!! getter procedures.
!!

#include "f90_assert.fpp"

module bndry_face_group_builder_type

  use base_mesh_class
  implicit none
  private

  type :: array_box
    integer, allocatable :: array(:)
  end type

  type, public :: bndry_face_group_builder
    private
    class(base_mesh), pointer :: mesh => null() ! reference only -- not owned
    logical :: bndry_only = .true., omit_offp = .false., no_overlap = .true.
    logical, allocatable :: tag(:)
    logical, allocatable :: mask(:) ! work space for add_group
    type(array_box), allocatable :: glist(:)
  contains
    procedure :: init
    procedure :: add_face_group
    procedure :: get_face_groups
    procedure :: get_unique_face_groups
  end type bndry_face_group_builder

contains

  !! Initialize an empty builder for MESH. Boundary-only validation and overlap
  !! rejection are enabled by default. Off-process faces are included by
  !! default. The optional arguments override those policies.
  subroutine init(this, mesh, bndry_only, omit_offp, no_overlap)
    class(bndry_face_group_builder), intent(out) :: this
    class(base_mesh), target :: mesh
    logical, intent(in), optional :: bndry_only, omit_offp, no_overlap
    this%mesh => mesh
    if (present(bndry_only)) this%bndry_only = bndry_only
    if (present(omit_offp))  this%omit_offp  = omit_offp
    if (present(no_overlap)) this%no_overlap = no_overlap
    allocate(this%mask(merge(mesh%nface_onP, mesh%nface, this%omit_offp)))
    if (this%no_overlap) allocate(this%tag(mesh%nface_onP), source=.false.)
    allocate(this%glist(0))
  end subroutine

  !! Add one group containing the union of the face sets in SETIDS. This is a
  !! collective operation. Errors from face-set lookup are propagated; the
  !! builder assigns STAT=1 for overlap with an existing group and STAT=2 for
  !! non-boundary faces when the corresponding checks are enabled.
  subroutine add_face_group(this, setids, stat, errmsg)

    use bitfield_type
    use parallel_communication, only: global_sum
    use string_utilities, only: i_to_c

    class(bndry_face_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n, ngface, nbndry, novrlp
    type(bitfield) :: bitmask

    !! Identify the faces specified by SETIDS.
    call this%mesh%get_face_set_bitmask(setids, bitmask, stat, errmsg)
    if (stat /= 0) return
    ngface = 0  ! count of faces in this group
    nbndry = 0  ! count of non-boundary faces (on-process)
    novrlp = 0  ! count of overlapping faces (on-process)
    do j = 1, size(this%mask)
      this%mask(j) = (popcnt(iand(bitmask, this%mesh%face_set_mask(j))) /= 0)
      if (.not.this%mask(j)) cycle
      if (j <= this%mesh%nface_onP) then
        if (.not.btest(this%mesh%face_set_mask(j),0)) nbndry = nbndry + 1
        if (allocated(this%tag)) then
          if (this%tag(j)) novrlp = novrlp + 1
          this%tag(j) = .true.
        end if
      end if
      ngface = ngface + 1
    end do

    !! If requested, verify these faces do not overlap a previous group.
    if (this%no_overlap) then
      n = global_sum(novrlp)
      if (n /= 0) then
        stat = 1
        errmsg = i_to_c(n) // ' faces belong to an existing group'
        return
      end if
    end if

    !! If requested, verify that these faces are boundary faces.
    if (this%bndry_only) then
      n = global_sum(nbndry)
      if (n /= 0) then
        stat = 2
        errmsg = i_to_c(n) // ' faces not on boundary'
        return
      end if
    end if

    !! Make space to store the list of faces for this group.
    block ! resize glist
      type(array_box), allocatable :: tmp(:)
      call move_alloc(this%glist, tmp)
      allocate(this%glist(size(tmp)+1))
      do j = 1, size(tmp)
        call move_alloc(tmp(j)%array, this%glist(j)%array)
      end do
      deallocate(tmp)
    end block

    !! Store the list of faces for this group.
    allocate(this%glist(size(this%glist))%array(ngface))
    associate (array => this%glist(size(this%glist))%array)
      n = 0
      do j = 1, size(this%mask)
        if (this%mask(j)) then
          n = n + 1
          array(n) = j
        end if
      end do
    end associate

  end subroutine add_face_group

  !! Return the direct group representation. The faces in group N are
  !! FACE(XGROUP(N):XGROUP(N+1)-1). Faces are ordered and unique within a group,
  !! but may occur in multiple groups when overlap is allowed.
  subroutine get_face_groups(this, ngroup, xgroup, face)

    class(bndry_face_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), face(:)

    integer :: n, j

    ngroup = size(this%glist)
    n = 0
    do j = 1, ngroup
      n = n + size(this%glist(j)%array)
    end do
    allocate(face(n), xgroup(ngroup+1))

    !! Face indices in group N are FACE(XGROUP(N):XGROUP(N+1)-1).
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + size(this%glist(n)%array)
      face(xgroup(n):xgroup(n+1)-1) = this%glist(n)%array
    end do

  end subroutine get_face_groups

  !! Return the mapped group representation. FACE is ordered and unique across
  !! all groups, and the faces in group N are
  !! FACE(MAP(XGROUP(N):XGROUP(N+1)-1)). MAP may contain repeated entries when
  !! overlap is allowed.
  subroutine get_unique_face_groups(this, ngroup, xgroup, map, face)

    class(bndry_face_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), map(:), face(:)

    integer, allocatable :: slot(:)
    integer :: entity, i, j, n

    ngroup = size(this%glist)
    allocate(xgroup(ngroup+1))
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + size(this%glist(n)%array)
    end do

    allocate(slot(size(this%mask)), source=0)
    do n = 1, ngroup
      do j = 1, size(this%glist(n)%array)
        entity = this%glist(n)%array(j)
        ASSERT(entity >= 1 .and. entity <= size(slot))
        slot(entity) = -1
      end do
    end do

    allocate(face(count(slot /= 0)))
    i = 0
    do entity = 1, size(slot)
      if (slot(entity) /= 0) then
        i = i + 1
        face(i) = entity
        slot(entity) = i
      end if
    end do

    allocate(map(xgroup(ngroup+1)-1))
    i = 0
    do n = 1, ngroup
      do j = 1, size(this%glist(n)%array)
        i = i + 1
        map(i) = slot(this%glist(n)%array(j))
      end do
    end do

    ASSERT(size(map) == xgroup(ngroup+1)-1)
    ASSERT(all(map >= 1 .and. map <= size(face)))
    if (size(face) > 1) then
      ASSERT(all(face(2:) > face(:size(face)-1)))
    end if

  end subroutine get_unique_face_groups

end module bndry_face_group_builder_type
