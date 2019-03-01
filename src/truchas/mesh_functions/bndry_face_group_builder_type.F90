!!
!! BNDRY_FACE_GROUP_BUILDER_TYPE
!!
!! This module defines an auxiliary object that constructs a grouped list of
!! boundary faces specified incrementally using face set IDs.  Its principal
!! use is in instantiating boundary condition objects.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The BNDRY_FACE_GROUP_BUILDER type has the following type bound procedures.
!!
!!  INIT(MESH) initializes the object. MESH is of type UNSTR_MESH.
!!
!!  ADD_FACE_GROUP(SETIDS, STAT, ERRMSG) defines a group of boundary faces
!!    comprising those faces that belong to a face set with ID in the integer
!!    array SETIDS. It is an error if any of the specified faces belong to
!!    a previously defined group, or are not a boundary face.  STAT returns 1
!!    in the former case and 2 in the latter, otherwise STAT returns 0. An
!!    message is returned in the deferred-length allocatable character ERRMSG
!!    if an error occurs.
!!
!!  GET_FACE_GROUPS(NGROUP, XGROUP, INDEX) returns the face groups defined by
!!    the previous calls to ADD_FACE_GROUP.  The array INDEX returns the list
!!    of face indices. INDEX(XGROUP(n):XGROUP(n)-1) are the face indices that
!!    belong to group n. The groups are sequentially numbered from 1 in the
!!    order of the calls to ADD_FACE_GROUP.  NGROUP is the number of groups.
!!    Both XGROUP and INDEX are allocatable arrays and are allocated by this
!!    procedure.
!!

module bndry_face_group_builder_type

  use base_mesh_class
  implicit none
  private

  type, public :: bndry_face_group_builder
    private
    class(base_mesh), pointer :: mesh => null() ! reference only -- not owned
    integer :: ngroup = 0
    integer, allocatable :: tag(:)
    logical, allocatable :: mask(:) ! work space for add_group
  contains
    procedure :: init
    procedure :: add_face_group
    procedure :: get_face_groups
  end type bndry_face_group_builder

contains

  subroutine init(this, mesh)
    class(bndry_face_group_builder), intent(out) :: this
    class(base_mesh), target :: mesh
    this%mesh => mesh
    this%ngroup = 0
    allocate(this%tag(mesh%nface), this%mask(mesh%nface))
    this%tag = 0
  end subroutine init

  subroutine add_face_group(this, setids, stat, errmsg)

    use bitfield_type
    use parallel_communication, only: global_count
    use string_utilities, only: i_to_c

    class(bndry_face_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(bitfield) :: bitmask

    !! Identify the faces specified by SETIDS.
    call this%mesh%get_face_set_bitmask(setids, bitmask, stat, errmsg)
    if (stat /= 0) return
    this%mask = (popcnt(iand(bitmask, this%mesh%face_set_mask)) /= 0)

    !! Check that these faces don't overlap those from preceding calls.
    n = global_count(this%mask .and. this%tag /= 0)
    if (n /= 0) then
      stat = 1
      errmsg = i_to_c(n) // ' faces belong to an existing group'
      return
    end if

    !! Verify that these faces are boundary faces.
    n = global_count(this%mask .and. .not.btest(this%mesh%face_set_mask,0))
    if (n /= 0) then
      stat = 2
      errmsg = i_to_c(n) // ' faces not on boundary'
      return
    end if

    !! Set the tag array.
    this%ngroup = 1 + this%ngroup
    where (this%mask) this%tag = this%ngroup

  end subroutine add_face_group

  subroutine get_face_groups(this, ngroup, xgroup, index)

    class(bndry_face_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:)

    integer :: n, j

    ngroup = this%ngroup
    n = count(this%tag > 0)
    allocate(index(n), xgroup(ngroup+1))

    !! Prepare XGROUP: indices in group N will be INDEX(XGROUP(N):XGROUP(N+1)-1).
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + count(this%tag == n)
    end do

    !! Fill the INDEX array; XGROUP(N) stores the next free location for group N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      index(xgroup(n)) = j
      xgroup(n) = 1 + xgroup(n)
    end do

    !! Restore XGROUP; XGROUP(N) is now the start of group N+1 instead of N.
    do n = ngroup, 1, -1
      xgroup(n+1) = xgroup(n)
    end do
    xgroup(1) = 1

  end subroutine get_face_groups

end module bndry_face_group_builder_type
