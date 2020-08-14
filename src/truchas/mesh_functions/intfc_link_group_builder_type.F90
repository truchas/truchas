!!
!! INTFC_LINK_GROUP_BUILDER_TYPE
!!
!! This module defines an auxiliary object that constructs a grouped list
!! of interface link faces specified incrementally using link set IDs.
!! Its principal use is in instantiating interface condition objects.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2018. Adapted from bndry_face_group_builder.F90
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The INTFC_LINK_GROUP_BUILDER type has the following type bound procedures.
!!
!!  INIT(MESH) initializes the object. MESH is of UNSTR_BASE_MESH class.
!!
!!  ADD_LINK_GROUP(SETIDS, STAT, ERRMSG) defines a group of interface links
!!    comprising those links that belong to a link set with ID in the integer
!!    array SETIDS. It is an error if any of the specified links belong to
!!    a previously defined group, or are not a boundary face.  STAT returns 1
!!    in the former case and 2 in the latter, otherwise STAT returns 0. An
!!    message is returned in the deferred-length allocatable character ERRMSG
!!    if an error occurs.
!!
!!  GET_LINK_GROUPS(NGROUP, XGROUP, INDEX) returns the link groups defined by
!!    the previous calls to ADD_LINK_GROUP.  The array INDEX returns the list
!!    of pairs of link face indices. INDEX(:,XGROUP(n):XGROUP(n)-1) are those
!!    that belong to group n. The groups are sequentially numbered from 1 in
!!    the order of the calls to ADD_FACE_GROUP. NGROUP is the number of groups.
!!    Both XGROUP and INDEX are allocatable arrays and are allocated by this
!!    procedure.
!!

module intfc_link_group_builder_type

  use unstr_base_mesh_class
  implicit none
  private

  type, public :: intfc_link_group_builder
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only -- not owned
    integer :: ngroup = 0
    integer, allocatable :: tag(:)
    logical, allocatable :: mask(:) ! work space for add_group
  contains
    procedure :: init
    procedure :: add_link_group
    procedure :: get_link_groups
  end type intfc_link_group_builder

contains

  subroutine init(this, mesh)
    class(intfc_link_group_builder), intent(out) :: this
    class(unstr_base_mesh), target :: mesh
    this%mesh => mesh
    this%ngroup = 0
    allocate(this%tag(mesh%nlink), this%mask(mesh%nlink))
    this%tag = 0
  end subroutine init

  subroutine add_link_group(this, setids, stat, errmsg)

    use bitfield_type
    use parallel_communication, only: global_count
    use string_utilities, only: i_to_c

    class(intfc_link_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(bitfield) :: bitmask

    !! Identify the links specified by SETIDS.
    call this%mesh%get_link_set_bitmask(setids, bitmask, stat, errmsg)
    if (stat /= 0) return
    this%mask = (popcnt(iand(bitmask, this%mesh%link_set_mask)) /= 0)

    !! Check that these links don't overlap those from preceding calls.
    n = global_count(this%mask .and. this%tag /= 0)
    if (n /= 0) then
      stat = 1
      errmsg = i_to_c(n) // ' links belong to an existing group'
      return
    end if

    !! Set the tag array.
    this%ngroup = 1 + this%ngroup
    where (this%mask) this%tag = this%ngroup

  end subroutine add_link_group

  subroutine get_link_groups(this, ngroup, xgroup, index)

    class(intfc_link_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:,:)

    integer :: n, j

    ngroup = this%ngroup
    n = count(this%tag > 0)
    allocate(index(2,n), xgroup(ngroup+1))

    !! Prepare XGROUP: indices in group N will be INDEX(XGROUP(N):XGROUP(N+1)-1).
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + count(this%tag == n)
    end do

    !! Fill the INDEX array; XGROUP(N) stores the next free location for group N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      index(:,xgroup(n)) = this%mesh%lface(:,j)
      xgroup(n) = 1 + xgroup(n)
    end do

    !! Restore XGROUP; XGROUP(N) is now the start of group N+1 instead of N.
    do n = ngroup, 1, -1
      xgroup(n+1) = xgroup(n)
    end do
    xgroup(1) = 1

  end subroutine get_link_groups

end module intfc_link_group_builder_type
