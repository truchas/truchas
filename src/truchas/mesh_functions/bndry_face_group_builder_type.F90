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
!!  INIT(MESH [,BNDRY_ONLY]) initializes the object. MESH is of type BASE_MESH.
!!    Faces will be required to be boundary faces unless the option BNDRY_ONLY
!!    is present with value .false.
!!
!!  ADD_FACE_GROUP(SETIDS, STAT, ERRMSG) defines a group of boundary faces
!!    comprising those faces that belong to a face set with ID in the integer
!!    array SETIDS. It is an error if any of the specified faces belong to
!!    a previously defined group, or are not a boundary face and BNDRY_ONLY
!!    was not specified with value .false.  STAT returns 1 in the former case
!!    and 2 in the latter, otherwise STAT returns 0. A message is returned in
!!    the deferred-length allocatable character ERRMSG if an error occurs.
!!
!!  GET_FACE_GROUPS(NGROUP, XGROUP, INDEX [,OMIT_OFFP]) returns the face groups
!!    defined by the previous calls to ADD_FACE_GROUP.  The array INDEX returns
!!    the list of face indices: INDEX(XGROUP(n):XGROUP(n)-1) are the face
!!    indices that belong to group n. The groups are sequentially numbered
!!    from 1 in the order of the calls to ADD_FACE_GROUP.  NGROUP is the number
!!    of groups. Both XGROUP and INDEX are allocatable arrays and are allocated
!!    by this procedure. If the option OMIT_OFFP is specified with value .true.
!!    then the face groups will only include on-process faces; otherwise, all
!!    specified faces are included irrespective of whether they are on or off-
!!    process, which is the default.
!!

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
  end type bndry_face_group_builder

contains

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

  subroutine get_face_groups(this, ngroup, xgroup, index)

    class(bndry_face_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:)

    integer :: n, j

    ngroup = size(this%glist)
    n = 0
    do j = 1, ngroup
      n = n + size(this%glist(j)%array)
    end do
    allocate(index(n), xgroup(ngroup+1))

    !! Face indices in group N are INDEX(XGROUP(N):XGROUP(N+1)-1).
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + size(this%glist(n)%array)
      index(xgroup(n):xgroup(n+1)-1) = this%glist(n)%array
    end do

  end subroutine get_face_groups

end module bndry_face_group_builder_type
