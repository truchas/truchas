!!
!! BNDRY_EDGE_GROUP_BUILDER_TYPE
!!
!! This module defines an auxiliary object that constructs a grouped list of
!! mesh edges specified incrementally using face set IDs. Its principal use is
!! in the instantiation of boundary condition objects.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! January 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bndry_edge_group_builder_type

  use simpl_mesh_type
  use bndry_face_group_builder_type
  implicit none
  private

  type, public :: bndry_edge_group_builder
    private
    type(simpl_mesh), pointer :: mesh => null() ! reference only
    type(bndry_face_group_builder) :: builder
    logical :: no_overlap = .true.
  contains
    procedure :: init
    procedure :: add_face_group
    procedure :: get_face_groups
    procedure :: get_edge_groups
  end type

contains

  subroutine init(this, mesh, no_overlap)
    class(bndry_edge_group_builder), intent(out) :: this
    type(simpl_mesh), target :: mesh
    logical, intent(in), optional :: no_overlap
    this%mesh => mesh
    if (present(no_overlap)) this%no_overlap = no_overlap
    call this%builder%init(mesh, bndry_only=.true., omit_offp=.false., no_overlap=no_overlap)
  end subroutine

  subroutine add_face_group(this, setids, stat, errmsg)
    class(bndry_edge_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:) ! NB: face set IDs
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
  end subroutine

  subroutine get_face_groups(this, ngroup, xgroup, index)
    class(bndry_edge_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:)
    call this%builder%get_face_groups(ngroup, xgroup, index)
  end subroutine

  subroutine get_edge_groups(this, ngroup, xgroup, index, stat, omit_edge_list)

    use parallel_communication, only: global_any

    class(bndry_edge_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:)
    integer, intent(out) :: stat
    integer, intent(in), optional :: omit_edge_list(:)

    integer :: i, j, n, gsize
    logical, allocatable :: gmask(:), mask(:)

    type :: array_box
      integer, allocatable :: array(:)
    end type
    type(array_box), allocatable :: glist(:)

    stat = 0

    call this%builder%get_face_groups(ngroup, xgroup, index)

    allocate(gmask(this%mesh%nedge), glist(ngroup))
    do i = 1, ngroup
      !! Tag all edges that belong to a face in the group
      gmask = .false.
      associate (gface => index(xgroup(i):xgroup(i+1)-1))
        do j = 1, size(gface)
          associate (fedge => abs(this%mesh%fedge(:,gface(j))))
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
      call this%mesh%edge_imap%gather_offp(gmask)

      !! Generate the list of edge indices
      gsize = count(gmask(:this%mesh%nedge_onp))
      allocate(glist(i)%array(gsize))
      n = 0
      do j = 1, size(gmask)
        if (n >= gsize) exit ! no more to be found
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
  end subroutine get_edge_groups

end module bndry_edge_group_builder_type
