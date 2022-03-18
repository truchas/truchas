!!
!! UNSTR_BASE_MESH_CLASS
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module unstr_base_mesh_class

  use kinds, only: r8
  use base_mesh_class
  use parallel_communication
  use index_map_type
  use bitfield_type
  implicit none
  private

  type, abstract, extends(base_mesh), public :: unstr_base_mesh
    !! Mesh interface links.
    integer :: nlink=0, nlink_onP=0
    integer, allocatable :: lface(:,:)      ! pointer due to localize_index_array
    integer, allocatable :: link_set_id(:)  ! user-assigned ID for each link block
    type(bitfield), allocatable :: link_set_mask(:)  ! link block index
    type(index_map) :: link_imap
  contains
    procedure :: get_link_set_bitmask
    procedure :: get_link_set_ids
    procedure(conn_list), deferred :: cell_node_list_view
    procedure(conn_list), deferred :: cell_face_list_view
    procedure(conn_list), deferred :: face_node_list_view
  end type unstr_base_mesh

  abstract interface
    !! Returns connectivity data of the given mesh element
    function conn_list(this, n) result(view)
      import unstr_base_mesh
      class(unstr_base_mesh), intent(in), target :: this
      integer, intent(in) :: n
      integer, pointer, contiguous :: view(:)
    end function
  end interface

contains

  !! Returns a scalar bit mask for use in bit operations with the link_set_mask
  !! array component.  The corresponding bit is set for each link set ID given
  !! in the array SETIDS.  STAT returns a non-zero value if an unknown link set
  !! ID is specified, and the optional allocatable deferred-length character
  !! ERRMSG is assigned an explanatory message if present.

  subroutine get_link_set_bitmask (this, setids, bitmask, stat, errmsg)
    use string_utilities, only: i_to_c
    class(unstr_base_mesh), intent(in) :: this
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

    class(unstr_base_mesh), intent(in) :: this
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

end module unstr_base_mesh_class
