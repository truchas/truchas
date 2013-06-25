!!
!! DATA_LAYOUT_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!! Frequently the degrees of freedom of a system will consist of a collection
!! of separately identifiable variables (e.g., cell temperatures, cell
!! enthalpies, face temperatures, boundary face radiosities, etc.) that need
!! to be packed into a single contiguous array of unknowns.  This module
!! provides a simple means of managing the packed layout of such data within
!! a contiguous array and providing access to the data segments of the array.
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines the derived data type DATA_LAYOUT (private components)
!! and the following methods that operate on instances of the type passed as
!! the THIS argument.  An instance of this type is defined by one or more calls
!! to ALLOC_SEGMENT, culminating in a call to ALLOC_COMPLETE.
!!
!! No defined assignment is provided; do not use instances of this type in an
!! assignment statement unless you really know what the default assignment is
!! doing.
!!
!!  ALLOC_SEGMENT(THIS, SIZE) allocates a contiguous segment of length SIZE
!!    in the layout and returns the integer handle used to access this segment.
!!
!!  ALLOC_COMPLETE(THIS) finalizes the layout after all the desired calls to
!!    ALLOC_SEGMENT have been made.  Once called, no further calls to
!!    ALLOC_SEGMENT are permitted.
!!
!!  CALL DELETE_LAYOUT (THIS) deallocates all allocated storage associated
!!    with the layout, returning the object to its default initialization state.
!!
!! Once the layout has been defined, the following methods are available.
!!
!!  LAYOUT_SIZE(THIS) returns the total size of the layout.  The ARRAY argument
!!    to SEGMENT_VIEW, SEGMENT_COPY and DATA_PTR are expected to have this size.
!!
!!  SEGMENT_SIZE(THIS, SEGID) returns the size of the segment with handle SEGID;
!!    this is the size specified in the corresponding call to ALLOC_SEGMENT.
!!
!!  LAYOUT_INDEX(THIS, SEGID, INDEX) returns the index into the layout that
!!    corresponds to the index INDEX of the segment SEGID of the layout.  If
!!    the pointer VIEW is defined by a call to SEGMENT_VIEW then VIEW(INDEX)
!!    and ARRAY(LAYOUT_INDEX(THIS,SEGID,INDEX)) reference the same storage
!!    location.
!!
!!  CALL GET_SEGMENT_VIEW (THIS, ARRAY, SEGID, VIEW) associates the pointer
!!    VIEW with the segment of the rank-1 array ARRAY specified by the segment
!!    handle SEGID of the layout.  ARRAY is either of default real or double
!!    precision type, and VIEW is a pointer to a rank-1 array of the same type.
!!    The size of ARRAY is should equal LAYOUT_SIZE(THIS).
!!    Note very carefully that this call aliases VIEW to (a portion of) ARRAY,
!!    and ARRAY must either be a pointer itself or have the TARGET attribute.
!!
!!  CALL GET_SEGMENT_COPY (THIS, ARRAY, SEGID, COPY) copies the values from the
!!    segment of ARRAY specified by the segment handle SEGID into the array
!!    COPY.  ARRAY is a rank-1 array of default real or double precision type,
!!    and COPY an array of the same type and rank with sufficient size to
!!    hold the copied data.  The size of ARRAY is should equal LAYOUT_SIZE(THIS).
!!
!!  SEGMENT_PTR(THIS, ARRAY, SEGID) returns a pointer to the segment of ARRAY
!!    specified by the segment handle SEGID of the layout.
!!    This is a function variant of GET_SEGMENT_VIEW, and the caution noted
!!    there applies here too.  In addition, note that the function reference
!!    may appear directly in an expression or as the actual argument for an
!!    intent-in dummy argument, but in order to pass the target of the result
!!    to an intent-out or intent-inout dummy argument, a pointer must be
!!    associated with the function result and that pointer passed instead.
!!
!!  CALL SET_SEGMENT (THIS, SOURCE, ARRAY, SEGID) copies values from SOURCE
!!    into the segment of ARRAY specified by the segment handle SEGID.  ARRAY
!!    is either of default real or double precision type, and SOURCE is an
!!    array of the same type and rank.  The size of SOURCE must be at least
!!    SEGMENT_SIZE(THIS, SEGID).  The size of ARRAY is should equal
!!    LAYOUT_SIZE(THIS).
!!

#include "f90_assert.fpp"

module data_layout_type

  implicit none
  private

  public :: data_layout, alloc_segment, alloc_complete, delete_layout
  public :: layout_size, segment_size, layout_index
  public :: get_segment_view, get_segment_copy, set_segment
  public :: segment_ptr  !TODO: should this be deprecated?

  type :: seg_desc
    integer :: lb = 0, ub = 0, size = 0
  end type

  type :: seg_node
    type(seg_desc) :: seg
    type(seg_node), pointer :: next => null()
  end type

  type :: data_layout
    private
    integer :: nseg = 0
    integer :: size = 0
    type(seg_desc), pointer :: seg(:) => null()
    type(seg_node), pointer :: list => null()
  end type data_layout

  interface segment_ptr
    module procedure data_ptr_sp, data_ptr_dp
  end interface

  interface get_segment_view
    module procedure get_segment_view_sp, get_segment_view_dp
  end interface

  interface get_segment_copy
    module procedure get_segment_copy_sp, get_segment_copy_dp
  end interface

  interface set_segment
    module procedure set_segment_sp, set_segment_dp
  end interface

contains

  integer function alloc_segment (this, size) result (segid)

    type(data_layout), intent(inout) :: this
    integer, intent(in) :: size

    type(seg_node), pointer :: new

    ASSERT(size >= 0)
    ASSERT(.not.associated(this%seg))

    this%nseg  = this%nseg + 1
    segid      = this%nseg

    !! Create a new segment descriptor.
    allocate(new)
    new%seg%size = size
    new%seg%lb = this%size + 1
    new%seg%ub = this%size + size
    this%size  = new%seg%ub

    !! Prepend it to the list.
    new%next  => this%list
    this%list => new

  end function alloc_segment


  subroutine alloc_complete (this)

    type(data_layout), intent(inout) :: this

    integer :: id
    type(seg_node), pointer :: next

    ASSERT(.not.associated(this%seg))

    !! Convert the linked list into a static array for quick lookup.
    !! The segment IDs were assigned in sequence and we pop things
    !! off the list in LIFO fashion.  (Perhaps it would be safer to
    !! have actually stored the ID in the seg_desc rather than just
    !! 'know' what it is implicitly.)
    allocate(this%seg(this%nseg))
    do id = this%nseg, 1, -1
      ASSERT(associated(this%list))
      this%seg(id) = this%list%seg
      next => this%list%next
      deallocate(this%list)
      this%list => next
    end do
    ASSERT(.not.associated(this%list))

  end subroutine alloc_complete


  elemental subroutine delete_layout (this)
    type(data_layout), intent(inout) :: this
    type(seg_node), pointer :: next
    if (associated(this%seg)) deallocate(this%seg)
    do while (associated(this%list))
      next => this%list%next
      deallocate(this%list)
      this%list => next
    end do
  end subroutine delete_layout

  pure integer function layout_size (this)
    type(data_layout), intent(in) :: this
    layout_size = this%size
  end function layout_size

  integer function segment_size (this, segid)
    type(data_layout), intent(in) :: this
    integer, intent(in) :: segid
    ASSERT(segid > 0 .and. segid <= this%nseg)
    segment_size = this%seg(segid)%size
  end function segment_size

  integer function layout_index (this, segid, index)
    type(data_layout), intent(in) :: this
    integer, intent(in) :: segid, index
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(index > 0 .and. index <= this%seg(segid)%size)
    layout_index = index + this%seg(segid)%lb - 1
  end function layout_index

  subroutine get_segment_view_sp (this, array, segid, view)
    type(data_layout), intent(in) :: this
    real, intent(in), target :: array(:)
    integer, intent(in) :: segid
    real, pointer :: view(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    view => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine get_segment_view_sp

  subroutine get_segment_view_dp (this, array, segid, view)
    type(data_layout), intent(in) :: this
    double precision, intent(in), target :: array(:)
    integer, intent(in) :: segid
    double precision, pointer :: view(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    view => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine get_segment_view_dp

  subroutine get_segment_copy_sp (this, array, segid, copy)
    type(data_layout), intent(in) :: this
    real, intent(in) :: array(:)
    integer, intent(in) :: segid
    real, intent(inout) :: copy(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    ASSERT(size(copy) >= this%seg(segid)%size)
    copy(:this%seg(segid)%size) = array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine get_segment_copy_sp

  subroutine get_segment_copy_dp (this, array, segid, copy)
    type(data_layout), intent(in) :: this
    double precision, intent(in) :: array(:)
    integer, intent(in) :: segid
    double precision, intent(inout) :: copy(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    ASSERT(size(copy) >= this%seg(segid)%size)
    copy(:this%seg(segid)%size) = array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine get_segment_copy_dp

  function data_ptr_sp (this, array, segid) result (ptr)
    type(data_layout), intent(in) :: this
    real, target, intent(in) :: array(:)
    integer, intent(in) :: segid
    real, pointer :: ptr(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    ptr => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end function data_ptr_sp

  function data_ptr_dp (this, array, segid) result (ptr)
    type(data_layout), intent(in) :: this
    double precision, target, intent(in) :: array(:)
    integer, intent(in) :: segid
    double precision, pointer :: ptr(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    ptr => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end function data_ptr_dp

  subroutine set_segment_sp (this, source, array, segid)
    type(data_layout), intent(in) :: this
    real, intent(in) :: source(:)
    real, intent(inout) :: array(:)
    integer, intent(in) :: segid
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    ASSERT(size(source) >= this%seg(segid)%size)
    array(this%seg(segid)%lb:this%seg(segid)%ub) = source(:this%seg(segid)%size)
  end subroutine set_segment_sp

  subroutine set_segment_dp (this, source, array, segid)
    type(data_layout), intent(in) :: this
    double precision, intent(in) :: source(:)
    double precision, intent(inout) :: array(:)
    integer, intent(in) :: segid
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(associated(this%seg))
    ASSERT(size(source) >= this%seg(segid)%size)
    array(this%seg(segid)%lb:this%seg(segid)%ub) = source(:this%seg(segid)%size)
  end subroutine set_segment_dp

end module data_layout_type
