!!
!! SETS
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 27 Feb 2007; last revised 4 Apr 2007.
!!
!! This module provides the derived type INTEGER_SET for representing
!! a collection of integers that contains no duplicates, and procedures
!! that operate on instances of this type.
!!
!! PROGRAMMING INTERFACE
!!
!!  The type INTEGER_SET is opaque, having private components.  An instance of
!!  this type in its default initialization state represents the empty set, with
!!  IS_EMPTY returning true.
!!
!!  CALL ADD (THIS, N) adds the default integer N to the set THIS if it does not
!!    already contain N.
!!
!!  CALL CLEAR (THIS) removes all elements from the set THIS, returning it to
!!    its default initialization state (empty set).
!!
!!  IS_EMPTY(THIS) returns the value true if the set THIS contains no elements;
!!    otherwise it returns the value false.
!!
!!  SIZE(THIS) returns the number of elements in the set THIS (its cardinality).
!!
!!  TO_ARRAY(THIS) returns a pointer to a rank-1 default integer array containing
!!    all the elements of the set THIS in increasing order.  The size of the
!!    result equals SIZE(THIS).
!!
!!    N.B.  The only proper use of this function is as the target of a pointer
!!    assignment; if used otherwise (in an expression, e.g.) a memory leak will
!!    result.  The caller is responsible for deallocating the array.
!!
!! IMPLEMENTATION NOTES
!!
!!  1. The set elements are stored in increasing order using a recursive linked-
!!    list data structure.  This should probably be reimplemented using a
!!    self-balanced binary search tree.
!!
!!  2. I would much prefer TO_ARRAY to return an array rather than a pointer to
!!    an array.  But to do this SET_SIZE would need to be a pure procedure and
!!    this isn't possible because its argument is an implicit target of a
!!    (completely innocuous) pointer assignment.  Having TO_ARRAY return an
!!    allocatable array result is an alternative, but this language feature is
!!    not yet widely supported.
!!

module sets

  implicit none
  private

  public :: is_empty, clear, add, to_array, size

  type, public :: integer_set
    private
    type(integer_set_member), pointer :: first => null()
  end type

  type :: integer_set_member
    integer :: value
    type(integer_set) :: rest
  end type

  interface size
    module procedure set_size
  end interface

contains

  logical function is_empty (this)
    type(integer_set), intent(in) :: this
    is_empty = .not.associated(this%first)
  end function is_empty

  recursive subroutine clear (this)
    type(integer_set), intent(inout) :: this
    if (associated(this%first)) then
      call clear (this%first%rest)
      deallocate(this%first)
    end if
  end subroutine clear

  recursive subroutine add (this,  n)
    type(integer_set), intent(inout) :: this
    integer, intent(in) :: n
    type(integer_set) :: rest
    if (.not.associated(this%first)) then
      allocate(this%first)
      this%first%value = n
    else if (n < this%first%value) then
      rest = this
      allocate(this%first)
      this%first%value = n
      this%first%rest = rest
    else if (n > this%first%value) then
      call add (this%first%rest, n)
    end if
  end subroutine add

  function to_array (this) result (array)
    type(integer_set), intent(in) :: this
    integer, pointer :: array(:)
    integer :: n
    type(integer_set) :: l
    allocate(array(set_size(this)))
    n = 0
    l = this
    do while (associated(l%first))
      n = n + 1
      array(n) = l%first%value
      l = l%first%rest
    end do
  end function to_array

  integer function set_size (this)
    type(integer_set), intent(in) :: this
    type(integer_set) :: l
    set_size = 0
    l = this
    do while (associated(l%first))
      set_size = set_size + 1
      l = l%first%rest
    end do
  end function set_size

end module sets
