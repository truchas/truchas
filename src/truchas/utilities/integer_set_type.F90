!!
!! INTEGER_SET_TYPE
!!
!! This module defines a set container that stores unique integer values.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for Fortran 2008, March 2015.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type INTEGER_SET implements a set container that stores unique
!!  integer values.  It has the following elemental type bound procedures.
!!
!!    IS_EMPTY() returns true if the set contains no values.
!!
!!    SIZE() returns the number of values in the set.
!!
!!    ADD (VALUE) adds the integer VALUE to the set.  If the set already
!!      stores the value it is not added again.
!!
!!    ADD (SET) adds the elements of the INTEGER_SET argument SET to the
!!      object.  SET is unmodified by the call.
!!
!!    REMOVE (VALUE) removes the integer VALUE from the set if it contains it.
!!
!!  The following defined assignment is provided:
!!
!!    TYPE(INTEGER_SET) :: SET
!!    INTEGER, ALLOCATABLE :: ARRAY(:)
!!    ARRAY = SET
!!
!!  The sorted values stored in SET are copied, in increasing order, to the
!!  elements of ARRAY. The assignment operation ensures that the size of ARRAY
!!  equals SET%SIZE() by automatically allocating or reallocating the array
!!  when necessary.
!!
!!  Scalar and array objects of this type are properly finalized when they
!!  are deallocated or otherwise cease to exist.
!!
!! IMPLEMENTATION NOTES
!!
!!  The values are stored in increasing order using a simple linked-list
!!  structure and are added using a linear search.  This is most appropriate
!!  for relatively small sets only.
!!

module integer_set_type

  implicit none
  private

  type, public :: integer_set
    private
    type(set_item), pointer :: head => null()
    integer :: n = 0  ! number of elements in the set
  contains
    procedure :: is_empty => set_is_empty
    procedure :: size => set_size
    procedure, private :: set_add
    procedure, private :: set_add_set
    procedure, private :: set_add_array
    generic   :: add => set_add, set_add_set, set_add_array
    procedure :: remove => set_remove
    procedure :: copy_to_array
    procedure, private, pass(rhs) :: set_to_array
    generic :: assignment(=) => set_to_array
    final :: integer_set_delete
    procedure :: begin => set_begin
  end type integer_set

  type :: set_item
    integer :: value
    type(set_item), pointer :: next => null()
  end type set_item

  type, public :: integer_set_iterator
    class(set_item), pointer, private :: item => null()
  contains
    procedure :: next => iter_next
    procedure :: at_end => iter_at_end
    procedure :: value => iter_value
  end type integer_set_iterator

  !! Defined INTEGER_SET_ITERATOR structure constructor
  interface integer_set_interator
    procedure set_begin
  end interface

contains

  !! Final subroutine for INTEGER_SET objects
  elemental subroutine integer_set_delete (this)
    type(integer_set), intent(inout) :: this
    type(set_item), pointer :: p, item
    p => this%head
    do while (associated(p))
      item => p
      p => p%next
      deallocate(item)
    end do
  end subroutine integer_set_delete

  !! Returns true if the set is empty.
  elemental logical function set_is_empty (this)
    class(integer_set), intent(in) :: this
    set_is_empty = .not.associated(this%head)
  end function set_is_empty

  !! Returns the number of elements in the set.
  elemental integer function set_size (this)
    class(integer_set), intent(in) :: this
    set_size = this%n
  end function set_size

  !! Adds a value to the set; no duplicates.
  elemental subroutine set_add (this, value)
    class(integer_set), intent(inout) :: this
    integer, intent(in) :: value
    type(set_item), pointer :: item, p, t
    p => this%head; t => null()
    do while (associated(p))
      if (value == p%value) return
      if (value < p%value) exit
      t => p
      p => p%next
    end do
    allocate(item)
    item%value = value
    item%next => p
    if (associated(t)) then
      t%next => item
    else
      this%head => item
    end if
    this%n = this%n + 1
  end subroutine set_add

  !! Adds the values from the given set.
  subroutine set_add_set (this, set)
    class(integer_set), intent(inout) :: this
    class(integer_set), intent(in) :: set
    type(set_item), pointer :: item, p, t, q
    p => this%head; t => null(); q => set%head
    nextq: do while (associated(q))
      do while (associated(p))
        if (q%value == p%value) then
          q => q%next
          t => p
          p => p%next
          cycle nextq
        end if
        if (q%value < p%value) exit
        t => p
        p => p%next
      end do
      allocate(item)
      item%value = q%value
      item%next => p
      if (associated(t)) then
        t%next => item
      else
        this%head => item
      end if
      this%n = this%n + 1
      t => item
      q => q%next
    end do nextq
  end subroutine set_add_set

  !! Adds values from an array to the set; no duplicates.
  subroutine set_add_array (this, values)
    class(integer_set), intent(inout) :: this
    integer, intent(in) :: values(:)
    integer :: j
    do j = 1, size(values)
      call set_add (this, values(j))
    end do
  end subroutine set_add_array

  !! Removes a value from the set.
  elemental subroutine set_remove (this, value)
    class(integer_set), intent(inout) :: this
    integer, intent(in) :: value
    type(set_item), pointer :: p, t
    p => this%head; t => null()
    do while (associated(p))
      if (value == p%value) then
        if (associated(t)) then
          t%next => p%next
        else
          this%head => p%next
        end if
        deallocate(p)
        this%n = this%n - 1
        return
      end if
      if (value < p%value) return
      t => p
      p => p%next
    end do
  end subroutine set_remove

  !! Defined assignment of the set to an allocatable rank-1 array.
  !! The array is allocated/reallocated to the correct size, and the
  !! sorted values in the set are written to the array.  If the set
  !! is empty the array is allocated with 0 size.
  subroutine set_to_array (lhs, rhs)
    integer, allocatable, intent(inout) :: lhs(:)
    class(integer_set), intent(in) :: rhs
    type(set_item), pointer :: p
    integer :: n
    if (allocated(lhs)) then
      if (size(lhs) /= rhs%n) deallocate(lhs)
    end if
    if (.not.allocated(lhs)) allocate(lhs(rhs%n))
    n = 0
    p => rhs%head
    do while (associated(p))
      n = n + 1
      lhs(n) = p%value
      p => p%next
    end do
  end subroutine set_to_array

  !! Copy the sorted values in the set to the user provided array.  Its size
  !! must be sufficiently large to store the values.  If larger, the values
  !! of the unused elements are left unchanged.
  subroutine copy_to_array (this, array)
    integer, intent(inout) :: array(:)
    class(integer_set), intent(in) :: this
    type(set_item), pointer :: p
    integer :: n
    n = 0
    p => this%head
    do while (associated(p))
      n = n + 1
      array(n) = p%value
      p => p%next
    end do
  end subroutine copy_to_array

  !! Returns an INTEGER_SET_ITERATOR positioned to the first element of the set.
  function set_begin (this) result (iter)
    class(integer_set), intent(in) :: this
    type(integer_set_iterator) :: iter
    iter%item => this%head
  end function set_begin

!!!! INTEGER_SET_TYPE_ITERATOR TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Advances the iterator to the next element in the set.
  pure subroutine iter_next (this)
    class(integer_set_iterator), intent(inout) :: this
    if(associated(this%item)) this%item => this%item%next
  end subroutine iter_next

  !! Returns true if the iterator has reached the end; that is, it has
  !! gone past the last element of the set.
  pure logical function iter_at_end (this)
    class(integer_set_iterator), intent(in) :: this
    iter_at_end = .not.associated(this%item)
  end function iter_at_end

  !! Returns the value of the current element.
  pure integer function iter_value (this)
    class(integer_set_iterator), intent(in) :: this
    iter_value = this%item%value
  end function iter_value

end module integer_set_type
