!!
!! STRING_SET_TYPE
!!
!! This module defines a set container that stores unique character string
!! values. This is very rudimentary with only the capabilities needed for the
!! moment. Only to be used for small sets and not in perfomance sensitive
!! contexts.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module string_set_type

  implicit none
  private

  type, public :: string_set
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: add
    procedure :: has
    procedure :: copy
    final :: string_set_delete
  end type

  type :: list_item
    character(:), allocatable :: string
    type(list_item), pointer :: next => null()
  contains
    final :: list_item_delete
  end type

  type, public :: string_set_iterator
    private
    class(list_item), pointer :: item => null()
  contains
    procedure :: next => iter_next
    procedure :: at_end => iter_at_end
    procedure :: string => iter_string
  end type

  !! User-defined STRING_SET_ITERATOR structure constructor
  interface string_set_iterator
    procedure string_set_begin
  end interface

contains

  subroutine string_set_delete(this)
    type(string_set), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine

  recursive subroutine list_item_delete(this)
    type(list_item), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  subroutine add(this, string)
    class(string_set), intent(inout) :: this
    character(*), intent(in) :: string
    type(list_item), pointer :: item
    item => find_list_item(this, string)
    if (associated(item)) return
    allocate(item)
    item%string = string
    item%next => this%first
    this%first => item
  end subroutine

  logical function has(this, string)
    class(string_set), intent(in) :: this
    character(*), intent(in) :: string
    has = associated(find_list_item(this, string))
  end function

  function find_list_item(this, string) result(item)
    class(string_set), intent(in) :: this
    character(*), intent(in) :: string
    type(list_item), pointer :: item
    item => this%first
    do while (associated(item))
      if (item%string == string) return
      item => item%next
    end do
  end function

  !! Copy the map SRC to DEST.
  subroutine copy(src, dest)
    class(string_set), intent(in) :: src
    type(string_set), intent(inout) :: dest
    type(string_set_iterator) :: iter
    iter = string_set_iterator(src)
    do while (.not.iter%at_end())
      call dest%add(iter%string())
      call iter%next
    end do
  end subroutine

  !!!! STRING_SET_ITERATOR TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Defined STRING_SET_ITERATOR constructor that is positioned
  !! to the beginning element of the specified STRING_SET object.
  function string_set_begin(this) result(iter)
    class(string_set), intent(in) :: this
    type(string_set_iterator) :: iter
    iter%item => this%first
  end function

  !! Advances the iterator to the next element in the set.
  subroutine iter_next(this)
    class(string_set_iterator), intent(inout) :: this
    if (associated(this%item)) this%item => this%item%next
  end subroutine

  !! Returns true if the iterator has reached the end; that is, it has
  !! gone past the last element of the set.
  pure logical function iter_at_end(this)
    class(string_set_iterator), intent(in) :: this
    iter_at_end = .not.associated(this%item)
  end function

  !! Returns the current element.
  function iter_string(this)
    class(string_set_iterator), intent(in) :: this
    character(:), allocatable :: iter_string
    iter_string = this%item%string
  end function

end module string_set_type
