!!
!! UMAP
!!
!! This module defines a simple associative array (or map) data structure that
!! stores (key, value) pairs. The keys are unique character strings and the
!! values are CLASS(*) objects. This implementation is derived from the
!! MAP_ANY_TYPE module from the MIT-licensed Petaca library.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!! December 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines the uMAP derived type that is an associative
!! array with character keys and VECTOR_FUNC class values. It has the following
!! type bound procedures.
!!
!!  INSERT(KEY, VALUE) adds the specified key and associated value to the map.
!!    If the mapping already exists, its value is replaced with the specifed
!!    one. VALUE is an allocatable variable. Its allocation is taken by the map
!!    and it is returned unallocated.
!!
!!  REMOVE(KEY) removes the specified key from the map and deallocates the
!!    associated value.  If the mapping does not exist, the map is unchanged.
!!
!!  MAPPED(KEY) returns the value .TRUE. if a mapping for the specified key
!!    exists; otherwise it returns .FALSE.
!!
!!  LOOKUP(KEY, VALUE) returns a CLASS(*) pointer to the mapped value for the
!!    specified key if it exists; otherwise a null pointer is returned.
!!
!!  CLEAR() removes all elements from the map, leaving it with a size of 0.
!!

module umap_type

  implicit none
  private

  type :: list_item
    character(:), allocatable :: key
    class(*), allocatable :: value
    type(list_item), pointer :: next => null(), prev => null()
  contains
    final :: list_item_delete
  end type

  type, public :: umap
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: insert
    procedure :: remove
    procedure :: lookup
    procedure :: mapped
    procedure :: clear
    procedure :: copy
    final :: umap_delete
  end type

  type, public :: umap_iterator
    private
    class(list_item), pointer :: item => null()
  contains
    procedure :: next => iter_next
    procedure :: at_end => iter_at_end
    procedure :: name => iter_name
    procedure :: get_value => iter_get_value
    procedure :: value_ref => iter_value_ref
  end type

  !! User-defined UMAP_ITERATOR structure constructor
  interface umap_iterator
    procedure umap_begin
  end interface

contains

  !! Final procedure for UMAP objects.
  subroutine umap_delete(this)
    type(umap), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine

  !! Final procedure for LIST_ITEM objects.  This recursively follows the
  !! NEXT pointer. When deallocating a linked-list structure only the root
  !! needs to be explicitly deallocated. When the desire is to deallocate a
  !! single LIST_ITEM object, first nullify the NEXT point to prevent the
  !! recursive finalization from possibly deallocating more than it should.

  recursive subroutine list_item_delete(this)
    type(list_item), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  !!!! AUXILLARY ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns a pointer to a new initialized (but unlinked) LIST_ITEM.
  function new_list_item(key, value)
    character(*), intent(in) :: key
    class(*), allocatable, intent(inout) :: value
    type(list_item), pointer :: new_list_item
    allocate(new_list_item)
    new_list_item%key = key
    call move_alloc(value, new_list_item%value)
    new_list_item%prev => new_list_item
  end function

  !! Returns a pointer to the LIST_ITEM having the specified key,
  !! or a null pointer of none was found.
  function find_list_item(this, key) result(item)
    class(umap), intent(in) :: this
    character(*), intent(in) :: key
    type(list_item), pointer :: item
    item => this%first
    do while (associated(item))
      if (item%key == key) exit
      item => item%next
    end do
  end function

  !! Blindly links the given LIST_ITEM (as made by NEW_LIST_ITEM) to the end
  !! of the list; it does not check that the key is unique (someone else must).
  subroutine append_list_item(this, item)
    class(umap), intent(inout) :: this
    type(list_item), pointer, intent(in) :: item
    type(list_item), pointer :: tail
    if (associated(this%first)) then
      tail => this%first%prev
      tail%next => item
      item%prev => tail
      this%first%prev => item
    else
      item%prev => item
      this%first => item
    end if
  end subroutine

  !!!! UMAP TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns a CLASS(*) pointer to the mapped value for KEY,
  !! or a null pointer if KEY is not mapped.
  subroutine lookup(this, key, value)
    class(umap), intent(in) :: this
    character(*), intent(in) :: key
    !class(*), allocatable, intent(out) :: value
    !type(list_item), pointer :: item
    !item => find_list_item(this, key)
    !if (associated(item)) allocate(value, source=item%value)
    class(*), pointer :: value
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    value => null()
    if (associated(item)) value => item%value
  end subroutine

  !! Inserts the (KEY, VALUE) pair into the map.  If the mapping already
  !! exists its value is replaced with the specified value.
  subroutine insert(this, key, value)
    class(umap), intent(inout) :: this
    character(*), intent(in) :: key
    class(*), allocatable, intent(inout) :: value
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    if (associated(item)) then
      call move_alloc(value, item%value)
    else
      call append_list_item(this, new_list_item(key, value))
    end if
  end subroutine

  !! Removes KEY from the map and deallocates the mapped value.
  !! If the mapping does not exist the map is unchanged.
  subroutine remove(this, key)
    class(umap), intent(inout) :: this
    character(*), intent(in) :: key
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    if (associated(item)) then
      if (associated(item%prev, item)) then ! single item list
        this%first => null()
      else if (associated(this%first, item)) then ! first item of multiple
        this%first => item%next
        item%next%prev => item%prev
      else if (.not.associated(item%next)) then ! last item of multiple
        item%prev%next => item%next
        this%first%prev => item%prev
      else ! interior item of multiple
        item%prev%next => item%next
        item%next%prev => item%prev
      end if
      item%next => null() ! stop recursive finalization when item is deallocated
      deallocate(item)
    end if
  end subroutine

  !! Removes all elements from the map.
  subroutine clear(this)
    class(umap), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine

  !! Returns true if a mapping for KEY exists; otherwise returns false.
  logical function mapped(this, key)
    class(umap), intent(in) :: this
    character(*), intent(in) :: key
    mapped = associated(find_list_item(this, key))
  end function

  !! Copy the elements of map SRC to DEST.
  subroutine copy(src, dest)
    class(umap), intent(in) :: src
    type(umap), intent(inout) :: dest
    type(umap_iterator) :: iter
    class(*), allocatable :: val
    iter = umap_iterator(src)
    do while (.not.iter%at_end())
      call iter%get_value(val)
      call dest%insert(iter%name(), val)
      call iter%next
    end do
  end subroutine

!!!! UMAP_ITERATOR TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Defined UMAP_ITERATOR constructor that is positioned
  !! to the beginning element of the specified UMAP object.
  function umap_begin(map) result(iter)
    class(umap), intent(in) :: map
    type(umap_iterator) :: iter
    iter%item => map%first
  end function

  !! Advances iterator to the next element in the map.
  subroutine iter_next(this)
    class(umap_iterator), intent(inout) :: this
    if (associated(this%item)) this%item => this%item%next
  end subroutine

  !! Returns true if the iterator as reached the end; that is, it has
  !! gone past the last element of the map.
  pure logical function iter_at_end(this)
    class(umap_iterator), intent(in) :: this
    iter_at_end = .not.associated(this%item)
  end function

  !! Returns the name of the current element.
  function iter_name(this)
    class(umap_iterator), intent(in) :: this
    character(:), allocatable :: iter_name
    iter_name = this%item%key
  end function

  !! Returns a copy of the value for the current element.
  subroutine iter_get_value(this, val)
    class(umap_iterator), intent(in) :: this
    class(*), allocatable, intent(out) :: val
    if (allocated(this%item%value)) allocate(val, source=this%item%value)
  end subroutine

  !! Returns a reference to the value for the current element.
  function iter_value_ref(this) result(ref)
    class(umap_iterator), intent(in) :: this
    class(*), pointer :: ref
    ref => this%item%value
  end function

end module umap_type
