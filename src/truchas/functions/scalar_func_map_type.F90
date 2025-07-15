!!
!! SCALAR_FUNC_MAP
!!
!! This module defines an associative array (or map) data structure that stores
!! (key, value) pairs. The keys are unique character strings and the values are
!! SCALAR_FUNC objects.  This implementation is derived from the MAP_ANY_TYPE
!! module from the MIT-licensed Petaca library.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines the SCALAR_FUNC_MAP derived type that is an associative
!! array with character keys and SCALAR_FUNC class values. It has the following
!! type bound procedures.
!!
!!  INSERT(KEY, VALUE) adds the specified key and associated value to the map.
!!    If the mapping already exists, its value is replaced with the specifed
!!    one.  VALUE is an allocatable variable of class SCALAR_FUNC.  Its
!!    allocation is taken by the map and is returned unallocated.
!!
!!  REMOVE(KEY) removes the specified key from the map and deallocates the
!!    associated value.  If the mapping does not exist, the map is unchanged.
!!
!!  MAPPED(KEY) returns the value .TRUE. if a mapping for the specified key
!!    exists; otherwise it returns .FALSE.
!!
!!  LOOKUP(KEY, VALUE) returns the mapped value for the specified key.  VALUE
!!    is an allocatable variable of class SCALAR_FUNC.  It is allocated and
!!    assigned a copy of the mapped value if it exists; otherwise it is
!!    is returned unallocated.
!!
!!  CLEAR() removes all elements from the map, leaving it with a size of 0.
!!
!! NB: The polymorphic SCALAR_FUNC class values in the interface are all
!! allocatable.  The functions inserted into the map are handed off to the
!! the map using the MOVE_ALLOC intrinsic subroutine; no copies are made.
!! On the other hand, the functions returned by LOOKUP are copies of the
!! stored value as created by sourced-allocation.  These are shallow copies.
!! For pointer components this means that a copy of the pointer is made
!! but not a copy of its target; the original pointer and its copy will
!! have the same target.  Currently, none of the extensions of SCALAR_FUNC
!! have pointer components (except for DL_SCALAR_FUNC which holds a C_PTR
!! component that is the handle to the library -- here a shallow copy is
!! fine), so these are deep copies.  But this may change in the future and
!! these copies may not be what is needed.
!!

module scalar_func_map_type

  use scalar_func_class
  implicit none
  private

  type :: list_item
    character(:), allocatable :: key
    class(scalar_func), allocatable :: value
    type(list_item), pointer :: next => null(), prev => null()
  contains
    final :: list_item_delete
  end type list_item

  type, public :: scalar_func_map
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: insert
    procedure :: remove
    procedure :: lookup
    procedure :: mapped
    procedure :: clear
    procedure :: copy
    final :: scalar_func_map_delete
    procedure :: dump
  end type scalar_func_map

  type, public :: scalar_func_map_iterator
    private
    class(list_item), pointer :: item => null()
  contains
    procedure :: next => iter_next
    procedure :: at_end => iter_at_end
    procedure :: name => iter_name
    procedure :: get_func => iter_get_func
  end type scalar_func_map_iterator

  !! User-defined SCALAR_FUNC_MAP_ITERATOR structure constructor
  interface scalar_func_map_iterator
    procedure scalar_func_map_begin
  end interface

contains

  !! Final procedure for MAP_ANY objects.
  subroutine scalar_func_map_delete(this)
    type(scalar_func_map), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine scalar_func_map_delete

  !! Final procedure for LIST_ITEM objects.  This recursively follows the
  !! NEXT pointer.  When deallocating a linked-list structure only the root
  !! needs to be explicitly deallocated.  When the desire is to deallocate a
  !! single LIST_ITEM object, first nullify the NEXT point to prevent the
  !! recursive finalization from possibly deallocating more than it should.
  recursive subroutine list_item_delete(this)
    type(list_item), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine list_item_delete

  !!!! AUXILLARY ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns a pointer to a new initialized (but unlinked) LIST_ITEM.
  function new_list_item(key, value)
    character(*), intent(in) :: key
    class(scalar_func), allocatable, intent(inout) :: value
    type(list_item), pointer :: new_list_item
    allocate(new_list_item)
    new_list_item%key = key
    call move_alloc(value, new_list_item%value)
    new_list_item%prev => new_list_item
  end function new_list_item

  !! Returns a pointer to the LIST_ITEM having the specified key,
  !! or a null pointer of none was found.
  function find_list_item(this, key) result(item)
    class(scalar_func_map), intent(in) :: this
    character(*), intent(in) :: key
    type(list_item), pointer :: item
    item => this%first
    do while (associated(item))
      if (item%key == key) exit
      item => item%next
    end do
  end function find_list_item

  !! Blindly links the given LIST_ITEM (as made by NEW_LIST_ITEM) to the end
  !! of the list; it does not check that the key is unique (someone else must).
  subroutine append_list_item(this, item)
    class(scalar_func_map), intent(inout) :: this
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
  end subroutine append_list_item

  !!!! SCALAR_FUNC_MAP TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns a CLASS(*) pointer to the mapped value for KEY,
  !! or a null pointer if KEY is not mapped.
  subroutine lookup(this, key, value)
    class(scalar_func_map), intent(in) :: this
    character(*), intent(in) :: key
    class(scalar_func), allocatable, intent(out) :: value
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    if (associated(item)) allocate(value, source=item%value)
  end subroutine lookup

  !! Inserts the (KEY, VALUE) pair into the map.  If the mapping already
  !! exists its value is replaced with the specified value.
  subroutine insert(this, key, value)
    class(scalar_func_map), intent(inout) :: this
    character(*), intent(in) :: key
    class(scalar_func), allocatable, intent(inout) :: value
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    if (associated(item)) then
      call move_alloc(value, item%value)
    else
      call append_list_item(this, new_list_item(key, value))
    end if
  end subroutine insert

  !! Removes KEY from the map and deallocates the mapped value.
  !! If the mapping does not exist the map is unchanged.
  subroutine remove(this, key)
    class(scalar_func_map), intent(inout) :: this
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
  end subroutine remove

  !! Removes all elements from the map.
  subroutine clear(this)
    class(scalar_func_map), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine clear

  !! Returns true if a mapping for KEY exists; otherwise returns false.
  logical function mapped(this, key)
    class(scalar_func_map), intent(in) :: this
    character(*), intent(in) :: key
    mapped = associated(find_list_item(this, key))
  end function mapped

  !! Copy the elements of map SRC to DEST.
  subroutine copy(src, dest)
    class(scalar_func_map), intent(in) :: src
    type(scalar_func_map), intent(inout) :: dest
    type(scalar_func_map_iterator) :: iter
    class(scalar_func), allocatable :: f
    iter = scalar_func_map_iterator(src)
    do while (.not.iter%at_end())
      call iter%get_func(f)
      call dest%insert(iter%name(), f)
      call iter%next
    end do
  end subroutine

  !!!! SCALAR_FUNC_MAP_ITERATOR TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Defined SCALAR_FUNC_MAP_ITERATOR constructor that is positioned
  !! to the beginning element of the specified SCALAR_FUNC_MAP object.
  function scalar_func_map_begin(map) result(iter)
    class(scalar_func_map), intent(in) :: map
    type(scalar_func_map_iterator) :: iter
    iter%item => map%first
  end function scalar_func_map_begin

  !! Advances the iterator to the next element in the map.
  subroutine iter_next(this)
    class(scalar_func_map_iterator), intent(inout) :: this
    if (associated(this%item)) this%item => this%item%next
  end subroutine iter_next

  !! Returns true if the iterator has reached the end; that is, it has
  !! gone past the last element of the map.
  pure logical function iter_at_end(this)
    class(scalar_func_map_iterator), intent(in) :: this
    iter_at_end = .not.associated(this%item)
  end function iter_at_end

  !! Returns the name for the current element.
  function iter_name(this)
    class(scalar_func_map_iterator), intent(in) :: this
    character(:), allocatable :: iter_name
    iter_name = this%item%key
  end function iter_name

  !! Returns a copy of the scalar function for the current element.
  subroutine iter_get_func(this, func)
    class(scalar_func_map_iterator), intent(in) :: this
    class(scalar_func), allocatable, intent(out) :: func
    if (associated(this%item)) allocate(func, source=this%item%value)
  end subroutine iter_get_func

  subroutine dump(this)
    class(scalar_func_map), intent(in) :: this
    type(list_item), pointer :: item, last, tail
    item => this%first
    if (associated(item)) then
      write(*,'(a)',advance='no') trim(item%key)
      tail => item%prev
      last => item
      item => item%next
      do while (associated(item))
        if (.not.associated(last,item%prev)) write(*,'(1x,a)',advance='no') 'X<'
        write(*,'(1x,a)',advance='no') trim(item%key)
        last => item
        item => item%next
      end do
      if (.not.associated(last,tail)) write(*,'(1x,a)',advance='no') 'XT'
      write(*,*)
    end if
  end subroutine

end module scalar_func_map_type
