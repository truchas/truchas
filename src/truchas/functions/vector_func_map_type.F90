!!
!! VECTOR_FUNC_MAP
!!
!! This module defines an associative array (or map) data structure that stores
!! (key, value) pairs. The keys are unique character strings and the values are
!! VECTOR_FUNC objects.  This implementation is derived from the MAP_ANY_TYPE
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
!! This module defines the VECTOR_FUNC_MAP derived type that is an associative
!! array with character keys and VECTOR_FUNC class values. It has the following
!! type bound procedures.
!!
!!  INSERT(KEY, VALUE) adds the specified key and associated value to the map.
!!    If the mapping already exists, its value is replaced with the specifed
!!    one.  VALUE is an allocatable variable of class VECTOR_FUNC.  Its
!!    allocation is taken by the map and is returned unallocated.
!!
!!  REMOVE(KEY) removes the specified key from the map and deallocates the
!!    associated value.  If the mapping does not exist, the map is unchanged.
!!
!!  MAPPED(KEY) returns the value .TRUE. if a mapping for the specified key
!!    exists; otherwise it returns .FALSE.
!!
!!  LOOKUP(KEY, VALUE) returns the mapped value for the specified key.  VALUE
!!    is an allocatable variable of class VECTOR_FUNC.  It is allocated and
!!    assigned a copy of the mapped value if it exists; otherwise it is
!!    is returned unallocated.
!!
!!  CLEAR() removes all elements from the map, leaving it with a size of 0.
!!
!! NB: The polymorphic VECTOR_FUNC class values in the interface are all
!! allocatable.  The functions inserted into the map are handed off to the
!! the map using the MOVE_ALLOC intrinsic subroutine; no copies are made.
!! On the other hand, the functions returned by LOOKUP are copies of the
!! stored value as created by sourced-allocation.  These are shallow copies.
!! For pointer components this means that a copy of the pointer is made
!! but not a copy of its target; the original pointer and its copy will
!! have the same target.  Currently, none of the extensions of VECTOR_FUNC
!! have pointer components (except for DL_VECTOR_FUNC which holds a C_PTR
!! component that is the handle to the library -- here a shallow copy is
!! fine), so these are deep copies.  But this may change in the future and
!! these copies may not be what is needed.
!!

module vector_func_map_type

  use vector_func_class
  implicit none
  private

  public :: dump

  type :: list_item
    character(:), allocatable :: key
    class(vector_func), allocatable :: value
    type(list_item), pointer :: next => null(), prev => null()
  contains
    final :: list_item_delete
  end type list_item

  type, public :: vector_func_map
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: insert
    procedure :: remove
    procedure :: lookup
    procedure :: mapped
    procedure :: clear
    final :: vector_func_map_delete
  end type vector_func_map

contains

  !! Final procedure for MAP_ANY objects.
  subroutine vector_func_map_delete (this)
    type(vector_func_map), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine vector_func_map_delete

  !! Final procedure for LIST_ITEM objects.  This recursively follows the
  !! NEXT pointer.  When deallocating a linked-list structure only the root
  !! needs to be explicitly deallocated.  When the desire is to deallocate a
  !! single LIST_ITEM object, first nullify the NEXT point to prevent the
  !! recursive finalization from possibly deallocating more than it should.
  recursive subroutine list_item_delete (this)
    type(list_item), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine list_item_delete

  !!!! AUXILLARY ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns a pointer to a new initialized (but unlinked) LIST_ITEM.
  function new_list_item (key, value)
    character(*), intent(in) :: key
    class(vector_func), allocatable, intent(inout) :: value
    type(list_item), pointer :: new_list_item
    allocate(new_list_item)
    new_list_item%key = key
    call move_alloc (value, new_list_item%value)
    new_list_item%prev => new_list_item
  end function new_list_item

  !! Returns a pointer to the LIST_ITEM having the specified key,
  !! or a null pointer of none was found.
  function find_list_item (this, key) result (item)
    class(vector_func_map), intent(in) :: this
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
  subroutine append_list_item (this, item)
    class(vector_func_map), intent(inout) :: this
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

  !!!! VECTOR_FUNC_MAP TYPE-BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns a CLASS(*) pointer to the mapped value for KEY,
  !! or a null pointer if KEY is not mapped.
  subroutine lookup (this, key, value)
    class(vector_func_map), intent(in) :: this
    character(*), intent(in) :: key
    class(vector_func), allocatable, intent(out) :: value
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    if (associated(item)) allocate(value, source=item%value)
  end subroutine lookup

  !! Inserts the (KEY, VALUE) pair into the map.  If the mapping already
  !! exists its value is replaced with the specified value.
  subroutine insert (this, key, value)
    class(vector_func_map), intent(inout) :: this
    character(*), intent(in) :: key
    class(vector_func), allocatable, intent(inout) :: value
    type(list_item), pointer :: item
    item => find_list_item(this, key)
    if (associated(item)) then
      call move_alloc (value, item%value)
    else
      call append_list_item(this, new_list_item(key, value))
    end if
  end subroutine insert

  !! Removes KEY from the map and deallocates the mapped value.
  !! If the mapping does not exist the map is unchanged.
  subroutine remove (this, key)
    class(vector_func_map), intent(inout) :: this
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
  subroutine clear (this)
    class(vector_func_map), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine clear

  !! Returns true if a mapping for KEY exists; otherwise returns false.
  logical function mapped (this, key)
    class(vector_func_map), intent(in) :: this
    character(*), intent(in) :: key
    mapped = associated(find_list_item(this, key))
  end function mapped

  subroutine dump (this)
    type(vector_func_map), intent(in) :: this
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

end module vector_func_map_type
