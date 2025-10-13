!!
!! VECTOR_FUNC_MAP
!!
!! This module defines an associative array (or map) data structure that stores
!! (key, value) pairs. The keys are unique character strings and the values are
!! VECTOR_FUNC objects.
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

  use umap_type
  use vector_func_class
  implicit none
  private

  type, public :: vector_func_map
    private
    type(umap) :: map
  contains
    procedure :: insert
    procedure :: remove
    procedure :: lookup
    procedure :: mapped
    procedure :: clear
  end type

contains

  subroutine lookup(this, key, value)
    class(vector_func_map), intent(in) :: this
    character(*), intent(in) :: key
    class(vector_func), allocatable, intent(out) :: value
    class(*), pointer :: uptr
    call this%map%lookup(key, uptr)
    if (.not.associated(uptr)) return
    select type (uptr)
    class is (vector_func)
      allocate(value, source=uptr)
    end select
  end subroutine

  subroutine insert(this, key, value)
    class(vector_func_map), intent(inout) :: this
    character(*), intent(in) :: key
    class(vector_func), allocatable, intent(inout) :: value
    class(*), allocatable :: tmp
    call move_alloc(value, tmp)
    call this%map%insert(key, tmp) ! tmp returned unallocated
  end subroutine

  subroutine remove(this, key)
    class(vector_func_map), intent(inout) :: this
    character(*), intent(in) :: key
    call this%map%remove(key)
  end subroutine

  logical function mapped(this, key)
    class(vector_func_map), intent(in) :: this
    character(*), intent(in) :: key
    mapped = this%map%mapped(key)
  end function

  subroutine clear(this)
    class(vector_func_map), intent(inout) :: this
    call this%map%clear
  end subroutine

end module vector_func_map_type
