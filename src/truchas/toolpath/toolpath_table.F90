!!
!! TOOLPATH_TABLE
!!
!! This module manages the global table of TOOLPATH objects.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL INSERT_TOOLPATH(NAME, TP) adds the TOOLPATH object TP to the table
!!    and associates it with NAME.  if a toolpath with this name already exists,
!!    its associated object is replaced by this one; the original object is
!!    deallocated.  TP is allocatable and its allocation is moved to the table
!!    and it is returned unallocated.
!!
!!  TOOLPATH_PTR(NAME) returns a pointer to the TOOLPATH object associated with
!!    NAME; a null pointer is returned if no such toolpath exists.
!!
!!  KNOWN_TOOLPATH(NAME) returns the value .TRUE. if the table contains a
!!    toolpath with the given NAME; otherwise it returns .FALSE.
!!
!! IMPLEMENTATION NOTES
!!
!!  The table is just a map (or associative array) of (name, toolpath) pairs.
!!  Rather than make yet-another-copy of a map data structure implementation
!!  with TOOLPATH type values (no support for generic programming in Fortran),
!!  I opted to leverage the MAP_ANY type which holds CLASS(*) values.  This
!!  container holds shallow copies of values inserted, and likewise returns
!!  shallow copies of stored values.  Unfortunately this isn't the semantics
!!  desired in this case.  Toolpath objects are large and so we want to move
!!  inserted values into the container and return pointers to the values.
!!  To get this behavior:
!!
!!  1. We wrap the pointer in a type OUTER_BOX.  A shallow copy of a variable
!!     of this type copies the pointer, not its target, which is what we want.
!!
!!  2. The allocation of allocatable variable can only be moved to another
!!     allocatable variable.  Thus we wrap the allocatable variable in a type
!!     INNNER_BOX, and have OUTER_BOX wrap a pointer to this type.
!!

#include "f90_assert.fpp"

module toolpath_table

  use toolpath_type
  use map_any_type
  implicit none
  private

  public :: insert_toolpath, toolpath_ptr, known_toolpath

  type(map_any) :: table

  !! See Note 1
  type :: outer_box
    type(inner_box), pointer :: ibox
  contains
    final :: outer_box_delete
  end type

  !! See Note 2
  type :: inner_box
    type(toolpath), allocatable :: tp
  end type

contains

  subroutine outer_box_delete(this)
    type(outer_box), intent(inout) :: this
    if (associated(this%ibox)) deallocate(this%ibox)
  end subroutine outer_box_delete

  subroutine insert_toolpath(name, tp)
    character(*), intent(in) :: name
    type(toolpath), allocatable, intent(inout) :: tp
    type(outer_box) :: obox
    allocate(obox%ibox)
    call move_alloc(tp, obox%ibox%tp)
    call table%insert(name, obox)
    obox%ibox => null() ! so finalization of obox will not deallocate %ibox
  end subroutine insert_toolpath

  function toolpath_ptr(name) result(tp)
    character(*), intent(in) :: name
    type(toolpath), pointer :: tp
    class(*), pointer :: uptr
    uptr => table%value(name)
    if (associated(uptr)) then
      select type (uptr)
      type is (outer_box)
        tp => uptr%ibox%tp
      class default
        INSIST(.false.)
      end select
    else
      tp => null()
    end if
  end function toolpath_ptr

  logical function known_toolpath(name)
    character(*), intent(in) :: name
    known_toolpath = table%mapped(name)
  end function known_toolpath

end module toolpath_table
