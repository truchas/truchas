!!
!! TOOLHEAD_TABLE_TYPE
!!
!! This module defines the TOOLHEAD_TABLE derived type which stores a table
!! of TOOLHEAD objects that are referenced by name.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! IMPLEMENTATION NOTES
!!
!!  The table is just a map (or associative array) of (name, toolhead) pairs.
!!  Rather than make yet-another-copy of a map data structure implementation
!!  with TOOLHEAD type values (no support for generic programming in Fortran),
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

module toolhead_table_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use map_any_type
  use toolhead_type
  use sim_event_queue_type, only: sim_event_queue
  implicit none
  private

  type, extends(map_any), public :: toolhead_table
  contains
    procedure :: add_toolhead
    procedure :: known_toolhead
    procedure :: set_initial_state
    procedure :: toolhead_ptr
    procedure :: add_events
  end type

  !! See Note 1
  type :: outer_box
    type(inner_box), pointer :: ibox
  contains
    final :: outer_box_delete
  end type

  !! See Note 2
  type :: inner_box
    type(toolhead), allocatable :: th
  end type

contains

  subroutine outer_box_delete(this)
    type(outer_box), intent(inout) :: this
    if (associated(this%ibox)) deallocate(this%ibox)
  end subroutine

  subroutine add_toolhead(this, name, th)
    class(toolhead_table), intent(inout) :: this
    character(*), intent(in) :: name
    type(toolhead), allocatable, intent(inout) :: th
    type(outer_box) :: obox
    allocate(obox%ibox)
    call move_alloc(th, obox%ibox%th)
    call this%insert(name, obox)
    obox%ibox => null() ! so finalization of obox will not deallocate %ibox
  end subroutine

  logical function known_toolhead(this, name)
    class(toolhead_table), intent(in) :: this
    character(*), intent(in) :: name
    known_toolhead = this%mapped(name)
  end function

  function toolhead_ptr(this, name) result(th)
    class(toolhead_table), intent(in) :: this
    character(*), intent(in) :: name
    type(toolhead), pointer :: th
    class(*), pointer :: uptr
    uptr => this%value(name)
    if (associated(uptr)) then
      select type (uptr)
      type is (outer_box)
        th => uptr%ibox%th
      class default
        INSIST(.false.)
      end select
    else
      th => null()
    end if
  end function

  subroutine set_initial_state(this, t)
    class(toolhead_table), intent(in) :: this
    real(r8), intent(in) :: t
    type(map_any_iterator) :: iter
    class(*), pointer :: uptr
    iter = map_any_iterator(this)
    do while (.not.iter%at_end())
      uptr => iter%value()
      select type (uptr)
      type is (outer_box)
        associate (th => uptr%ibox%th)
          call th%set_initial_state(t)
        end associate
      class default
        INSIST(.false.)
      end select
      call iter%next
    end do
  end subroutine


  subroutine add_events(this, eventq)
    use toolhead_event_type, only: add_toolhead_events
    class(toolhead_table), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    type(map_any_iterator) :: iter
    class(*), pointer :: uptr
    iter = map_any_iterator(this)
    do while (.not.iter%at_end())
      uptr => iter%value()
      select type (uptr)
      type is (outer_box)
        associate (th => uptr%ibox%th)
          call add_toolhead_events(th, eventq)
        end associate
      end select
      call iter%next
    end do
  end subroutine

end module toolhead_table_type
