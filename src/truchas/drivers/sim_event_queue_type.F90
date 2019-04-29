!!
!! SIM_EVENT_QUEUE_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! SYNOPSIS
!!
!!  QUEUE CREATION:
!!
!!    type, extends(event_action) :: my_action1
!!    type, extends(event_action) :: my_action2
!!    ...
!!
!!    type(my_action1) :: foo
!!    type(my_action2) :: bar
!!    foo = ...
!!    bar = ...
!!
!!    type(sim_event_queue) :: eventq
!!    call eventq%add_event(2.0_r8, foo, rank=3)
!!    call eventq%add_event(1.0_r8, bar)  ! default rank is 0
!!    call eventq%add_event(2.0_r8, bar, rank=1)
!!    ...
!!
!!  QUEUE CONSUMPTION DURING A SIMULATION:
!!
!!    t = eventq%next_time()
!!    integrate until reaching time t (exactly)
!!
!!    type(action_list), allocatable :: actions
!!    call eventq%pop_actions(actions)  ! this pops the time and actions from the queue
!!
!!    do
!!      class(event_action), allocatable :: action
!!      call actions%next_action(action) ! this removes action from the list
!!      if (.not.allocated(action)) exit
!!      select type (action)
!!      type is (my_action1)
!!        ...
!!      type is (my_action2)
!!        ...
!!      end select
!!    end do

module sim_event_queue_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: event_action
  end type event_action

  type, public :: action_list
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: get_next_action
    procedure, private :: add_action
    generic :: assignment(=) => action_list_copy
    procedure, private :: action_list_copy
    final :: action_list_delete
  end type

  type :: list_item
    class(event_action), allocatable :: action
    integer :: rank = 0
    type(list_item), pointer :: next => null()
  end type

  type, public :: sim_event_queue
    private
    type(queue_item), pointer :: top => null()
    real(r8) :: h = 0.0_r8 ! time resolution
  contains
    procedure :: add_event
    procedure :: next_time
    procedure :: pop_actions
    procedure :: fast_forward
    procedure :: is_empty
    procedure :: set_time_resolution
    generic :: assignment(=) => copy
    procedure, private :: copy
    final :: sim_event_queue_delete
  end type sim_event_queue

  type :: queue_item
    real(r8) :: time
    type(action_list), allocatable :: actions
    type(queue_item), pointer :: next => null()
  end type queue_item

contains

  !! Final subroutine for SIM_EVENT_QUEUE objects.
  subroutine sim_event_queue_delete(this)
    type(sim_event_queue), intent(inout) :: this
    type(queue_item), pointer :: item
    do while (associated(this%top))
      item => this%top
      this%top => item%next
      deallocate(item)
    end do
  end subroutine sim_event_queue_delete

  pure logical function is_empty(this)
    class(sim_event_queue), intent(in) :: this
    is_empty = .not.associated(this%top)
  end function is_empty

  pure function next_time(this)
    class(sim_event_queue), intent(in) :: this
    real(r8) :: next_time
    next_time = this%top%time
  end function next_time

  pure subroutine add_event(this, time, action, rank)
    class(sim_event_queue), intent(inout) :: this
    real(r8), intent(in) :: time
    class(event_action), intent(in) :: action
    integer, intent(in), optional :: rank
    type(queue_item), pointer :: item, new_item
    call find_queue_item(this, time, item)
    if (associated(item)) then
      if (associated(item%next)) then
        if (time <= min(item%time+this%h,(item%time+item%next%time)/2)) then
          call item%actions%add_action(action, rank)
          return
        else if (time >= max(item%next%time-this%h,(item%time+item%next%time)/2)) then
          call item%next%actions%add_action(action, rank)
          return
        end if
      else if (time <= item%time+this%h) then
        call item%actions%add_action(action, rank)
        return
      end if
      !! Insert new after item
      allocate(new_item)
      new_item%time = time
      allocate(new_item%actions)
      call new_item%actions%add_action(action, rank)
      new_item%next => item%next
      item%next => new_item
    else
      if (associated(this%top)) then
        if (time >= this%top%time-this%h) then
          call this%top%actions%add_action(action, rank)
          return
        end if
      end if
      !! Insert new at the top
      allocate(new_item)
      new_item%time = time
      allocate(new_item%actions)
      call new_item%actions%add_action(action, rank)
      new_item%next => this%top
      this%top => new_item
    end if
  end subroutine add_event

  subroutine pop_actions(this, actions)
    class(sim_event_queue), intent(inout) :: this
    type(action_list), allocatable, intent(out) :: actions
    type(queue_item), pointer :: item
    if (associated(this%top)) then
      call move_alloc(this%top%actions, actions)
      item => this%top
      this%top => item%next
      deallocate(item)
    end if
  end subroutine pop_actions

  subroutine fast_forward(this, time)
    class(sim_event_queue), intent(inout) :: this
    real(r8), intent(in) :: time
    type(queue_item), pointer :: item
    do while (associated(this%top))
      if (this%top%time > time) exit
      item => this%top
      this%top => item%next
      deallocate(item)
    end do
  end subroutine fast_forward

  subroutine set_time_resolution(this, dt)
    class(sim_event_queue), intent(inout) :: this
    real(r8), intent(in) :: dt
    this%h = dt
    !TODO: should scan the current queue to ensure that none of the
    !events are spaced closer that DT.
  end subroutine set_time_resolution

  subroutine copy(lhs, rhs)
    class(sim_event_queue), intent(inout) :: lhs  ! inout in case lhs is rhs
    class(sim_event_queue), intent(in) :: rhs
    if (associated(lhs%top, rhs%top)) return ! lhs and rhs are same
    call sim_event_queue_delete(lhs)
    if (associated(rhs%top)) lhs%top => item_copy(rhs%top)
    lhs%h = rhs%h
  end subroutine copy

  recursive function item_copy(item) result(copy)
    type(queue_item), intent(in) :: item
    type(queue_item), pointer :: copy
    allocate(copy, mold=item)
    copy%time = item%time
    if (allocated(item%actions)) then
      allocate(copy%actions)
      copy%actions = item%actions
    end if
    if (associated(item%next)) copy%next => item_copy(item%next)
  end function item_copy

  !! This auxiliary procedure returns a pointer to the QUEUE_ITEM that is the
  !! insertion point in the ordered list for a new item with the given time.
  !! Note that the time of the QUEUE_ITEM may equal the given time.  A null
  !! pointer is returned if the insertion point is at the head of the list.

  pure subroutine find_queue_item(this, time, item)
    class(sim_event_queue), intent(inout) :: this
    real(r8), intent(in) :: time
    type(queue_item), intent(out), pointer :: item
    item => this%top
    if (.not.associated(item)) return
    if (time < item%time) then
      item => null()
      return
    end if
    do while (associated(item%next))
      if (time < item%next%time) return
      item => item%next
    end do
  end subroutine find_queue_item

  !!!! ACTION_LIST TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Final subroutine for ACTION_LIST objects.
  subroutine action_list_delete(this)
    type(action_list), intent(inout) :: this
    type(list_item), pointer :: item => null()
    do while (associated(this%first))
      item => this%first
      this%first => item%next
      deallocate(item)
    end do
  end subroutine action_list_delete

  subroutine get_next_action(this, action)
    class(action_list), intent(inout) :: this
    class(event_action), allocatable, intent(out) :: action
    type(list_item), pointer :: item
    if (associated(this%first)) then
      call move_alloc(this%first%action, action)
      item => this%first
      this%first => item%next
      deallocate(item)
    end if
  end subroutine get_next_action

  pure subroutine add_action(this, action, rank)
    class(action_list), intent(inout) :: this
    class(event_action), intent(in) :: action
    integer, intent(in), optional :: rank
    type(list_item), pointer :: item, new_item
    allocate(new_item)
    allocate(new_item%action, source=action)
    if (present(rank)) new_item%rank = rank
    call find_list_item(this, new_item%rank, item)
    if (associated(item)) then
      new_item%next => item%next
      item%next => new_item
    else
      new_item%next => this%first
      this%first => new_item
    end if
  end subroutine add_action

  subroutine action_list_copy(lhs, rhs)
    class(action_list), intent(inout) :: lhs
    class(action_list), intent(in) :: rhs
    if (associated(lhs%first, rhs%first)) return  ! lhs and rhs are same
    call action_list_delete(lhs)
    if (associated(rhs%first)) lhs%first => list_item_copy(rhs%first)
  end subroutine action_list_copy

  recursive function list_item_copy(item) result(copy)
    type(list_item), intent(in) :: item
    type(list_item), pointer :: copy
    allocate(copy, source=item)
    if (associated(item%next)) copy%next => list_item_copy(item%next)
  end function list_item_copy

  !! This auxiliary procedure returns a pointer to the LIST_ITEM that is the
  !! insertion point in the ordered list for a new item with the given RANK.
  !! The insertion point is the last item having equal or lower RANK. A null
  !! pointer is returned if the insertion point is at the head of the list.

  pure subroutine find_list_item(this, rank, item)
    class(action_list), intent(inout) :: this
    integer, intent(in) :: rank
    type(list_item), intent(out), pointer :: item
    item => this%first
    if (.not.associated(item)) return
    if (rank < item%rank) then
      item => null()
      return
    end if
    do while (associated(item%next))
      if (rank < item%next%rank) return
      item => item%next
    end do
  end subroutine find_list_item

end module sim_event_queue_type
