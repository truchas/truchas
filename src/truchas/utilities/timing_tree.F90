!!===================================================================================
!!TIMING TREE
!! 
!!PURPOSE(S):
!! 
!!    Provide a system of timers that can be accessed through three simple 
!!    subroutines. 
!!    All timers are contained in a nested structure. When a timer is started for the
!!    first time, it is created as a child of the current timer, and then becomes the
!!    current timer. When a timer is stopped, the current timer becomes its parent. A
!!    timer can have multiple children, resulting in a tree like structure,  When the
!!    print_timers subroutine is called, all timers and times are printed in a tree
!!    structure:
!!
!!    >Timer A
!!    > >Timer B
!!    > > >Timer C
!!    > > >Timer D
!!    > >Timer E
!!
!!    All timers record both the elapsed wall clock time and cpu usage.
!!
!!CONTAINS:
!!    start_timer(key)      !Creates (or moves to) a child timer named key.
!!    stop_timer(key)       !Stops the running timer and moves to parent.
!!    write_timer_tree      !Print information for all existing timers.
!!
!!    custom_timer_tree     !Execute custom code for all existing timers
!!
!!    walk_tree             !Initializes a sequence of methods to retrieve timer
!!                          !information in a logical order.
!!    has_next              !Check to see if another timer is available. If so, move
!!                          !to it.
!!    next_timer            !Get the next timer and all its information.
!!
!!    flatten_timing_tree   !Return all the tree information stored in arrays.
!!===================================================================================

#include "f90_assert.fpp"

module timing_tree
implicit none
private
public :: start_timer, stop_timer, write_timing_tree, destroy_timing_tree, flatten_timing_tree

integer,  parameter, public :: MAX_KEY_LEN = 31

type :: timer
   character(len=MAX_KEY_LEN) :: key                     ! Name of the clock--Maximum 31 Characters
   integer                    :: wall_start = 0          ! The time when the wall clock was started.
   real                       :: wall_clock = 0.0        ! The running total for the wall clock.
   real                       :: cpu_start  = 0.0        ! The time when the cpu clock was started.
   real                       :: cpu_clock  = 0.0        ! The running total for the cpu clock.
   integer                    :: depth                   ! Depth of position in the tree.
   integer                    :: marker     = 0
   type(timer),    pointer    :: parent     => null()    ! Point to the parent timer
   type(timer),    pointer    :: sibling    => null()    ! Point to the next sibling timer
   type(timer),    pointer    :: child      => null()    ! Point to the child timer
end type timer

type(timer), save, pointer :: root_timer    => null() ! The root timer: Not actually used as a timer.
type(timer), save, pointer :: current_timer => null() ! Pointer to the current timer
integer,     save          :: maxdepth                ! Keep track of the meximum depth of tree
integer,     save          :: markers = 1


contains

!!==============================================
!!START_TIMER(key)
!!
!!PURPOSE(S):
!!    Set the current timer to be a timer with 
!!    name [key]. Either restart an already 
!!    existing timer with that name, or create a
!!    new timer with that name.
!!    Timers can optionally be given a marker
!!    for later access and partial tree
!!    printing.
!!
!!PARAMETER(S):
!!    key    !Name of the timer to be started.
!!    mark   !Optional place marker for partial
!!           !printing of the tree.
!!==============================================
subroutine start_timer(key, mark)
  character(len=*),          intent(in)  :: key
  integer,         optional, intent(out) :: mark
  type(timer),     pointer               :: new_timer

  ASSERT(len(key) <= MAX_KEY_LEN)

  !If first time, create root and add timer.
  if (.not.associated(current_timer)) then
     call add_new_timer(key)

  else

     !Search to see if a timer named [key] exists.
     new_timer => current_timer%child
     do while (associated(new_timer))
        if (new_timer%key == key) exit
        new_timer => new_timer%sibling
     end do
     
     !If [key] does not exist, create it.
     if (.not. associated(new_timer)) then 
        call add_new_timer(key)
     
     !Move to the new timer
     else 
        current_timer => new_timer
     end if

  end if

  !To be used if is passed a marker.
  if (present(mark)) then
     current_timer%marker = markers
     mark = markers
     markers = markers + 1
     print *, 'Marker: ', mark
  end if

  call system_clock(current_timer%wall_start)
  call cpu_time(current_timer%cpu_start)
end subroutine start_timer


!!==============================================
!!STOP_TIMER(key)
!!
!!PURPOSE(S):
!!    Stop the current timer and moves to the
!!    parent.
!!
!!PARAMETER(S):
!!    key    !Name of timer to be stopped.
!!==============================================
subroutine stop_timer(key)
  character(len=*), intent(in)    :: key
  integer                         :: wall_buffer
  real                            :: cpu_buffer
  
  !Check that key matches name of current timer.
  ASSERT(associated(current_timer))
  ASSERT(current_timer%key == key)
  
  !Get current time from FORTRAN subs.
  call system_clock(wall_buffer)
  call cpu_time(cpu_buffer)
  
  !Set the values of the current timer.
  current_timer%wall_clock = current_timer%wall_clock + (wall_buffer - current_timer%wall_start)
  current_timer%cpu_clock = current_timer%cpu_clock + (cpu_buffer - current_timer%cpu_start)

  !Move to the parent timer.
  current_timer => current_timer%parent 
end subroutine stop_timer


!!================================================
!!DESTROY_TIMING_TREE
!!
!!PURPOSE(S):
!!    Destroy all existing timers, 
!!
!!CONTAINS:
!!    printSubTimers  !Destroy all child timers
!!================================================
subroutine destroy_timing_tree()
  
  ASSERT(associated(current_timer))

  call destroySubTimers(root_timer)

  contains

  recursive subroutine destroySubTimers(root)
    type(timer), pointer :: root
    type(timer), pointer :: root_destroy
    type(timer), pointer :: temp
    type(timer), pointer :: temp2

    root_destroy => root

    !Start the same recursive subroutine
    temp => root_destroy%child
    do while (associated(temp))
       temp2 => temp%sibling
       call destroySubTimers(temp)
       temp => temp2
    end do

    deallocate(root)

  end subroutine destroySubTimers
end subroutine destroy_timing_tree




!!================================================
!!WRITE_TIMING_TREE
!!
!!PURPOSE(S):
!!    Print all created timers and their recorded
!!    times in an outline structure.    
!!    Optional paramater for partial tree printing
!!
!!CONTAINS:
!!    printSubTimers  !Loop through all available
!!                     !timers and print times.
!!================================================
subroutine write_timing_tree(file, mark)
  integer, optional, intent(in)   :: file
  integer, optional, intent(in)   :: mark
  
  character(len=(maxdepth*2) + 9) :: seperator
  character(len=8)                :: prefix
  integer                         :: unit
  

  ASSERT(associated(current_timer))

  if (present(file)) then
     unit = file
  else 
     unit = 6
  end if
     

  call accumulate_elapsed_time()

  !Initialize the spacer values.
  seperator = " "
  prefix = " "
  
  !Print the heading
  write(unit,fmt = "(a)") " "
  write(unit,fmt = "(4a)") prefix, " Timer Structure", seperator, "            WALL         CPU    "
  write(unit,fmt = "(4a)") prefix, " ---------------", seperator, "          ---------    ---------"
  

  !Start the recursive subroutine
  if (present(mark)) then
     call printSubTimers("", root_timer, maxdepth*2, .false.)
  else
     call printSubTimers("", root_timer, maxdepth*2, .true.)
  end if
  
  !Give the output a little spare room.
  write(unit,fmt = "(a)") " "

  contains

  recursive subroutine printSubTimers(indent, root, dep, print_data)
    character(len=*),   intent(in)  :: indent
    type(timer),        intent(in)  :: root
    integer,            intent(in)  :: dep
    logical,            intent(in)  :: print_data
    
    type(timer),        pointer     :: temp
    character(len=dep)              :: space
    logical                         :: param

    !Initialize the spacer
    space = " "

    param = print_data

    if (present(mark)) then
       if (root%marker == mark) param = .true.
    end if
    
    !Unless the current timer is the root, print the report.
    if ((.not.(root%key == "___ROOT")).and.(param)) then
       write(unit, fmt = "(4a, es13.3, es13.3, i3)") prefix, indent, root%key, space, root%wall_clock/1000, root%cpu_clock
    end if
   
    
    !Start the same recursive subroutine
    temp => root%child
    do while (associated(temp))
       call printSubTimers(indent//"  ", temp, dep-2, param)
       temp => temp%sibling
    end do

  end subroutine printSubTimers

end subroutine write_timing_tree


!!================================================
!!FLATTEN TIMING TREE
!!
!!PURPOSE(S):
!!    Print all created timers and their recorded
!!    times in an outline structure.    
!!
!!CONTAINS:
!!    printSubTimers  !Loop through all available
!!                    !timers and print times.
!!================================================
subroutine flatten_timing_tree(order, keys, walls, cpus)
  integer,          pointer :: order(:)
  character(len=*), pointer :: keys(:)
  real,             pointer :: walls(:), cpus(:)
 
  integer :: num
  integer                                          :: current
  integer                                          :: counter

  ASSERT(associated(current_timer))

  !call accumulate_elapsed_time()

  counter = 0
  current = 0

  call count(num);
  allocate(order(num*2))
  allocate(keys(num))
  allocate(walls(num))
  allocate(cpus(num))
  

  call assign_values (root_timer, num)

  contains

  recursive subroutine assign_values(root, dep)
    type(timer),                intent(in) :: root
    type(timer),       pointer             :: temp
    integer,                    intent(in) :: dep
    integer                                :: index

    if (.not.(root%key=="___ROOT")) then
       current = current + 1;
       index = current;
       counter = counter + 1
       order(counter) = index
       
       keys(index) = root%key
       walls(index) = root%wall_clock
       cpus(index) = root%cpu_clock
    end if

    temp => root%child
    do while (associated(temp))
       call assign_values(temp, dep-2)
       temp => temp%sibling
    end do

    if (.not.(root%key=="___ROOT")) then
       counter = counter +  1
       order(counter) = index
    end if
 
  end subroutine assign_values
end subroutine flatten_timing_tree


!!=================================================
!!ADD_NEW_TIMER(key)
!!
!!PURPOSE(S):
!!    Create a new timer with name [key] as a child
!!    of the current timer. If no timers exist, it
!!    will also create the root.
!!    This subroutine is private.   
!!
!!PARAMETER(S):
!!    key   !Name of the new timer
!!=================================================
subroutine add_new_timer(key)
  type(timer),       pointer             :: new_timer
  type(timer),       pointer             :: temp
  character(len=*),           intent(in) :: key

  !Create root timer if it does not exist.
  if (.not. associated(root_timer)) then
     allocate(root_timer)
     root_timer%key = "___ROOT"
     root_timer%depth = 0
     current_timer => root_timer
  end if

  !Create new timer.
  allocate(new_timer);
  new_timer%key = key
  new_timer%parent => current_timer

  if (associated(current_timer%child)) then
     temp => current_timer%child
     do while (associated(temp%sibling))
        temp => temp%sibling
     end do
     temp%sibling => new_timer
  end if

  if (.not.(associated(current_timer%child))) current_timer%child => new_timer


  current_timer => new_timer
  current_timer%depth = current_timer%parent%depth + 1

  if (maxdepth < current_timer%depth) maxdepth = current_timer%depth
  
end subroutine add_new_timer



!!================================================
!!COUNT
!!
!!PURPOSE(S):
!!    Count the number of existing timers.
!!================================================
subroutine count(num)
  integer, intent(out) :: num
  num = 0
  call printSubTimers (root_timer, 0)
  num = num-1
contains
  recursive subroutine printSubTimers(root, mode)
    type(timer),        intent(in)  :: root
    integer,            intent(in)  :: mode
    type(timer),        pointer     :: temp
    num = num +1;
    temp => root%child
    do while (associated(temp))
       call printSubTimers(temp, mode)
       temp => temp%sibling
    end do
  end subroutine printSubTimers
end subroutine count


!!================================================
!!ACCUMULATE ELAPSED TIME
!!
!!PURPOSE(S):
!!    Stop current timer and all parent timers if running.
!!================================================
subroutine accumulate_elapsed_time()
  type(timer), pointer :: temp
  integer              :: wall_buffer
  real                 :: cpu_buffer

  temp => current_timer

  do while (associated(temp))
     call system_clock(wall_buffer)
     call cpu_time(cpu_buffer)
     temp%wall_clock = temp%wall_clock + (wall_buffer - temp%wall_start)
     temp%cpu_clock = temp%cpu_clock + (cpu_buffer - temp%cpu_start)
     temp => temp%parent
  end do

end subroutine accumulate_elapsed_time


end module timing_tree
