MODULE PGSLib_Timing_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Support for timing/instrumenting codes
  !
  ! This module provides the PGSLib_Timer.  Timers should be declared
  ! with this type.  Routines in this module which take a timer
  ! as an argument expect the argument to be of type PGSLib_Timer.
  !
  ! Most routines can be called with or without a timer.  If a
  ! timer is not passed, then an internal timer stack is used.  The
  ! behavior depends on the routine.
  ! In general, if there is an active current stack timer, then it
  ! is used.  Otherwise a new timer is pushed onto the stack.
  !    

  ! Routines in this module:
  !    T = argument of PGSLib_Timer
  !   *T = optional argument 
  !
  !    PGSLib_Reset_Timer(*T)        - clear the elapsed time, or the current stack
  !                                    timer for T not present.
  !    PGSLib_Start_Timer(*T)        - start the timer.  If T is not present, then
  !                                    allocate a new stack timer if one is not running.
  !                                    If the current_timer is running, then get a new 
  !                                    stack timer.
  !    PGSLib_Stop_Timer(*T)         - Stop timer, accumulate elapsed time.  If
  !                                    T is not present, stop the current_timer, if any.
  !PGSLib_Print_Elapsed_Time(*T, *S) - Print the elapsed time.  If string S is present,
  !                                    print that string ahead of the elapsed time.
  !                                    If T is not present, then print the current_timer,
  !                                    if any.  If current_timer is not running, pop
  !                                    it off the stack, otherwise leave it running.
  !     PGSLib_Read_Elapsed_Time(*T) - Return the elapsed time, with the same
  !                                    behavior as PGSLib_Print_Elapsed_Time.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_timing_module.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $
  use pgslib_type_module

  implicit none
  save
  private
  public :: PGSLib_Start_Timer, PGSLib_Stop_Timer
  public :: PGSLib_Increment_Timer
  public :: PGSLib_Push_Timer
  public :: PGSLib_Read_Elapsed_Time, PGSLib_Print_Elapsed_Time
  public :: PGSLib_Reset_Timer
  public :: PGSLib_Timer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Generic procedures in this module
  
  INTERFACE PGSLib_Read_Elapsed_Time
     module procedure PGSLib_Read_Elapsed_Time_Timer
  END INTERFACE


  type Timer_Stack_Item
     type (PGSLib_Timer),       Pointer :: Timer
     type (Timer_Stack_Item), Pointer :: Previous
  end type Timer_Stack_Item

  integer :: Timers_Allocated = 0
  type (Timer_Stack_Item), POINTER :: Current_Timer 

  logical :: initialized = .FALSE.

! This is the proper definition of the parameters
  integer, parameter :: MSec = 8, Sec = 7, Min = 6, Hour = 5, Day = 3, Mnth = 2, Year = 1

! This is the definition that must be used on the DEC, with multi-proc HPF
!  integer, parameter :: MSec = 9, Sec = 8, Min = 7, Hour = 6, Day = 4, Mnth = 3, Year = 2

  integer, parameter :: MIN_PRINT_LENGTH = 35

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Timer list.
! This list is all the global timers used in the code.  


CONTAINS
  subroutine New_Timer_From_Stack()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate a new timer, push it onto the stack.

    implicit none
    
    ! Local variables
    integer :: memerror
    type (Timer_Stack_Item), POINTER  :: New_Timer_Item
    
    if (.not. initialized) then
       initialized = .true.
       NULLIFY(Current_Timer)
    end if


    ! Allocate a new timer stack item, and its associated timer
    ALLOCATE(New_Timer_Item, stat=memerror)
    if (memerror /= 0) then
       NULLIFY(New_Timer_Item)
       print *, 'WARNING: Could not allocate a new timer.'
       RETURN
    endif
    ALLOCATE(New_Timer_Item%Timer, stat=memerror)
    if (memerror /= 0) then
       NULLIFY(New_Timer_Item%Timer)
       print *, 'WARNING: Could not allocate a new timer.'
       RETURN
    endif

    ! Initialize the new timer.
    Timers_Allocated = Timers_Allocated + 1
    New_Timer_Item%Timer%Timer_N = Timers_Allocated
    New_Timer_Item%Timer%Running = .FALSE.
    New_Timer_Item%Timer%Start_Time  = -HUGE(0)
    New_Timer_Item%Timer%Stop_Time   = -HUGE(0)
    call PGSLib_Reset_Timer(New_Timer_Item%Timer)
    New_Timer_Item%Timer%Start_Clock  = -HUGE(0)
    New_Timer_Item%Timer%Stop_Clock   = -HUGE(0)
    New_Timer_Item%Timer%Timer_String = ""

    ! Point Current_Timer to the new stack item

    if (.NOT. ASSOCIATED(Current_Timer)) then
       Current_Timer => New_Timer_Item
       NULLIFY(Current_Timer%Previous)
    else
       New_Timer_Item%Previous => Current_Timer
       Current_Timer => New_Timer_Item
    end if

    RETURN
  END subroutine New_Timer_From_Stack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Decrement_Timer_Stack()
    implicit none

    ! Local Variables
    type (Timer_Stack_Item), POINTER :: Item_To_Pop

    ! If there is no current stack item, then no stack to decrement
    if (.NOT. ASSOCIATED(Current_Timer))  RETURN
    
    ! If the current timer is pointing to a timer, then deallocate it
    if (ASSOCIATED(Current_Timer%Timer)) then
       DEALLOCATE(Current_Timer%Timer)
       NULLIFY(Current_Timer%Timer)
    endif

    ! If the current timer is not the end of the stack, reset
    ! current timer to point to previous item
    ! If there is no previous, the nullify Current_Timer.
    if (ASSOCIATED(Current_Timer%Previous)) then
       Item_To_Pop => Current_Timer
       Current_Timer => Current_Timer%Previous
       DEALLOCATE(Item_To_Pop)
       NULLIFY(Item_To_Pop)
       Timers_Allocated = Timers_Allocated - 1
    else
       DEALLOCATE(Current_Timer)
       NULLIFY(Current_Timer)
       Timers_Allocated = Timers_Allocated - 1
    end if

    RETURN
  END SUBROUTINE Decrement_Timer_Stack
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Start_Timer(T, String)
    implicit none
    type (PGSLib_Timer), OPTIONAL, TARGET :: T
    character (LEN=*), OPTIONAL, TARGET :: String


    ! If Timer argument T is present, then start timer T.  
    ! If Timer argument T is not present, the start a stack timer.  If
    ! the current stack timer is not running, then start it.  If the
    ! current stack timer is running, then push a new timer onto the stack
    ! and start it.
    ! In all cases, only set new start times if the ref_counter
    ! for the timer == 0.  In other cases do nothing.
    ! In all cases increment the ref counter.

    ! Local Variables
    type (PGSLib_Timer), POINTER :: Local_T

    Local_T => Get_Timer(T)

    ! If we are using a stack timer, may need a new one.
    if (.NOT. PRESENT(T)) then
       ! If we do not have a current timer, get a new one from the stack.
       if (.NOT. ASSOCIATED(Local_T)) then
          Call New_Timer_From_Stack()       
          Local_T => Current_Timer%Timer
       else
          ! If the current stack timer is running, then get a new one from the stack
          if (ASSOCIATED(Current_Timer)) then
             if (Current_Timer%Timer%Running) then
                Call New_Timer_From_Stack()
                Local_T => Current_Timer%Timer
             end if
          endif
       end if
    end if
       
    if (Local_T%T_Ref_Counter == 0) then
       call DATE_AND_TIME(VALUES=Local_T%Start_Clock)
       call cpu_time (Local_T%Start_Time)
    end if
    ! Increment timer counter
    Local_T%T_Ref_Counter = Local_T%T_Ref_Counter + 1
    Local_T%Running = .TRUE.

    if (PRESENT(String)) then
       Local_T%Timer_String = String
    end if

    RETURN
  END SUBROUTINE PGSLib_Start_Timer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Push_Timer(String)
    implicit none
    character (LEN=*), OPTIONAL, TARGET :: String

    ! Push a new timer onto the stack and start it.

    ! New timer
    Call New_Timer_From_Stack()

    ! Start it.
    call PGSLib_Start_Timer(STRING=String)

    RETURN
  END SUBROUTINE PGSLib_Push_Timer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Stop_Timer(T)
    implicit none
    type (PGSLib_Timer), OPTIONAL, TARGET :: T

    ! Stop timer T.  If no timer is given, then stop current stack
    ! timer.  If there is no current stack timer, then do nothing.

    ! Local Variables
    type (PGSLib_Timer), POINTER :: Local_T
    
    Local_T => Get_Timer(T)

    ! If there is no timer to stop, do nothing.
    if (.NOT. ASSOCIATED(Local_T)) then
       RETURN
    END if

    ! Decrement the reference counter
    Local_T%T_Ref_Counter = Local_T%T_Ref_Counter - 1

    ! If this is the last reference, find the stop time, and stop it

    if (Local_T%T_Ref_Counter == 0) then
       if (Local_T%Running) then
          call DATE_AND_TIME(VALUES=Local_T%Stop_Clock)
          Local_T%Elapsed_Clock = Local_T%Elapsed_Clock + Compute_Elapsed_Clock(Local_T)
          call cpu_time (Local_T%Stop_Time)
          Local_T%Elapsed_Time = Local_T%Elapsed_Time + (Local_T%Stop_Time - Local_T%Start_Time)
       endif
       Local_T%Running = .FALSE.
    end if

    RETURN
  END SUBROUTINE PGSLib_Stop_Timer
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Increment_Timer(T, Increment)
    implicit none
    type (PGSLib_Timer) :: T
    real                :: Increment

    ! Increment the Elapsed_Time field of the Timer

    ! Local Variables

    T%Elapsed_Time = T%Elapsed_Time + Increment

    RETURN
  END SUBROUTINE PGSLib_Increment_Timer
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Reset_Timer(T, String)
    implicit none
    type (PGSLib_Timer), OPTIONAL, TARGET :: T
    character (LEN=*), OPTIONAL, TARGET :: String

    ! Reset timer T (elapsed_time = 0.0).  If no timer is given, then stop current stack
    ! timer.  If there is no current stack timer, then do nothing.

    ! Local Variables
    type (PGSLib_Timer), POINTER :: Local_T
    
    Local_T => Get_Timer(T)

    ! If there is no timer to stop, do nothing.
    if (.NOT. ASSOCIATED(Local_T)) RETURN

    ! If there is a timer, reset it
    Local_T%Elapsed_Clock = 0.0
    Local_T%Elapsed_Time = 0.0
    Local_T%Internal_Time = 0.0
    Local_T%Running = .FALSE.
    Local_T%T_Ref_Counter = 0

    ! If String is present, set it
    if (PRESENT(String)) then
       Local_T%Timer_String = String
    else
       Local_T%Timer_String = ""
    end if

    
    RETURN
  END SUBROUTINE PGSLib_Reset_Timer
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function PGSLib_Read_Elapsed_Time_Timer(T)
    implicit none
    real PGSLib_Read_Elapsed_Time_Timer
    type (PGSLib_Timer), intent(IN), OPTIONAL, TARGET :: T

    ! Return the elapsed time, in seconds, as a real number
    ! If timer T is passed, then read elapsed time for that timer.
    ! If no timer is passed, then the current stack timer, Current_Timer%Timer
    ! is read.
    ! In either case, if the timer is running, it remains running
    ! If Current_Timer%Timer is read, and if it is stopped, then
    ! Current_Timer is popped off the time stack.
    ! If T not present, and there is no stack timer, do nothing,
    ! and return 0.0.

    ! Local Variables
    type (PGSLib_Timer), POINTER :: Local_T

    Local_T => Get_Timer(T)

    ! If there is no timer to read, then return 0.0.
    if (.NOT. ASSOCIATED(Local_T)) then
       PGSLib_Read_Elapsed_Time_Timer = 0.0
       RETURN
    end if

    ! If the timer is running, then we have to read the current time.
    if (Local_T%Running) then
       call DATE_AND_TIME(VALUES=Local_T%Stop_Clock)
       call cpu_time (Local_T%Stop_Time)
       PGSLib_Read_Elapsed_Time_Timer = Local_T%Elapsed_Time + (Local_T%Stop_Time - Local_T%Start_Time)
    else
       PGSLib_Read_Elapsed_Time_Timer = Local_T%Elapsed_Time
    end if
       
    ! If reading the stack timer, and it is stopped, then pop it off the stack
    if (.NOT. PRESENT(T)) then
       if (.NOT. Local_T%Running) then
          Call Decrement_Timer_Stack()
       end if
    end if

    
    RETURN
  END function PGSLib_Read_Elapsed_Time_Timer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PGSLib_Print_Elapsed_Time(T, STRING)
    implicit none
    type (PGSLib_Timer), intent(IN),         optional :: T
    character (LEN=*), intent(IN), TARGET, optional :: STRING

    ! Print the elapsed time for timer T, if passed in, or for
    ! the current stack pointer.  Uses PGSLib_Read_Elapsed_Time, so
    ! refer to that routine for details.
    ! Uses the Timer_String field for printing unless
    ! an optional string is passed in, in which case that string
    ! is used.

    ! Local variables
    real :: elapsed
    integer :: length, string_length
    character (LEN=1024) :: local_string
    type (PGSLib_Timer), POINTER :: Local_T

    Local_T => Get_Timer(T)

    if (.NOT. ASSOCIATED(Local_T)) THEN
       RETURN
    END if

    IF (PRESENT(STRING) ) then
       local_string = TRIM(STRING)
    else
       local_string = TRIM(Local_T%Timer_String)
    end IF

    elapsed = PGSLib_Read_Elapsed_Time(Local_T)
    Length = LEN_TRIM(local_string)
    string_length = MAX(MIN_PRINT_LENGTH,Length)
    
    WRITE(*, 39) local_string(1:Length), REPEAT(" ",NCOPIES=MAX(0,MIN_PRINT_LENGTH - Length)), elapsed
39  FORMAT(1X,a,a,f11.3)

    RETURN
  END SUBROUTINE PGSLib_Print_Elapsed_Time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function Get_Timer(T) RESULT(Timer_Ptr)
    implicit none
    type (PGSLib_Timer),                      POINTER :: Timer_Ptr
    type (PGSLib_Timer), intent(IN), OPTIONAL, TARGET :: T

    ! Return a pointer to a timer.
    ! If T is present, then return a pointer to T.
    ! If T is not present, then return a pointer to the Current_Timer%Timer,
    ! if that is associated.  If not, then return a null pointer.

    ! Local Variables

    if (PRESENT(T)) then
       Timer_Ptr => T
    else
       IF (ASSOCIATED(Current_Timer)) then
          if (ASSOCIATED(Current_Timer%Timer)) then
             Timer_Ptr => Current_Timer%Timer
          else
             NULLIFY(Timer_Ptr)
          endif
       else
          NULLIFY(Timer_Ptr)
       end IF
    end if
    RETURN
  end function Get_Timer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function Compute_Elapsed_Clock(Local_T)
    ! Returns elapsed time (difference between stop_clock and start_clock)
    ! in units of seconds.

    implicit none
    real                          :: Compute_Elapsed_Clock
    type (PGSLib_Timer), intent(IN) :: Local_T

    ! Local variables
    real MilliSeconds

    MilliSeconds = (Local_T%Stop_Clock(MSec) - Local_T%Start_Clock(MSec))  +  & ! Milli Seconds
         &  1000.*( (Local_T%Stop_Clock(Sec)  - Local_T%Start_Clock(Sec))   +  & ! Seconds
         &    60.*( (Local_T%Stop_Clock(Min)  - Local_T%Start_Clock(Min))   +  & ! Minutes
         &    60.*( (Local_T%Stop_Clock(Hour) - Local_T%Start_Clock(Hour))  +  & ! Hours
         &    24.*( (Local_T%Stop_Clock(Day)  - Local_T%Start_Clock(Day))   +  & ! Days
         &    12.*( (Local_T%Stop_Clock(Mnth) - Local_T%Start_Clock(Mnth))  +  & ! Months
         &   365.*( (Local_T%Stop_Clock(Year) - Local_T%Start_Clock(Year))     & ! Years
         &        ))))))

    Compute_Elapsed_Clock = MilliSeconds/1000.0
    RETURN
  END function Compute_Elapsed_Clock

     


END MODULE PGSLib_Timing_MODULE

