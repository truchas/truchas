MODULE PGSlib_Stats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines for instrumenting code
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_stats.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $

!  USE PGSLib_Misc_Utility
  use pgslib_error_module
  USE PGSLib_Type_MODULE
  USE PGSLib_Timing_MODULE

  PRIVATE
  SAVE
  PUBLIC :: Enter_Routine, Exit_Routine
  PUBLIC :: Initialize_Instrument_Array

  ! Routines for pointing to routine names
  PUBLIC :: PGSLib_Instrument_T
  PUBLIC :: SETUP_TRACE_STATISTICS
  PUBLIC :: SETUP_Supplement_STATISTICS
  PUBLIC :: Supplement_G_TO_L_Statistics
  PUBLIC :: Supplement_L_From_G_Statistics
  PUBLIC :: Add_Item_To_T_STATISTICS
  PUBLIC :: Item_Index_From_T_STATISTICS
  PUBLIC :: GS_INIT_TRACE_STATISTICS
  PUBLIC :: GS_RELEASE_TRACE_STATISTICS
  PUBLIC :: SETUP_BUFFERS_STATISTICS
  PUBLIC :: GATHER_BUFFER_STATISTICS
  PUBLIC :: SCATTER_BUFFER_STATISTICS
  PUBLIC :: GATHER_STATISTICS
  PUBLIC :: SCATTER_STATISTICS
  PUBLIC :: GLOBAL_ALL_STATISTICS
  PUBLIC :: GLOBAL_ANY_STATISTICS
  PUBLIC :: GLOBAL_DOT_PRODUCT_STATISTICS
  PUBLIC :: GLOBAL_MAX_STATISTICS
  PUBLIC :: GLOBAL_MIN_STATISTICS
  PUBLIC :: GLOBAL_SUM_STATISTICS

!!!!! Array holds the instrumentation for each routine
  integer, parameter :: Max_Instrument_Slots = 256
  type (PGSLib_Instrument_T), TARGET, dimension(Max_Instrument_Slots) :: Instrument_Array

CONTAINS
  subroutine Enter_Routine(Slot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !           Start timers, perform other instrumentation init
    !           required for start of a routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type (PGSLib_Instrument_T), POINTER :: Slot

    ! Stop stack timer from routine above
    call PGSLib_Stop_Timer()
    ! Start new stack timer
    call PGSLib_Push_Timer()
    ! Start routine timer
    call pgslib_start_timer(Slot%Timer)

  end subroutine Enter_Routine
  
  subroutine Exit_Routine(Slot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !           Stop timers, perform other instrumentation clean-up
    !           required for exit from a routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type (PGSLib_Instrument_T), POINTER :: Slot

    ! Stop the routine timer
    call pgslib_stop_timer(Slot%Timer)
    ! Stop stack timer
    Call PGSLib_Stop_Timer()
    ! Read stack timer, also pops stack.  This is time in this routine only, excluding
    ! routines below it (assuming they have been instrumented)
    Slot%Timer%Internal_Time = PGSLib_Read_Elapsed_Time()
    ! Start stack timer for routine above
    Call PGSLib_Start_Timer()

  end subroutine Exit_Routine
  
  subroutine Initialize_Instrument_Array()
    implicit none
    integer i
    do i = 1, Max_Instrument_Slots
       call pgslib_reset_timer(Instrument_Array(i)%Timer)
    end do
    return
  end subroutine Initialize_Instrument_Array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This section contains the mapping of Routine references to 
    ! instrumentation Slot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Next_Instrument_Item()
  implicit none
  type (PGSLib_Instrument_T), POINTER :: Next_Instrument_Item
  type (PGSLib_Instrument_slot), SAVE :: Next_Instrument_Item_Slot = PGSLib_Instrument_Slot(0)
  ! Increment timer slot, make sure it is not too large, if not return value
  Next_Instrument_Item_Slot%Slot = Next_Instrument_Item_Slot%Slot + 1
  if (Next_Instrument_Item_Slot%Slot > Max_Instrument_Slots) then
!     call pgslib_fatal_error("Using too many timers in Next_Instrument_Item")
  end if
  Next_Instrument_Item => Instrument_Array(Next_Instrument_Item_Slot%Slot)
  Next_Instrument_Item%Slot%Slot = Next_Instrument_Item_Slot%Slot
  return
end function Next_Instrument_Item

  
#define _FUNC_ SETUP_TRACE
#include "stats_name.fpp"

#define _FUNC_ SETUP_SUPPLEMENT
#include "stats_name.fpp"

#define _FUNC_ SUPPLEMENT_G_TO_L
#include "stats_name.fpp"

#define _FUNC_ SUPPLEMENT_L_FROM_G
#include "stats_name.fpp"

#define _FUNC_ Add_Item_To_T
#include "stats_name.fpp"

#define _FUNC_ Item_Index_From_T
#include "stats_name.fpp"

#define _FUNC_ GS_INIT_TRACE
#include "stats_name.fpp"

#define _FUNC_ GS_RELEASE_TRACE
#include "stats_name.fpp"

#define _FUNC_ SETUP_BUFFERS
#include "stats_name.fpp"

#define _FUNC_ GATHER_BUFFER
#include "stats_name.fpp"

#define _FUNC_ SCATTER_BUFFER
#include "stats_name.fpp"

#define _FUNC_ GATHER
#include "stats_name.fpp"

#define _FUNC_ SCATTER
#include "stats_name.fpp"

#define _FUNC_ GLOBAL_ALL
#include "stats_name.fpp"

#define _FUNC_ GLOBAL_ANY
#include "stats_name.fpp"

#define _FUNC_ GLOBAL_DOT_PRODUCT
#include "stats_name.fpp"

#define _FUNC_ GLOBAL_MAX
#include "stats_name.fpp"

#define _FUNC_ GLOBAL_MIN
#include "stats_name.fpp"

#define _FUNC_ GLOBAL_SUM
#include "stats_name.fpp"



END MODULE PGSlib_Stats
    


