MODULE PGSlib_Instrument
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines for instrumenting code
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_instrument.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $

  USE PGSLib_Misc_Utility
  use pgslib_error_module
  use pgslib_red_numeric_module
  use pgslib_stats
  USE PGSLib_Type_MODULE
  USE PGSLib_Timing_MODULE

  PRIVATE
  SAVE
  PUBLIC :: PGSLib_Read_Maximum_Time
  PUBLIC :: PGSLib_Read_Elapsed_Time

  ! Routines for getting stats from C stuff
  PUBLIC :: PGSLib_Trace_Degree
  PUBLIC :: PGSLib_Barrier_Time
  PUBLIC :: PGSLib_SR_Time

  ! Interfaces for Generic routines
  INTERFACE PGSLib_Read_Maximum_Time
     MODULE PROCEDURE PGSLib_Read_Maximum_Time_Timer
     MODULE PROCEDURE PGSLib_Read_Maximum_Time_Array
  END INTERFACE

  INTERFACE PGSLib_Read_Elapsed_Time
     MODULE PROCEDURE PGSLib_Read_Elapsed_Time_Array
  END INTERFACE

  interface PGSLib_Trace_Degree
     module procedure PGSLib_Trace_Degree_F
  end interface
  
  ! Interfaces to C routines
  interface PGSLib_Trace_Degree_C
     subroutine PGSLib_Trace_Degree_C(Scatter_Degree, Gather_Degree, GS_Trace)
       implicit none
       integer :: Scatter_Degree, Gather_Degree
       integer, DIMENSION(*) :: GS_Trace
     end subroutine PGSLib_Trace_Degree_C
  end interface

  interface PGSLib_Barrier_Time_C
     subroutine PGSLib_Barrier_Time_C(BT)
       USE PGSLib_Type_MODULE
       implicit none
       real (PGSLib_Real_Type) :: BT
     end subroutine PGSLib_Barrier_Time_C
  end interface

  interface PGSLib_Sr_Time_C
     subroutine PGSLib_Sr_Time_C(T)
       USE PGSLib_Type_MODULE
       implicit none
       real (PGSLib_Real_Type) :: T
     end subroutine PGSLib_Sr_Time_C
  end interface

CONTAINS
  function PGSLib_Read_Maximum_Time_Timer(Timer) RESULT(Max_Time)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !          Report the Maximum of all the times of a timer.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type (PGSLib_Timer)     :: Timer
    real (PGSLib_REAL_TYPE) :: MAx_Time

    Max_Time = PGSLib_Global_MAXVAL( PGSLib_Read_Elapsed_Time(Timer))
    return
  end function PGSLib_Read_Maximum_Time_Timer
  

  function PGSLib_Read_Maximum_Time_Array(Slot) RESULT(Max_Time)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !          Report the Maximum of all the times of Slot.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type (PGSLib_Instrument_T), POINTER :: Slot
    real (PGSLib_REAL_TYPE) :: Max_Time

    Max_Time = PGSLib_Global_MAXVAL( PGSLib_Read_Elapsed_Time(Slot%Timer) )
    return
  end function PGSLib_Read_Maximum_Time_Array
  

  function PGSLib_Read_Elapsed_Time_Array(Slot) RESULT(Elapsed_Time)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !          Report the Elapsed of all the times of Slot.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type (PGSLib_Instrument_T), POINTER :: Slot
    real (PGSLib_REAL_TYPE) :: Elapsed_Time

    Elapsed_Time = PGSLib_Read_Elapsed_Time(Slot%Timer) 
    return
  end function PGSLib_Read_Elapsed_Time_Array
  

  function PGSLib_Barrier_Time() RESULT(Max_Time)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !          Report the maximum time used by barriers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    real (PGSLib_Single_TYPE) :: Max_Time
    real (PGSLib_Single_Type) :: bt

    call pgslib_barrier_time_c(bt)

    Max_Time = PGSLib_Global_MAXVAL( bt )
    return
  end function PGSLib_Barrier_Time
  


  function PGSLib_Sr_Time() RESULT(Max_Time)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PURPOSE :
    !          Report maximum time used in send-receive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    real (PGSLib_Single_TYPE) :: Max_Time
    real (PGSLib_Single_Type) :: bt

    call pgslib_Sr_time_c(bt)

    Max_Time = PGSLib_Global_MAXVAL( bt )
    return
  end function PGSLib_Sr_Time
  



  subroutine PGSLib_Trace_Degree_F(Scatter_Degree, Gather_Degree, Trace)
    use PGSLib_Type_MODULE
    implicit none
    integer :: Scatter_Degree, Gather_Degree
    type (PGSLib_GS_Trace) :: Trace
    call pgslib_trace_degree_c(Scatter_Degree, Gather_Degree, Trace%GS_Trace)
    return
  end subroutine PGSLib_Trace_Degree_F


END MODULE PGSlib_Instrument
    


