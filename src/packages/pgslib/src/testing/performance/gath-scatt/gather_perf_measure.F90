! Measure performance of the pgslib gather routines
! $Id: gather_perf_measure.F,v 1.1.1.1 2000/10/11 22:44:25 ferrell Exp $
MODULE gather_perf_measure
  use pgslib_module
  use perf_types

  implicit none

  PRIVATE
  PUBLIC :: gather_performance

CONTAINS

  function gather_performance(INDEX, NDest, NLoop) RESULT(Timings)
    !======================================================================
    !          gather_performance(Index, NLoop)
    ! Measure performance of gather routine.  Index is input pattern
    ! NDest tells how big an array we are gathering from.
    ! NLoop is number of times to performan gather.
    ! Output is a GS_Perf_Timer_Type, which contains the timing info
    !======================================================================
    

    integer(PGSLib_Int_Type),               &
                              dimension(:), &
                              intent(IN)    :: Index
    integer(PGSLib_Int_Type), intent(IN)    :: NDest
    integer(PGSLib_Int_Type), intent(IN)    :: NLoop
    type(GS_Perf_Timer_Type)                :: Timings


    ! Local variables
    integer                                 :: l
    type (PGSLib_GS_Trace), POINTER         :: Trace
    integer(PGSLib_Int_Type),               &
                   dimension(SIZE(Index,1)) :: Local_Index
!!$    integer (PGSLib_Int_Type),              &
!!$                   dimension(SIZE(Index,1)) :: Dest
!!$    integer (PGSLib_Int_Type),              &
!!$                   dimension(NDest)         :: Source
    real (PGSLib_Single_Type),              &
                   dimension(SIZE(Index,1)) :: Dest
    real (PGSLib_Single_Type),              &
                   dimension(NDest)         :: Source
    
    integer (PGSLib_Int_Type)               :: TotalWords, TotalOffPEWords
    
    ! Setup the trace
    Local_Index = Index
    Trace => PGSLib_Setup_Trace(Local_Index, SIZE_OF_DEST=NDest) ! Got this backward because thinking only of gather

    ! Setup the source
    call Random_Number(HARVEST = Source)

    ! Initialize the dest
    Dest = 0
    ! Loop the proper number of times
    do l = 1, NLoop
       call pgslib_gather(Dest, Source, Local_Index, Trace)
    end do


    ! Grab performance data from trace.
    call pgslib_trace_degree(Timings%ScatterDegree, Timings%GatherDegree, Trace)

    Timings%SetupTime = PGSLib_Read_Maximum_Time(Trace%Setup_Timer)
    ! Each iterartion we gather a word for every element in Index.
    Timings%TotalData = NLoop * PGSLib_Global_SUM(Size(Local_Index,1))
    
    Timings%TotalTime = PGSLib_Read_Maximum_Time(Trace%GatherTotal_Timer)
    Timings%NetBW     = Timings%TotalData/Timings%TotalTime
    Timings%NetBWNorm = Timings%NetBW/PGSLib_Inquire_nPE()
    
    ! Each iteration, every item in supplement buffer gets moved between processors
    Timings%offPEData   = NLoop * PGSLib_Global_SUM(PGSLib_Size_Of_Sup(Trace))

    Timings%OffPETime   = PGSLib_Read_Maximum_Time(Trace%GatherBuffer_Timer)
    if (Timings%OffPETime > 0) then
       Timings%OffPEBW     = Timings%offPEData/Timings%OffPETime
    else
       Timings%OffPEBW = 0.0
    end if
    Timings%OffPEBWNorm = Timings%OffPEBW/PGSLib_Inquire_nPE()

    ! We are finished with the trace
    call PGSLib_Deallocate_Trace(Trace)

    ! We are completely done, we have got all the numbers we need, so go home
    RETURN
  end function gather_performance

end MODULE gather_perf_measure
