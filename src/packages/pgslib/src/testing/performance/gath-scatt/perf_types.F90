! The data and types used for performance measures
! $Id: perf_types.F,v 1.1.1.1 2000/10/11 22:44:25 ferrell Exp $
MODULE perf_types
  use pgslib_module
  implicit none

  PRIVATE
  PUBLIC :: GS_Perf_Timer_Type
  PUBLIC :: TotalTime, OffPETime, SetupTime
  PUBLIC :: NetBW, NetBWNorm, offPEBW, offPEBWNorm
  PUBLIC :: TotalData, offPEData, BytesPerWord
  PUBLIC :: GatherDegree, ScatterDegree

  type GS_Perf_Timer_Type
     real :: SetupTime
     real :: TotalTime, OffPETime
     real :: NetBW, NetBWNorm, offPEBW, offPEBWNorm
     integer :: TotalData, offPEData, BytesPerWord     
     integer :: GatherDegree, ScatterDegree
  end type GS_Perf_Timer_Type

CONTAINS
  function SetupTime(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: SetupTime
    SetupTime = T%SetupTime
    return
  end function SetupTime

  function TotalTime(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: TotalTime
    TotalTime = T%TotalTime
    return
  end function TotalTime

  function OffPETime(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: OffPETime
    OffPETime = T%OffPETime
    return
  end function OffPETime

  function NetBW(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: NetBW
    NetBW = T%NetBW
    return
  end function NetBW

  function NetBWNorm(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: NetBWNorm
    NetBWNorm = T%NetBWNorm
    return
  end function NetBWNorm

  function offPEBW(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: offPEBW
    offPEBW = T%offPEBW
    return
  end function offPEBW

  function offPEBWNorm(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: offPEBWNorm
    offPEBWNorm = T%offPEBWNorm
    return
  end function offPEBWNorm

  function TotalData(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: TotalData
    TotalData = T%TotalData
    return
  end function TotalData

  function offPEData(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    real :: offPEData
    offPEData = T%offPEData
    return
  end function offPEData

  function BytesPerWord(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    integer :: BytesPerWord
    BytesPerWord = T%BytesPerWord
    return
  end function BytesPerWord

  function GatherDegree(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    integer :: GatherDegree
    GatherDegree = T%GatherDegree
    return
  end function GatherDegree

  function ScatterDegree(T)
    type (GS_Perf_Timer_Type), intent(IN) :: T
    integer :: ScatterDegree
    ScatterDegree = T%ScatterDegree
    return
  end function ScatterDegree

end MODULE perf_types

  
