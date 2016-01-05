!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_CHARTS_Atlases
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition joint
  !   chart and atlas operations in Telluride.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  use bc_enum_types
  use bc_atlases, only : BC_Atlas, BC_Get_Values, BC_Get_ValueIndex, BC_Get_UseFunction, &
                         BC_Get_Positions, &
                         BC_Get_Face, BC_Get_Cell,&
                         BC_Get_Offset, BC_Get_Length, BC_Get_DOF, &
                         DIMENSIONALITY, BC_NumberOfCharts, APPEND,&
                         INITIALIZE, FREE
  use bc_charts

  Implicit None
  Private

  PUBLIC :: BC_Atlas, FREE
  PUBLIC :: BC_Chart_ID
  PUBLIC :: INITIALIZE, APPEND
  PUBLIC :: BC_Invalid_Chart_ID
  PUBLIC :: BC_Chart_ID_Is_Valid
  PUBLIC :: BC_Chart_Matches
  PUBLIC :: BC_Chart_Set_Chart
  PUBLIC :: BC_Chart_Next_Chart
  
  PUBLIC :: BC_Chart_Length
  PUBLIC :: BC_Chart_Values
  PUBLIC :: BC_Chart_Positions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! A chart ID, identifies a chart in an atlas
  type BC_CHART_ID
     PRIVATE
     integer :: AtlasIndex
     type(BC_Atlas), POINTER :: Atlas => NULL()
  end type BC_CHART_ID

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE InitChartID
  END INTERFACE

  INTERFACE APPEND
     MODULE PROCEDURE AppendChartAtlas
  END INTERFACE

  INTERFACE BC_Invalid_Chart_ID
     MODULE PROCEDURE InvalidChart_ID
  END INTERFACE

  INTERFACE BC_Chart_ID_Is_Valid
     MODULE PROCEDURE ChartIDIsValid
  END INTERFACE
     
  INTERFACE BC_CHART_Atlas
     MODULE PROCEDURE ChartIDAtlas
  END INTERFACE
  
  INTERFACE BC_CHART_AtlasIndex
     MODULE PROCEDURE ChartIDAtlasIndex
  END INTERFACE
  
  INTERFACE BC_CHART_LENGTH
     MODULE PROCEDURE ChartIDLength
  END INTERFACE
  
  INTERFACE BC_CHART_VALUES
     MODULE PROCEDURE GetChartIDValues
  END INTERFACE

  INTERFACE BC_Chart_Positions
     MODULE PROCEDURE GetChartIDPositions
  END INTERFACE

  INTERFACE BC_CHART_CELL
     MODULE PROCEDURE ChartIDCell
  END INTERFACE

  INTERFACE BC_CHART_FACE
     MODULE PROCEDURE ChartIDFace
  END INTERFACE
  
  INTERFACE BC_Chart_Matches
     MODULE PROCEDURE ChartIDMatches
  END INTERFACE

  INTERFACE BC_Chart_Set_Chart
     MODULE PROCEDURE SetChartID
  END INTERFACE

  INTERFACE BC_CHART_NEXT_CHART
     MODULE PROCEDURE SetNextChartID
  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate a Chart_ID
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function InvalidChart_ID() RESULT(ChartID)
    ! Return a pointer to an invalid chart ID
    type(BC_Chart_ID), TARGET, SAVE :: Invalid_Chart_ID
    logical, SAVE :: Initialized = .FALSE.
    type(BC_Chart_ID), POINTER :: ChartID

    if (.NOT. Initialized) then
       ! No valid Atlas
       NULLIFY(Invalid_Chart_ID%Atlas)
       ! No valid Atlas pointer
       Invalid_Chart_ID%AtlasIndex = BC_INVALID_ATLASINDEX
       Initialized = .TRUE.
    end if
    ChartID => Invalid_Chart_ID
  end function InvalidChart_ID

  subroutine InitChartID(ChartID)
    type(BC_Chart_ID), target, intent(OUT) :: ChartID
    ChartID%AtlasIndex = BC_INVALID_ATLASINDEX
    NULLIFY(ChartID%Atlas)
  end subroutine InitChartID

  function ChartIDIsValid(ChartID)
    ! Returns .TRUE. if ChartID is a valid chart ID.
    ! Otherwise returns .FALSE.
    ! The test is on the AtlasIndex field.
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    logical                        :: ChartIDIsValid
    ChartIDIsValid = ChartID%AtlasIndex /= BC_INVALID_ATLASINDEX
  end function ChartIDIsValid

  function ChartIDAtlas(ChartID)
    type(BC_Chart_ID), TARGET, intent(IN):: ChartID
    type(BC_Atlas), POINTER :: ChartIDAtlas
    ChartIDAtlas => ChartID%Atlas
  end function ChartIDAtlas

  function ChartIDAtlasIndex(ChartID)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    integer :: ChartIDAtlasIndex
    ChartIDAtlasIndex = ChartID%AtlasIndex
  end function ChartIDAtlasIndex

  function ChartIDCell(ChartID)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    integer :: ChartIDCell
    ChartIDCell = BC_Get_Cell(BC_Chart_Atlas(ChartID), BC_CHART_AtlasIndex(ChartID))
  end function ChartIDCell

  function ChartIDFace(ChartID)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    integer :: ChartIDFace
    ChartIDFace = BC_Get_Face(BC_Chart_Atlas(ChartID), BC_CHART_AtlasIndex(ChartID))
  end function ChartIDFace

  function ChartIDLength(ChartID)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    integer :: ChartIDLength
    ChartIDLength = BC_Get_Length(BC_Chart_Atlas(ChartID), BC_CHART_AtlasIndex(ChartID))
  end function ChartIDLength

  function ChartIDOffset(ChartID)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    integer :: ChartIDOffset
    ChartIDOffset = BC_Get_Offset(BC_Chart_Atlas(ChartID), BC_CHART_AtlasIndex(ChartID))
  end function ChartIDOffset

  function GetChartIDValues(ChartID) RESULT(Values)
    real(r8), dimension(:,:), POINTER :: Values
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    Values => BC_Get_Values(BC_Chart_Atlas(ChartID), BC_Chart_AtlasIndex(ChartID))
  end function GetChartIDValues

  function GetChartIDValueIndex(ChartID) RESULT(ValueIndex)
    integer,dimension(:), POINTER :: ValueIndex
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    ValueIndex => BC_Get_ValueIndex(BC_Chart_Atlas(ChartID), BC_Chart_AtlasIndex(ChartID))
  end function GetChartIDValueIndex

  function GetChartIDUseFunction(ChartID) RESULT(UseFunction)
    logical,dimension(:), POINTER :: UseFunction
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    UseFunction => BC_Get_UseFunction(BC_Chart_Atlas(ChartID), BC_Chart_AtlasIndex(ChartID))
  end function GetChartIDUseFunction

  function GetChartIDPositions(ChartID) RESULT(Positions)
    real(r8), dimension(:,:), POINTER :: Positions
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    Positions => BC_Get_Positions(BC_Chart_Atlas(ChartID), BC_Chart_AtlasIndex(ChartID))
  end function GetChartIDPositions

  subroutine SetChartID(ChartID, Atlas, AtlasIndex)
    ! Set the ChartID to point at the appropriate Atlas & AtlasIndex
    type(BC_Chart_ID), target, intent(OUT) :: ChartID
    type(BC_Atlas), TARGET, intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    ChartID%AtlasIndex = AtlasIndex
    ChartID%Atlas => Atlas
  end subroutine SetChartID

  subroutine SetChartChartID(Chart, ChartID)
    ! Set the chart according to ChartID information
    type(BC_Chart),    target, intent(OUT)    :: Chart
    type(BC_Chart_ID), target, intent(IN) :: ChartID

    ! Local variables
    type(BC_Atlas), POINTER :: Atlas
    real(r8), dimension(:,:), pointer :: Values
    integer :: AtlasIndex

    Atlas => BC_Chart_Atlas(ChartID)
    AtlasIndex = BC_CHART_AtlasIndex(ChartID)
    Values => BC_Get_Values(Atlas, AtlasIndex)

    Call BC_Chart_Set_Chart(Chart, LENGTH = BC_Get_Length(Atlas, AtlasIndex),       &
                                   DIMENSIONALITY = DIMENSIONALITY(Atlas), &
                                   DOF    = BC_Get_DOF(Atlas),             &
                                   CELL   = BC_Get_Cell(Atlas, AtlasIndex),       &
                                   FACE   = BC_Get_Face(Atlas, AtlasIndex),       &
                                   VALUES = Values, &
                                   VALUEINDEX = BC_Get_ValueIndex(Atlas, AtlasIndex), &
                                   USEFUNCTION = BC_Get_UseFunction(Atlas, AtlasIndex), &
                                   POSITIONS = BC_Get_Positions(Atlas, AtlasIndex))

  end subroutine SetChartChartID
    

  function SetNextChartID(ChartID, NextChartID)
    ! Get the chart ID following the input chart ID
    ! Returns the data into NextChart
    ! Returns .TRUE. if there is such a chart
    ! Returns .FALSE. if there is not a next chart
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    type(BC_Chart_ID), target, intent(OUT) :: NextChartID
    logical :: SetNextChartID

    ! Local variables
    integer :: CurrentIndex, NextIndex, NumberOfCharts
    type(BC_Atlas), POINTER :: Atlas
    
    ! Information about the input chart
    CurrentIndex = BC_Chart_AtlasIndex(ChartID)
    Atlas => BC_Chart_Atlas(ChartID)
    NumberOfCharts = BC_NumberOfCharts(Atlas)

    ! Get the next AtlasIndex
    NextIndex = CurrentIndex + 1
    if (NextIndex  <= NumberOfCharts) then
       SetNextChartID = .TRUE.
       Call BC_Chart_Set_Chart(NextChartID, Atlas, NextIndex)
    else
       SetNextChartID = .FALSE.
       NextChartID = BC_Invalid_Chart_ID()
    end if

  end function SetNextChartID

  function ChartIDMatches(ChartID, Cell, Face)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    integer, target, intent(IN) :: Cell
    integer, target, intent(IN), OPTIONAL :: Face
    logical :: ChartIDMatches 
    
    if (BC_Chart_ID_Is_Valid(ChartID)) then
       
       if (PRESENT(Face)) then
          ChartIDMatches = (Cell == BC_Chart_Cell(ChartID)) &
               .AND.       (Face == BC_Chart_Face(ChartID))
       else
          ChartIDMatches = (Cell == BC_Chart_Cell(ChartID))
       end if
    else
       ChartIDMatches = .FALSE.
    end if
  end function ChartIDMatches

  function ValidChart(ChartID)
    type(BC_Chart_ID), target, intent(IN) :: ChartID
    logical :: ValidChart

    ! To make this routine fast, we only check that the index is not 
    ! BC_INVALID_ATLASINDEX.  We assume any other index is valid.
    ! Obviously it is possible to violate this condition.  The rest 
    ! of the code is supposed to keep this condition true.
    

    ValidChart = ChartID%AtlasIndex /= BC_INVALID_ATLASINDEX
  end function ValidChart

  subroutine AppendChartAtlas(Atlas, Chart)
    ! Append the Chart data to the Atlas
    type(BC_Atlas), target, intent(INOUT) :: Atlas
    type(BC_Chart), target, intent(IN) :: Chart

    ! Local variables
    integer,dimension( : ),pointer :: ValueIndex
    logical,dimension( : ),pointer :: UseFunction
    real(r8), dimension(:,:), pointer :: Values
    real(r8), dimension(:,:), pointer :: Positions
    integer :: Cell, Face
    
    Cell = BC_Chart_Cell(Chart)
    Face = BC_Chart_Face(Chart)
    Values => BC_Chart_Values(Chart)
    ValueIndex => BC_Chart_ValueIndex(Chart)
    UseFunction=> BC_Chart_UseFunction(Chart)
    Positions => BC_Chart_Positions(Chart)
    call APPEND(Atlas, VALUES     = Values,     &
                       VALUEINDEX = ValueIndex, &
                       USEFUNCTION= UseFunction,&
                       POSITIONS  = Positions,  &
                       CELL       = Cell,       &
                       FACE       = Face)

  end subroutine AppendChartAtlas
  
END Module BC_CHARTS_Atlases

