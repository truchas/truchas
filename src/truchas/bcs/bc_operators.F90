!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_OPERATORS
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition operators
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use bc_enum_types
  use bc_atlases,         only: DIMENSIONALITY, BC_NumberOfCharts, COLLATE, BC_Get_DOF
  use bc_charts_atlases,  only: BC_Chart_ID, BC_Chart_Set_Chart, &
                                BC_Chart_Next_Chart, BC_Atlas,   &
                                BC_Invalid_Chart_ID, INITIALIZE, FREE
  use bc_regions, only: BC_Region, INITIALIZE, FREE, DIMENSIONALITY, COLLATE
  Implicit None
  Private

  PUBLIC :: BC_Operator
  PUBLIC :: INITIALIZE, FREE
  PUBLIC :: BC_OP_Get_Chart
  PUBLIC :: BC_OP_Get_Region
  PUBLIC :: BC_OP_Get_Atlas
  PUBLIC :: BC_OP_Get_ID
  PUBLIC :: BC_Op_Start_Search
  PUBLIC :: BC_Get_Chart
  PUBLIC :: BC_Get_Region
  PUBLIC :: BC_Get_Atlas
  PUBLIC :: DIMENSIONALITY
  PUBLIC :: COLLATE
  PUBLIC :: BC_Op_Set_State
  PUBLIC :: BC_Op_Is_Active
  PUBLIC :: BC_Get_Name
  PUBLIC :: BC_Set_Name
  PUBLIC :: BC_Get_DOF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Type to store the data for the BC operator
  ! Type to identify an operator
  type BC_Operator
     PRIVATE
     integer            :: OP_ID
     type (BC_Region)   :: Region
     type (BC_Atlas)    :: Atlas
     type (BC_Chart_ID), POINTER :: CurrentChart_ID
     type (BC_Chart_ID), POINTER :: NextChart_ID
     ! This is true if this operator is active.  
     integer            :: State
     character(LEN=BC_STRING_LEN) :: Name = BC_NO_NAME
  end type BC_Operator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE InitOperator
  END INTERFACE

  INTERFACE FREE
     MODULE PROCEDURE FreeOperator
  END INTERFACE

  INTERFACE BC_OP_GET_CHART
     MODULE PROCEDURE GetBCChartID
  END INTERFACE

  INTERFACE BC_GET_CHART
     MODULE PROCEDURE GetBCChartID
  END INTERFACE
  
  INTERFACE BC_OP_GET_REGION
     MODULE PROCEDURE GetRegion
  END INTERFACE

  INTERFACE BC_GET_REGION
     MODULE PROCEDURE GetRegion
  END INTERFACE

  INTERFACE BC_OP_GET_ATLAS
     MODULE PROCEDURE GetAtlas
  END INTERFACE

  INTERFACE BC_GET_ATLAS
     MODULE PROCEDURE GetAtlas
  END INTERFACE

  INTERFACE BC_OP_Get_ID
     MODULE PROCEDURE GetOpID
  END INTERFACE

  INTERFACE BC_OP_GET_CURRENT_CHART
     MODULE PROCEDURE GetCurrentChart_ID
  END INTERFACE
  
  INTERFACE BC_OP_SET_CURRENT_CHART
     MODULE PROCEDURE SetCurrentChart_ID
  END INTERFACE

  INTERFACE BC_OP_GET_NEXT_CHART
     MODULE PROCEDURE GetNextChart_ID
  END INTERFACE
  
  INTERFACE BC_OP_SET_NEXT_CHART
     MODULE PROCEDURE SetNextChart_ID
  END INTERFACE

  INTERFACE BC_Op_Start_Search
     MODULE PROCEDURE StartChartSearch
  END INTERFACE
       
  INTERFACE BC_Op_Step_Chart
     MODULE PROCEDURE StepChartIDs
  END INTERFACE

  INTERFACE BC_Op_Set_State
     MODULE PROCEDURE SetState
  END INTERFACE

  INTERFACE BC_Op_Is_Active
     MODULE PROCEDURE ActiveP
  END INTERFACE

  INTERFACE DIMENSIONALITY
     MODULE PROCEDURE GetOpDimensionality
  END INTERFACE

  INTERFACE COLLATE
     MODULE PROCEDURE OpCollate
  END INTERFACE
  
  INTERFACE BC_Get_Name
     MODULE PROCEDURE OpGetName
  END INTERFACE
  INTERFACE BC_Set_Name
     MODULE PROCEDURE OpSetName
  END INTERFACE

  INTERFACE BC_Get_DOF
     MODULE PROCEDURE OpGetDOF
  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  subroutine InitOperator(BC_Op, OP_ID) 
    type(BC_Operator), target, intent(INOUT) :: BC_Op
    integer,           intent(IN   ) :: OP_ID

    BC_Op%OP_ID = OP_ID
    call INITIALIZE(BC_Op%Region)
    call INITIALIZE(BC_Op%Atlas)
    ALLOCATE(BC_Op%CurrentChart_ID)
    call INITIALIZE(BC_Op%CurrentChart_ID)
    ALLOCATE(BC_Op%NextChart_ID)
    call INITIALIZE(BC_Op%NextChart_ID)
    call BC_Op_Set_State(BC_Op, BC_OP_INVALID_STATE)
    call BC_Set_Name(BC_Op, BC_NO_NAME)

  end subroutine InitOperator

  subroutine FreeOperator(BC_OP)
    ! Free all memory used by the operator
    type(BC_Operator), target, intent(INOUT) :: BC_Op

    BC_Op%OP_ID = BC_INVALID_ID
    call BC_Op_Set_State(BC_Op, BC_OP_INVALID_STATE)
    call BC_Set_Name(BC_Op, BC_NO_NAME)
    call FREE(BC_Op%Region)
    call FREE(BC_Op%Atlas)
    if (ASSOCIATED(BC_Op%CurrentChart_ID)) then
       DEALLOCATE(BC_Op%CurrentChart_ID)
    end if
    if (ASSOCIATED(BC_Op%NextChart_ID)) then
       DEALLOCATE(BC_Op%NextChart_ID)
    end if
    

  end subroutine FreeOperator

  function GetRegion(BC_Op) RESULT(Region)
    type (BC_Operator), intent(IN), TARGET :: BC_Op
    type (BC_Region),   POINTER            :: Region

    Region => BC_Op%Region

  end function GetRegion
      
  function GetAtlas(BC_Op) RESULT(Atlas)
    type (BC_Operator), intent(IN), TARGET :: BC_Op
    type (BC_Atlas),   POINTER            :: Atlas

    Atlas => BC_Op%Atlas

  end function GetAtlas
      
  function GetOpDimensionality(BC_Op) RESULT(Dims)
    type (BC_Operator), target, intent(IN) :: BC_Op
    integer                        :: Dims

    type (BC_Atlas),   POINTER            :: Atlas
    
    Atlas => BC_OP_Get_Atlas(BC_Op)
    Dims = DIMENSIONALITY(Atlas)

  end function GetOpDimensionality

  function GetOpID(BC_Op) RESULT(ID)
    type (BC_Operator), target, intent(IN) :: BC_Op
    integer                        :: ID

    ID = BC_Op%OP_ID

  end function GetOpID

  subroutine SetState(BC_Op, State)
    type (BC_Operator), target, intent(INOUT) :: BC_Op
    integer                           :: State

    BC_Op%State = State
  end subroutine SetState

  function ActiveP(BC_Op)
    type (BC_Operator), target, intent(IN) :: BC_Op
    logical              :: ActiveP

    ActiveP = (BC_Op%State == BC_OP_ACTIVE)
  end function ActiveP

  function GetCurrentChart_ID(BC_Op) RESULT(Chart_ID)
    type (BC_Operator), intent(INOUT), TARGET :: BC_Op
    type (BC_Chart_ID),   POINTER            :: Chart_ID
    
    Chart_ID => BC_Op%CurrentChart_ID

  end function GetCurrentChart_ID

  subroutine SetCurrentChart_ID(BC_Op, Chart_ID)
    type (BC_Operator), target, intent(INOUT) :: BC_Op
    type (BC_Chart_ID),    TARGET,        &
                        intent(IN   ) :: Chart_ID
    
    BC_Op%CurrentChart_ID => Chart_ID
  end subroutine SetCurrentChart_ID

  function GetNextChart_ID(BC_Op) RESULT(Chart_ID)
    type (BC_Operator), intent(IN), TARGET :: BC_Op
    type (BC_Chart_ID),   POINTER            :: Chart_ID
    
    Chart_ID => BC_Op%NextChart_ID
  end function GetNextChart_ID

  subroutine SetNextChart_ID(BC_Op, Chart_ID)
    type (BC_Operator), target, intent(INOUT) :: BC_Op
    type (BC_Chart_ID),    TARGET,        &
                        intent(IN   ) :: Chart_ID
    
    BC_Op%NextChart_ID => Chart_ID
    
  end subroutine SetNextChart_ID

  subroutine StartChartSearch(BC_Op)
    ! Initialize Chart ID's inside BC_Op to prepare to step through all charts
    type (BC_Operator), target, intent(INOUT) :: BC_Op

    ! Local variables
    type (BC_Chart_ID), POINTER :: CurrentID, NextID
    logical :: foundNext

    ! Get the current chart
    CurrentID => BC_Op_Get_Current_Chart(BC_Op)
    ! Get the next chart
    NextID => BC_Op_Get_Next_Chart(BC_Op)

    ! If there are no charts, then can never search, so set both current and next ID's to invalide

    if (BC_NumberOfCharts(BC_OP_Get_Atlas(BC_Op)) <= 0) then
       CurrentID   = BC_Invalid_Chart_ID()
       NextID      = BC_Invalid_Chart_ID()
    else
       ! Set the currentID to point to the first entry, which has index = 1
       call BC_Chart_Set_Chart(CurrentID, ATLAS=BC_OP_Get_Atlas(BC_Op), AtlasIndex = 1)
       ! Set the next chart ID
       FoundNext = BC_Chart_Next_Chart(CurrentID, NextID)
    end if

  end subroutine StartChartSearch

  subroutine StepChartIDs(BC_Op)
    ! Advance the chartIDs.  The data in the "NextChartID" becomes the "CurrentChartID",
    ! and we get new data for the NextChartID
    use bc_charts_atlases,  only: BC_Chart_ID, BC_Chart_Next_Chart
    type (BC_OPERATOR), target, intent(INOUT) :: BC_Op
    
    ! Local variables
    type (BC_Chart_ID), POINTER :: CurrentChartID, NextChartID, TempChartID
    logical :: FoundNextChartID

    ! We want to use the existing memory, which we assume is already allocated.
    ! so save a pointer to the current chart id memory, which we will reuse
    TempChartID => BC_OP_Get_Current_Chart(BC_Op)
    
    ! Point the CurrentChartID field of BC_Op to the data for NextChartID.
    CurrentChartID => BC_Op_Get_Next_Chart(BC_Op)
    call BC_OP_Set_Current_Chart(BC_Op, CurrentChartID)

    ! Now we want to use the memory space that was being used for the CurrentChartID,
    ! but isn't since we reset that field.  We need the data in the CurrentChartID
    ! so we can advance from there, but we already have that.  We need the space
    ! that was used by the current chart, before we re-pointered.  That space
    ! is pointed at by TempChartID.
    NextChartID => TempChartID

    ! Now we need to get new data into there.
    FoundNextChartID = BC_Chart_Next_Chart(CurrentChartID, NextChartID)
    
    ! Finally, poke this into the chartID
    call BC_Op_Set_Next_Chart(BC_Op, NextChartID)
  end subroutine StepChartIDs

  function GetBCChartID(ChartID, Operator, Cell, Face) RESULT(FoundChartID)
    use bc_charts_atlases, only: BC_Chart_ID, BC_Chart_Matches, &
                                 BC_Invalid_Chart_ID
    type (BC_Chart_ID),    POINTER    :: ChartID
    type (BC_Operator), TARGET,     &
                        intent(INOUT) :: Operator
    integer, intent(IN)           :: Cell
    integer, intent(IN), &
                          OPTIONAL :: Face
    logical :: FoundChartID

    ! Local variables
    type (BC_Chart_ID), POINTER :: NextChartID => NULL()
    type (BC_Chart_ID), POINTER :: CurrentChartID => NULL()

    ! Look up chart ID for the Cell or (Face, Cell).
    ! If a chart ID is found FoundChart is TRUE and the
    ! Chart ID argument points to the found chart.
    ! IF a chart matching Cell or (Face,Cell) is not
    ! found, result of the function is FALSE
    ! and the chart ID points to an invalid chart.

    ! Get the current chart from the operator
    CurrentChartID => BC_Op_Get_Current_Chart(Operator)
    
    ! Look for a match
    FoundChartID = .FALSE.
    ! First see if the Cell or (Cell,Face) is the current chart.
    if (BC_Chart_Matches(CurrentChartID, Cell, Face)) then
       FoundChartID = .TRUE.
       ChartID        => CurrentChartID
    else
       ! If that didn't match, look to see if the next chart matches
       NextChartID => BC_Op_Get_Next_Chart(Operator)
       if (BC_Chart_Matches(NextChartID, Cell, Face)) then
          FoundChartID = .TRUE.
          ChartID        => NextChartID 
          ! If found a match, step to the next chart
          call BC_Op_Step_Chart(Operator)
       end if
    end if
    
    ! If match was not found, then return an empty chart
    if (.NOT. FoundChartID) then
       ChartID => BC_Invalid_Chart_ID()
    end if

  end function GetBCChartID

  subroutine OpCollate(Collated_Op, Local_Op)
    ! Collate all the components of Local_Op into Collated_Op.
    ! This routine assumes that Collated_Op has been appropriately allocated,
    ! so Collated_Op is an INOUT argument.
    type (BC_Operator), target, intent(INOUT) :: Collated_Op
    type (BC_Operator), target, intent(IN   ) :: Local_Op

    ! Local variables
    Collated_Op%OP_ID = Local_Op%Op_ID

    call COLLATE(Collated_Op%Region, Local_Op%Region)

    call COLLATE(Collated_Op%Atlas,  Local_Op%Atlas )

    call BC_Set_Name(Collated_Op, BC_Get_NAME(Local_Op))
  end subroutine OpCollate
    
  subroutine OpSetName(BC_Op, NAME)
    type (BC_Operator), target, intent(INOUT) :: BC_Op
    character(LEN=*), intent(IN) :: NAME
    BC_Op%Name = TRIM(NAME)
  END subroutine OpSetName

  function OpGetName(BC_Op) RESULT(NAME)
    type (BC_Operator), target, intent(IN   ) :: BC_Op
    character(LEN=BC_STRING_LEN) ::NAME
    NAME = TRIM(BC_Op%Name)
  END function OpGetName

  function OpGetDOF(BC_Op) RESULT(DOF)
    type (BC_Operator), target, intent(IN   ) :: BC_Op
    integer :: DOF
    DOF = BC_Get_DOF(BC_Op%Atlas)
  END function OpGetDOF

END Module BC_OPERATORS
