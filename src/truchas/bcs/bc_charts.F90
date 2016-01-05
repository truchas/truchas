!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_CHARTS
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition chart
  !   operations in Telluride.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  use bc_enum_types
  Implicit None
  Private

  PUBLIC :: BC_Chart, BC_Chart_Data, BC_Chart_Spec
  PUBLIC :: BC_Chart_Length
  PUBLIC :: BC_Chart_Size
  PUBLIC :: BC_Chart_Empty_Chart
  PUBLIC :: BC_Chart_Matches
  PUBLIC :: BC_Chart_Set_Chart
  PUBLIC :: BC_CHART_VALUES
  PUBLIC :: BC_CHART_VALUEINDEX
  PUBLIC :: BC_CHART_USEFUNCTION
  PUBLIC :: BC_CHART_POSITIONS
  PUBLIC :: BC_Chart_Cell
  PUBLIC :: BC_Chart_Face
  PUBLIC :: INITIALIZE, ALLOC, FREE
  PUBLIC :: BC_Get_DOF
  PUBLIC :: BC_Set_DOF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Type to hold data for a chart.
  type BC_Chart_Data
     PRIVATE
     ! DataSize is the total size of the arrays.  NOT necessarily the length
     ! of the chart data stored in the array (DataSize may be > than valid length)
     integer :: DataSize
     ! Dimensionality is the size of the first dimension of Position
     integer :: Dimensionality
     integer :: DegreesOfFreedom
     integer,  pointer, dimension(:)   :: ValueIndex
     logical,  pointer, dimension(:)   :: UseFunction
     real(r8), pointer, dimension(:,:) :: Values
     real(r8), pointer, dimension(:,:) :: Position
  end type BC_Chart_Data

  ! Type to specify information about a chart.
  type BC_Chart_Spec
     PRIVATE
     integer :: Face
     integer :: Cell
     integer :: Offset
     ! Length is the number of data items in this chart.  
     integer :: Length
  end type BC_Chart_Spec

  ! A chart
  type BC_CHART
     PRIVATE
     type(BC_Chart_Spec) :: Spec
     type(BC_Chart_Data) :: Data
  end type BC_CHART
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE InitChartData
     MODULE PROCEDURE InitChartSpec
     MODULE PROCEDURE InitChart
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE AllocChartData
     MODULE PROCEDURE AllocChart
  END INTERFACE

  INTERFACE FREE
     MODULE PROCEDURE FreeChartData
     MODULE PROCEDURE FreeChart
  END INTERFACE

  INTERFACE BC_CHART_EMPTY_CHART
     MODULE PROCEDURE EmptyChart
  END INTERFACE
  
  INTERFACE BC_CHART_LENGTH
     MODULE PROCEDURE ChartSpecLength
     MODULE PROCEDURE ChartLength
  END INTERFACE

  INTERFACE BC_CHART_SIZE
     MODULE PROCEDURE ChartDataSize
     MODULE PROCEDURE ChartSize     
  END INTERFACE
  
  INTERFACE BC_CHART_CELL
     MODULE PROCEDURE ChartSpecCell
     MODULE PROCEDURE ChartCell
  END INTERFACE

  INTERFACE BC_CHART_FACE
     MODULE PROCEDURE ChartSpecFace
     MODULE PROCEDURE ChartFace
  END INTERFACE
  
  INTERFACE BC_CHART_VALUES
     MODULE PROCEDURE ChartDataValues
     MODULE PROCEDURE ChartValues
  END INTERFACE
  
  INTERFACE BC_CHART_VALUEINDEX
     MODULE PROCEDURE ChartDataValueIndex
     MODULE PROCEDURE ChartValueIndex
  END INTERFACE
  
  INTERFACE BC_CHART_USEFUNCTION
     MODULE PROCEDURE ChartDataUseFunction
     MODULE PROCEDURE ChartUseFunction
  END INTERFACE
  
  INTERFACE BC_CHART_POSITIONS
     MODULE PROCEDURE ChartDataPositions
     MODULE PROCEDURE ChartPositions
  END INTERFACE

  INTERFACE BC_Chart_Matches
     MODULE PROCEDURE ChartMatches
  END INTERFACE

  INTERFACE BC_Chart_Set_Chart
     MODULE PROCEDURE setChartData
     MODULE PROCEDURE setChartSpec
     MODULE PROCEDURE SetChart
  END INTERFACE

  INTERFACE BC_Get_DOF
     MODULE PROCEDURE ChartDataGetDOF
     MODULE PROCEDURE ChartGetDOF
  END INTERFACE
  INTERFACE BC_Set_DOF
     MODULE PROCEDURE ChartDataSetDOF
  END INTERFACE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(BC_Chart_Spec), PARAMETER :: Not_A_Chart = BC_Chart_Spec(-1,-1,-1,-1)
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate Chart_Data
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function ChartDataSize(ChartData)
    type(BC_Chart_Data), intent(IN) :: ChartData
    integer :: ChartDataSize
    ChartDataSize = ChartData%DataSize
  end function ChartDataSize

  subroutine InitChartData(ChartData)
    ! Initialize fields in ChartData
    type(BC_Chart_Data), intent(OUT) :: ChartData

    NULLIFY(ChartData%Values)
    NULLIFY(ChartData%ValueIndex)
    NULLIFY(ChartData%UseFunction)
    NULLIFY(ChartData%Position)
    ChartData%DataSize = BC_INVALID_SIZE
    ChartData%Dimensionality = BC_INVALID_DIMENSIONALITY
    ChartData%DegreesOfFreedom = BC_INVALID_SIZE

  end subroutine InitChartData

  subroutine AllocChartData(ChartData, SIZE, DIMENSIONALITY, DOF)
    ! Allocate fields in ChartData
    type(BC_Chart_Data), intent(OUT) :: ChartData
    integer, intent(IN) :: Size, Dimensionality, DOF

    ALLOCATE(ChartData%Values(DOF,Size))
    ALLOCATE(ChartData%ValueIndex(Size))
    ALLOCATE(ChartData%UseFunction(Size))
    ALLOCATE(ChartData%Position(Dimensionality,Size))
    ChartData%DataSize = Size
    ChartData%Dimensionality     = Dimensionality
    ChartData%DegreesOfFreedom = DOF
  end subroutine AllocChartData

  subroutine FreeChartData(ChartData)
    ! Free fields in ChartData
    type(BC_Chart_Data), intent(INOUT) :: ChartData

    if (ASSOCIATED(ChartData%Values)) then
       DEALLOCATE(ChartData%Values)
    end if
    if (ASSOCIATED(ChartData%ValueIndex)) then
       DEALLOCATE(ChartData%ValueIndex)
    end if
    if (ASSOCIATED(ChartData%UseFunction)) then
       DEALLOCATE(ChartData%UseFunction)
    end if
    if (ASSOCIATED(ChartData%Position)) then
       DEALLOCATE(ChartData%Position)
    end if
    call INITIALIZE(ChartData)

  end subroutine FreeChartData

  subroutine setChartData(ChartData, Values, ValueIndex, UseFunction, Positions)
    ! Set chart data, assumes space is already allocated
    type(BC_Chart_Data),  intent(INOUT) :: ChartData
    real(r8), dimension(:,:), intent(IN) :: Values
    integer,dimension(:), intent(IN) :: ValueIndex
    logical,dimension(:), intent(IN) :: UseFunction
    real(r8), dimension(:,:), intent(IN) :: Positions

#ifdef ADD_DEBUG_CODE_HERE
    ! Local variables
    integer :: DataSize

    ! Check that the ChartData Size is large enough
    DataSize = BC_Chart_Data_Size(ChartData)
    if (DataSize < SIZE(Values,2)) then
       print *, 'FATAL: ChartData not large enough in setChartData'
    end if
    ! Check that the ChartData has the right DOFs
    if (BC_Get_DOF(ChartData) /= SIZE(Values,1)) then
       print *, 'FATAL: ChartData does not have the right number of Degrees of Freedom in setChartData'
    end if
#endif

    ChartData%Values     = Values
    ChartData%ValueIndex   = ValueIndex
    ChartData%UseFunction  = UseFunction
    ChartData%Position     = Positions
  end subroutine setChartData
  
  function ChartDataValues(ChartData) RESULT(Values)
    ! Return a pointer to the Values field of chart data
    type(BC_Chart_Data), intent(IN), TARGET :: ChartData
    real(r8), dimension(:,:), pointer :: VALUES
    Values => ChartData%Values
  end function ChartDataValues
  
  function ChartDataValueIndex(ChartData) RESULT(ValueIndex)
    ! Return a pointer to the ValueIndex field of chart data
    type(BC_Chart_Data), intent(IN), TARGET :: ChartData
    integer, dimension(:), pointer :: VALUEINDEX
    ValueIndex => ChartData%ValueIndex
  end function ChartDataValueIndex
  
  function ChartDataUseFunction(ChartData) RESULT(UseFunction)
    ! Return a pointer to the UseFunction field of chart data
    type(BC_Chart_Data), intent(IN), TARGET :: ChartData
    logical, dimension(:), pointer :: UseFunction
    UseFunction => ChartData%UseFunction
  end function ChartDataUseFunction
  
  function ChartDataPositions(ChartData) RESULT(Positions)
    ! Return a pointer to the Positions field of chart data
    type(BC_Chart_Data), intent(IN), TARGET :: ChartData
    real(r8), dimension(:,:), pointer :: Positions
    Positions => ChartData%Position
  end function ChartDataPositions
  

  subroutine ChartDataSetDOF(ChartData, DOF)
    type(BC_Chart_Data), intent(INOUT) :: ChartData
    integer, intent(IN)    :: DOF
    ChartData%DegreesOfFreedom = DOF
  end subroutine ChartDataSetDOF

  function ChartDataGetDOF(ChartData) RESULT(DOF)
    type(BC_Chart_Data), intent(IN) :: ChartData
    integer :: DOF
    DOF = ChartData%DegreesOfFreedom
  end function ChartDataGetDOF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate Chart_Spec
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine InitChartSpec(ChartSpec)
    ! Initialize a chart spec to be nothing
    type(BC_Chart_Spec), intent(OUT) :: ChartSpec
    ChartSpec = Not_A_Chart
  end subroutine InitChartSpec

  function ChartSpecLength(ChartSpec)
    type(BC_Chart_Spec), intent(IN) :: ChartSpec
    integer :: ChartSpecLength
    ChartSpecLength = ChartSpec%Length
  end function ChartSpecLength
  
  function ChartSpecFace(ChartSpec)
    type(BC_Chart_Spec), intent(IN) :: ChartSpec
    integer :: ChartSpecFace
    ChartSpecFace = ChartSpec%Face
  end function ChartSpecFace

  function ChartSpecCell(ChartSpec)
    type(BC_Chart_Spec), intent(IN) :: ChartSpec
    integer :: ChartSpecCell
    ChartSpecCell = ChartSpec%Cell
  end function ChartSpecCell

  subroutine setChartSpec(ChartSpec, Length, Cell, Face)
    ! Set the fields in a Chart Spec
    ! Offset is not set, since only relevant with an Atlas
    type(BC_Chart_Spec), intent(OUT) :: ChartSpec
    integer,             intent(IN) :: Length, Cell
    integer, OPTIONAL,   intent(IN) :: Face

    ChartSpec%Length = Length
    ChartSpec%Offset = BC_INVALID_OFFSET
    ChartSpec%Cell   = Cell
    if (PRESENT(Face)) then
       ChartSpec%Face = Face
    else
       ChartSpec%Face = BC_INVALID_FACE
    end if

  end subroutine setChartSpec
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate a Chart
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function EmptyChart()
    ! Return a pointer to an empty chart
    type(BC_Chart), TARGET, SAVE :: Empty_Chart
    logical, SAVE :: Initialized = .FALSE.
    type(BC_Chart), POINTER      :: EmptyChart
    if (.NOT. Initialized) then
       Call INITIALIZE(Empty_Chart)
       Initialized = .TRUE.
    end if
    EmptyChart => Empty_Chart
  end function EmptyChart

  subroutine InitChart(Chart)
    ! Initialize a BC_Chart
    type(BC_Chart), intent(OUT) :: Chart
    Call INITIALIZE(Chart%Data)
    Call INITIALIZE(Chart%Spec)
  end subroutine InitChart
  
  subroutine AllocChart(Chart, SIZE, Dimensionality, DOF)
    ! Allocate the chart data fields of a chart
    type(BC_Chart), intent(INOUT) :: Chart
    integer, intent(IN) :: SIZE, Dimensionality, DOF
    Call ALLOC(Chart%Data, SIZE, DIMENSIONALITY, DOF)
  end subroutine AllocChart

  subroutine FreeChart(Chart)
    ! Allocate the chart data fields of a chart
    type(BC_Chart), intent(INOUT) :: Chart
    Call FREE(Chart%Data)
  end subroutine FreeChart

  function ChartCell(Chart)
    type(BC_Chart), intent(IN) :: Chart
    integer :: ChartCell
    ChartCell = BC_Chart_Cell(Chart%Spec)
  end function ChartCell

  function ChartFace(Chart)
    type(BC_Chart), intent(IN) :: Chart
    integer :: ChartFace
    ChartFace = BC_Chart_Face(Chart%Spec)
  end function ChartFace

  function ChartLength(Chart)
    ! Returns the length of the chart = number of valid data points
    type(BC_Chart), intent(IN) :: Chart
    integer :: ChartLength
    ChartLength = BC_Chart_Length(Chart%Spec)
  end function ChartLength

  function ChartSize(Chart)
    ! Returns the size of the chart data arrays, which may be larger
    ! than ChartLength.
    type(BC_Chart), intent(IN) :: Chart
    integer :: ChartSize
    ChartSize = BC_Chart_Size(Chart%Data)
  end function ChartSize

  function ChartGetDOF(Chart) RESULT(DOF)
    ! Returns the size of the chart data arrays, which may be larger
    ! than ChartLength.
    type(BC_Chart), intent(IN) :: Chart
    integer :: DOF
    DOF = BC_Get_DOF(Chart%Data)
  end function ChartGetDOF

  function ChartValues(Chart) RESULT(Values)
    ! Return a pointer to the VALUES field of a chart
    type(BC_Chart), intent(IN), TARGET :: Chart
    real(r8), dimension(:,:), pointer :: VALUES
    VALUES => BC_CHART_VALUES(Chart%Data)
  end function ChartValues
  
  function ChartValueIndex(Chart) RESULT(ValueIndex)
    ! Return a pointer to the VALUEINDEX field of a chart
    type(BC_Chart), intent(IN), TARGET :: Chart
    integer, dimension(:), pointer :: VALUEINDEX
    VALUEINDEX => BC_CHART_VALUEINDEX(Chart%Data)
  end function ChartValueIndex
  
  function ChartUseFunction(Chart) RESULT(UseFunction)
    ! Return a pointer to the USEFUNCTION field of a chart
    type(BC_Chart), intent(IN), TARGET :: Chart
    logical, dimension(:), pointer :: USEFUNCTION
    USEFUNCTION => BC_CHART_USEFUNCTION(Chart%Data)
  end function ChartUseFunction
  
  function ChartPositions(Chart) RESULT(Positions)
    ! Return a pointer to the POSITIONS field of a chart
    type(BC_Chart), intent(IN), TARGET :: Chart
    real(r8), dimension(:,:), pointer :: POSITIONS
    POSITIONS => BC_CHART_POSITIONS(Chart%Data)
  end function ChartPositions
  
  subroutine SetChart(Chart, Length, Dimensionality, DOF, Cell, Face, Values, ValueIndex, UseFunction, Positions)
    ! Set the chart according to given data.
    type(BC_Chart), intent(OUT)    :: Chart
    integer, intent(IN) :: Length, Dimensionality, DOF, Cell
    integer, OPTIONAL, intent(IN) :: Face
    real(r8), dimension(:,:),intent(IN) :: Values
    integer,dimension(:),intent(IN) :: ValueIndex
    logical,dimension(:),intent(IN) :: UseFunction
    real(r8), dimension(:,:),intent(IN) :: Positions

    Call BC_Chart_Set_Chart(Chart%Spec, LENGTH = Length, &
                                        CELL   = Cell,   &
                                        FACE   = Face)

    ! Allocate a chart with exactly the amount of space for this chart
    Call ALLOC(Chart%Data, SIZE=Length, &
                           DIMENSIONALITY = Dimensionality, &
                           DOF = DOF)
    Call BC_Chart_Set_Chart(Chart%Data, VALUES     = Values,     &
                                        VALUEINDEX = ValueIndex, &
                                        USEFUNCTION= UseFunction,&
                                        POSITIONS  = Positions)

  end subroutine SetChart

  function ChartMatches(Chart, Cell, Face)
    type(BC_Chart), intent(IN) :: Chart
    integer, intent(IN) :: Cell
    integer, intent(IN), OPTIONAL :: Face
    logical :: ChartMatches 
    
    if (PRESENT(Face)) then
       ChartMatches = (Cell == BC_Chart_Cell(Chart)) &
            .AND.     (Face == BC_Chart_Face(Chart))
    else
       ChartMatches = (Cell == BC_Chart_Cell(Chart))
    end if
  end function ChartMatches

END Module BC_CHARTS

