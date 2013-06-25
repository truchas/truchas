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
  use kind_module, only: int_kind, real_kind, log_kind
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
     integer( int_kind)                          :: DataSize
     ! Dimensionality is the size of the first dimension of Position
     integer( int_kind)                          :: Dimensionality
     integer( int_kind)                          :: DegreesOfFreedom
     integer( int_kind), pointer, dimension(:)   :: ValueIndex
     logical( log_kind), pointer, dimension(:)   :: UseFunction
     real(real_kind),    pointer, dimension(:,:) :: Values
     real(real_kind),    pointer, dimension(:,:) :: Position
  end type BC_Chart_Data

  ! Type to specify information about a chart.
  type BC_Chart_Spec
     PRIVATE
     integer( int_kind) :: Face
     integer( int_kind) :: Cell
     integer( int_kind) :: Offset
     ! Length is the number of data items in this chart.  
     integer( int_kind) :: Length
  end type BC_Chart_Spec

  ! A chart
  type BC_CHART
     PRIVATE
     type (BC_Chart_Spec) :: Spec
     type (BC_Chart_Data) :: Data
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

  type (BC_Chart_Spec), PARAMETER :: Not_A_Chart = BC_Chart_Spec(-1,-1,-1,-1)
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate Chart_Data
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function ChartDataSize(ChartData)
    implicit none
    type(BC_Chart_Data), intent(IN) :: ChartData
    integer( int_kind) :: ChartDataSize

    ChartDataSize = ChartData%DataSize
    return
  end function ChartDataSize

  subroutine InitChartData(ChartData)
    ! Initialize fields in ChartData
    type(BC_Chart_Data), intent(  OUT) :: ChartData

    NULLIFY(ChartData%Values)
    NULLIFY(ChartData%ValueIndex)
    NULLIFY(ChartData%UseFunction)
    NULLIFY(ChartData%Position)
    ChartData%DataSize = BC_INVALID_SIZE
    ChartData%Dimensionality     = BC_INVALID_DIMENSIONALITY
    ChartData%DegreesOfFreedom = BC_INVALID_SIZE

    return
  end subroutine InitChartData

  subroutine AllocChartData(ChartData, SIZE, DIMENSIONALITY, DOF)
    ! Allocate fields in ChartData
    type(BC_Chart_Data), intent(  OUT) :: ChartData
    integer( int_kind),             intent(IN   ) :: Size, Dimensionality, DOF

    ALLOCATE(ChartData%Values(DOF,Size))
    ALLOCATE(ChartData%ValueIndex(Size))
    ALLOCATE(ChartData%UseFunction(Size))
    ALLOCATE(ChartData%Position(Dimensionality,Size))
    ChartData%DataSize = Size
    ChartData%Dimensionality     = Dimensionality
    ChartData%DegreesOfFreedom = DOF
    return
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

    return
  end subroutine FreeChartData

  subroutine setChartData(ChartData, Values, ValueIndex, UseFunction, Positions)
    ! Set chart data, assumes space is already allocated
    implicit none
    type(BC_Chart_Data),  intent(INOUT) :: ChartData
    real(real_kind), dimension(:,:), intent(IN   ) :: Values
    integer( int_kind),dimension(:), intent(IN   ) :: ValueIndex
    logical( log_kind),dimension(:), intent(IN   ) :: UseFunction
    real(real_kind), dimension(:,:), intent(IN   ) :: Positions

#ifdef ADD_DEBUG_CODE_HERE
    ! Local variables
    integer( int_kind) :: DataSize

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
    return
  end subroutine setChartData
  
  function ChartDataValues(ChartData) RESULT(Values)
    ! Return a pointer to the Values field of chart data
    type (BC_Chart_Data), intent(IN), TARGET :: ChartData
    real(real_kind), dimension(:,:), pointer :: VALUES
    Values => ChartData%Values
    return
  end function ChartDataValues
  
  function ChartDataValueIndex(ChartData) RESULT(ValueIndex)
    ! Return a pointer to the ValueIndex field of chart data
    type (BC_Chart_Data), intent(IN), TARGET :: ChartData
    integer( int_kind),           dimension(:), pointer :: VALUEINDEX
    ValueIndex => ChartData%ValueIndex
    return
  end function ChartDataValueIndex
  
  function ChartDataUseFunction(ChartData) RESULT(UseFunction)
    ! Return a pointer to the UseFunction field of chart data
    type (BC_Chart_Data), intent(IN), TARGET :: ChartData
    logical( log_kind),           dimension(:), pointer :: UseFunction
    UseFunction => ChartData%UseFunction
    return
  end function ChartDataUseFunction
  
  function ChartDataPositions(ChartData) RESULT(Positions)
    ! Return a pointer to the Positions field of chart data
    type (BC_Chart_Data), intent(IN), TARGET :: ChartData
    real(real_kind), dimension(:,:), pointer :: Positions
    Positions => ChartData%Position
    return
  end function ChartDataPositions
  

  subroutine ChartDataSetDOF(ChartData, DOF)
    implicit none
    type (BC_Chart_Data), intent(INOUT) :: ChartData
    integer(int_kind),    intent(IN)    :: DOF
    ChartData%DegreesOfFreedom = DOF
    return
  end subroutine ChartDataSetDOF

  function ChartDataGetDOF(ChartData) RESULT(DOF)
    implicit none
    type (BC_Chart_Data), intent(IN) :: ChartData
    integer(int_kind)                :: DOF
    DOF = ChartData%DegreesOfFreedom
    return
  end function ChartDataGetDOF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate Chart_Spec
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine InitChartSpec(ChartSpec)
    implicit none
    ! Initialize a chart spec to be nothing
    type(BC_Chart_Spec), intent(  OUT) :: ChartSpec
    ChartSpec = Not_A_Chart
  end subroutine InitChartSpec

  function ChartSpecLength(ChartSpec)
    implicit none
    type(BC_Chart_Spec), intent(IN) :: ChartSpec
    integer( int_kind)                         :: ChartSpecLength
    ChartSpecLength = ChartSpec%Length
    return
  end function ChartSpecLength
  
  function ChartSpecFace(ChartSpec)
    implicit none
    type(BC_Chart_Spec), intent(IN) :: ChartSpec
    integer( int_kind)                         :: ChartSpecFace
    ChartSpecFace = ChartSpec%Face
    return
  end function ChartSpecFace

  function ChartSpecCell(ChartSpec)
    implicit none
    type(BC_Chart_Spec), intent(IN) :: ChartSpec
    integer( int_kind)                         :: ChartSpecCell
    ChartSpecCell = ChartSpec%Cell
    return
  end function ChartSpecCell

  subroutine setChartSpec(ChartSpec, Length, Cell, Face)
    ! Set the fields in a Chart Spec
    ! Offset is not set, since only relevant with an Atlas
    implicit none
    type(BC_Chart_Spec), intent(  OUT) :: ChartSpec
    integer( int_kind),             intent(IN   ) :: Length, Cell
    integer( int_kind), OPTIONAL,   intent(IN   ) :: Face

    ChartSpec%Length = Length
    ChartSpec%Offset = BC_INVALID_OFFSET
    ChartSpec%Cell   = Cell
    if (PRESENT(Face)) then
       ChartSpec%Face = Face
    else
       ChartSpec%Face = BC_INVALID_FACE
    end if

    return
  end subroutine setChartSpec
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Routines to access and manipulate a Chart
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function EmptyChart()
    ! Return a pointer to an empty chart
    implicit none
    type (BC_Chart), TARGET, SAVE :: Empty_Chart
    logical,                 SAVE :: Initialized = .FALSE.
    type (BC_Chart), POINTER      :: EmptyChart
    if (.NOT. Initialized) then
       Call INITIALIZE(Empty_Chart)
       Initialized = .TRUE.
    end if
    EmptyChart => Empty_Chart
    return
  end function EmptyChart

  subroutine InitChart(Chart)
    ! Initialize a BC_Chart
    implicit none
    type (BC_Chart), intent(OUT) :: Chart
    Call INITIALIZE(Chart%Data)
    Call INITIALIZE(Chart%Spec)
  end subroutine InitChart
  
  subroutine AllocChart(Chart, SIZE, Dimensionality, DOF)
    ! Allocate the chart data fields of a chart
    type(BC_Chart), intent(INOUT) :: Chart
    integer( int_kind),        intent(IN   ) :: SIZE, Dimensionality, DOF
    Call ALLOC(Chart%Data, SIZE, DIMENSIONALITY, DOF)
    return
  end subroutine AllocChart

  subroutine FreeChart(Chart)
    ! Allocate the chart data fields of a chart
    implicit none
    type(BC_Chart), intent(INOUT) :: Chart
    Call FREE(Chart%Data)
    return
  end subroutine FreeChart

  function ChartCell(Chart)
    implicit none
    type (BC_Chart), intent(IN)   :: Chart
    integer( int_kind)                       :: ChartCell
    ChartCell = BC_Chart_Cell(Chart%Spec)
    return
  end function ChartCell

  function ChartFace(Chart)
    implicit none
    type (BC_Chart), intent(IN)   :: Chart
    integer( int_kind)                       :: ChartFace
    ChartFace = BC_Chart_Face(Chart%Spec)
    return
  end function ChartFace

  function ChartLength(Chart)
    ! Returns the length of the chart = number of valid data points
    implicit none
    type (BC_Chart), intent(IN)   :: Chart
    integer( int_kind)                       :: ChartLength
    ChartLength = BC_Chart_Length(Chart%Spec)
    return
  end function ChartLength

  function ChartSize(Chart)
    ! Returns the size of the chart data arrays, which may be larger
    ! than ChartLength.
    implicit none
    type (BC_Chart), intent(IN)   :: Chart
    integer( int_kind)                       :: ChartSize
    ChartSize = BC_Chart_Size(Chart%Data)
    return
  end function ChartSize

  function ChartGetDOF(Chart) RESULT(DOF)
    ! Returns the size of the chart data arrays, which may be larger
    ! than ChartLength.
    implicit none
    type (BC_Chart), intent(IN)   :: Chart
    integer( int_kind)            :: DOF
    DOF = BC_Get_DOF(Chart%Data)
    return
  end function ChartGetDOF

  function ChartValues(Chart) RESULT(Values)
    ! Return a pointer to the VALUES field of a chart
    implicit none
    type(BC_Chart), intent(IN), TARGET :: Chart
    real(real_kind), dimension(:,:), pointer :: VALUES
    VALUES => BC_CHART_VALUES(Chart%Data)
    RETURN
  end function ChartValues
  
  function ChartValueIndex(Chart) RESULT(ValueIndex)
    ! Return a pointer to the VALUEINDEX field of a chart
    implicit none
    type(BC_Chart), intent(IN), TARGET :: Chart
    integer( int_kind),     dimension(:), pointer :: VALUEINDEX
    VALUEINDEX => BC_CHART_VALUEINDEX(Chart%Data)
    RETURN
  end function ChartValueIndex
  
  function ChartUseFunction(Chart) RESULT(UseFunction)
    ! Return a pointer to the USEFUNCTION field of a chart
    implicit none
    type(BC_Chart), intent(IN), TARGET :: Chart
    logical( log_kind),     dimension(:), pointer :: USEFUNCTION
    USEFUNCTION => BC_CHART_USEFUNCTION(Chart%Data)
    RETURN
  end function ChartUseFunction
  
  function ChartPositions(Chart) RESULT(Positions)
    ! Return a pointer to the POSITIONS field of a chart
    implicit none
    type(BC_Chart), intent(IN), TARGET :: Chart
    real(real_kind), dimension(:,:), pointer :: POSITIONS
    POSITIONS => BC_CHART_POSITIONS(Chart%Data)
    RETURN
  end function ChartPositions
  
  subroutine SetChart(Chart, Length, Dimensionality, DOF, Cell, Face, Values, ValueIndex, UseFunction, Positions)
    ! Set the chart according to given data.
    implicit none
    type (BC_Chart),    intent(  OUT)    :: Chart
    integer( int_kind),             intent(IN   ) :: Length, Dimensionality, DOF, Cell
    integer( int_kind), OPTIONAL,   intent(IN   ) :: Face
    real(real_kind), dimension(:,:),intent(IN   ) :: Values
    integer( int_kind),dimension(:),intent(IN   ) :: ValueIndex
    logical( log_kind),dimension(:),intent(IN   ) :: UseFunction
    real(real_kind), dimension(:,:),intent(IN   ) :: Positions

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

    return
  end subroutine SetChart

  function ChartMatches(Chart, Cell, Face)
    implicit none
    type (BC_Chart), intent(IN) :: Chart
    integer( int_kind),         intent(IN) :: Cell
    integer( int_kind),         intent(IN), &
                     OPTIONAL   :: Face
    logical :: ChartMatches 
    
    if (PRESENT(Face)) then
       ChartMatches = (Cell == BC_Chart_Cell(Chart)) &
            .AND.     (Face == BC_Chart_Face(Chart))
    else
       ChartMatches = (Cell == BC_Chart_Cell(Chart))
    end if
    RETURN
  end function ChartMatches

END Module BC_CHARTS

