Module BC_ATLASES_INTERNAL
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition atlas
  !   operations
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

  PUBLIC :: BC_Atlas_Max_Extents
  PUBLIC :: BC_Atlas_Current
  PUBLIC :: BC_Atlas_Data
  PUBLIC :: BC_Atlas_Spec
  PUBLIC :: SIZE
  PUBLIC :: INITIALIZE
  PUBLIC :: ALLOC, REALLOC, FREE, APPEND
  PUBLIC :: CLONE
  PUBLIC :: PERMUTE
  PUBLIC :: REDISTRIBUTE
  PUBLIC :: DIMENSIONALITY
  PUBLIC :: ASSIGNMENT(=)
  PUBLIC :: BC_Get_DOF
  PUBLIC :: BC_Set_DOF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Parameters which determine how much an Atlas can hold
  type BC_Atlas_Max_Extents
     ! These fields cannot change unless the atlas is re-allocated
     integer :: MaxNumberOfCharts
     integer :: MaxAtlasDataSize
     integer :: Dimensionality
     integer :: DegreesOfFreedom
  end type BC_Atlas_Max_Extents

  ! Parameters which are changed depending on what is in the Atlas
  type BC_Atlas_Current
     ! These fields change depending on data and use of the atlas
     integer :: NumberOfCharts
     integer :: AtlasDataSize
     integer :: Scope
  end type BC_Atlas_Current

  ! Type to hold Atlas data
  type BC_Atlas_Data
     integer( int_kind), pointer, dimension(:)   :: Face => null()
     integer( int_kind), pointer, dimension(:)   :: Cell => null()
     real(real_kind),    pointer, dimension(:,:) :: Values => null()
     integer( int_kind), pointer, dimension(:)   :: ValueIndex => null()
     logical(log_kind),  pointer, dimension(:)   :: UseFunction => null()
     real(real_kind),    pointer, dimension(:,:) :: Positions => null()
  end type BC_Atlas_Data

  ! Type to identify charts in an atlas
  type BC_Atlas_Spec
     integer( int_kind), pointer, dimension(:) :: Offset => null()
     integer( int_kind), pointer, dimension(:) :: Lengths => null()
  end type BC_Atlas_Spec

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE InitializeMaxExtents
     MODULE PROCEDURE InitializeCurrent
     MODULE PROCEDURE InitializeAtlasSpec
     MODULE PROCEDURE InitializeAtlasData
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE AllocAtlasSpec
     MODULE PROCEDURE AllocAtlasData
  END INTERFACE

  INTERFACE FREE
     MODULE PROCEDURE FreeAtlasSpec
     MODULE PROCEDURE FreeAtlasData
  END INTERFACE

  INTERFACE DIMENSIONALITY
     MODULE PROCEDURE AtlasDataDimensionality
  END INTERFACE
  
  INTERFACE SIZE
     MODULE PROCEDURE AtlasSpecSize
     MODULE PROCEDURE AtlasDataSize
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE AtlasSpecCopy
     MODULE PROCEDURE AtlasDataCopy
  END INTERFACE

  INTERFACE CLONE
     MODULE PROCEDURE AtlasSpecClone
     MODULE PROCEDURE AtlasDataClone
  END INTERFACE
  
  INTERFACE REALLOC
     MODULE PROCEDURE AtlasSpecRealloc
     MODULE PROCEDURE AtlasDataRealloc
  END INTERFACE
  
  INTERFACE APPEND
     MODULE PROCEDURE AtlasSpecAppend
     MODULE PROCEDURE AtlasDataAppend
  END INTERFACE

  INTERFACE PERMUTE
     MODULE PROCEDURE AtlasDataPermute
  END INTERFACE

  INTERFACE REDISTRIBUTE
     MODULE PROCEDURE AtlasDataRedistribute
  END INTERFACE
  
  INTERFACE BC_Get_DOF
     MODULE PROCEDURE AtlasMaxGetDOF
     MODULE PROCEDURE AtlasDataGetDOF
  END INTERFACE
  INTERFACE BC_Set_DOF
     MODULE PROCEDURE AtlasMaxSetDOF
  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routines to handle the Max_Extents and Current types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine InitializeMaxExtents(AtlasMax)
    implicit none
    type (BC_Atlas_Max_Extents), intent(INOUT) :: AtlasMax

    AtlasMax%MaxNumberOfCharts = BC_INVALID_COUNT
    AtlasMax%MaxAtlasDataSize  = BC_INVALID_COUNT
    AtlasMax%Dimensionality    = BC_INVALID_COUNT
    AtlasMax%DegreesOfFreedom  = BC_INVALID_COUNT
  end subroutine InitializeMaxExtents

  subroutine InitializeCurrent(AtlasCurrent)
    implicit none
    type (BC_Atlas_Current), intent(INOUT) :: AtlasCurrent
    AtlasCurrent%NumberOfCharts = BC_INVALID_COUNT
    AtlasCurrent%AtlasDataSize  = BC_INVALID_COUNT
    AtlasCurrent%Scope          = BC_Cells_No_Scope
  end subroutine InitializeCurrent

  subroutine AtlasMaxSetDOF(AtlasMax, DOF)
    implicit none
    type (BC_Atlas_Max_Extents), intent(INOUT) :: AtlasMax
    integer, intent(IN) :: DOF

    AtlasMax%DegreesOfFreedom  = DOF
    return
  end subroutine AtlasMaxSetDOF

  function AtlasMaxGetDOF(AtlasMax) RESULT(DOF)
    implicit none
    type (BC_Atlas_Max_Extents), intent(IN) :: AtlasMax
    integer :: DOF
    DOF = AtlasMax%DegreesOfFreedom
    return
  end function AtlasMaxGetDOF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routines to create and destroy Atlas_Spec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine InitializeAtlasSpec(AtlasSpec)
    implicit none
    type (BC_Atlas_Spec), intent(INOUT) :: AtlasSpec
    NULLIFY(AtlasSpec%Offset)
    NULLIFY(AtlasSpec%Lengths)
    return
  end subroutine InitializeAtlasSpec

  subroutine AllocAtlasSpec(AtlasSpec, SIZE)
    implicit none
    type (BC_Atlas_Spec), intent(INOUT) :: AtlasSpec
    integer( int_kind),              intent(IN   ) :: SIZE

    ALLOCATE(AtlasSpec%Offset(SIZE))
    ALLOCATE(AtlasSpec%Lengths(SIZE))
    
    AtlasSpec%Offset = BC_INVALID_OFFSET
    AtlasSpec%Lengths = BC_INVALID_LENGTH
    
    return
  end subroutine AllocAtlasSpec
  
  subroutine FreeAtlasSpec(AtlasSpec)
    implicit none
    type (BC_Atlas_Spec), intent(INOUT) :: AtlasSpec
    
    if (ASSOCIATED(AtlasSpec%Offset))  DEALLOCATE(AtlasSpec%Offset)
    if (ASSOCIATED(AtlasSpec%Lengths)) DEALLOCATE(AtlasSpec%Lengths)
    return
  end subroutine FreeAtlasSpec

  function AtlasSpecSize(AtlasSpec)
    implicit none
    type(BC_Atlas_Spec), intent(IN) :: AtlasSpec
    integer( int_kind) :: AtlasSpecSize
    AtlasSpecSize = SIZE(AtlasSpec%Offset)
    return
  END function AtlasSpecSize
    
  subroutine AtlasSpecCopy(AtlasSpecDest, AtlasSpecSrc)
    ! Assign AtlasSpecDest contents of AtlasSpecSrc.  Okay if AtlasSpecDest is larger.  This
    ! copies AtlasSpecSrc into low portion of AtlasSpecDest.
    ! This does not do any memory allocation.  Assumes that is already done.
    implicit none
    type(BC_Atlas_Spec), intent(INOUT) :: AtlasSpecDest
    type(BC_Atlas_Spec), intent(IN   ) :: AtlasSpecSrc

    ! Local variables
    integer :: DestMaxNumCharts, SrcNumCharts

    ! Check that dest has enough chart space
    SrcNumCharts     = SIZE(AtlasSpecSrc)
    DestMaxNumCharts = SIZE(AtlasSpecDest)
    if (DestMaxNumCharts < SrcNumCharts) then
       print *, 'MaxNumberOfCharts not large enough for LHS in AtlasSpecCopy'
    end if

    ! Assume that dest is already initialized as needed.
    AtlasSpecDest%Offset(1:SrcNumCharts) = AtlasSpecSrc%Offset(1:SrcNumCharts)
    AtlasSpecDest%Lengths(1:SrcNumCharts) = AtlasSpecSrc%Lengths(1:SrcNumCharts)
    return
  end subroutine AtlasSpecCopy
  
  subroutine AtlasSpecClone(AtlasSpecDest, AtlasSpecSrc)
    ! Given an empty AtlasSpecDest, and a non-empty AtlasSpecSrc,
    ! CLONE allocateds Dest to be the same size as Src, then copies
    ! all Src data into it.
    ! The cloned AtlasSpec has all of its own memory, and is fully
    ! indepenent of the src.  
    implicit none
    type(BC_Atlas_Spec), intent(  OUT) :: AtlasSpecDest
    type(BC_Atlas_Spec), intent(IN   ) :: AtlasSpecSrc

    ! Before we can do anything with Dest, we need to initialize it.
    call INITIALIZE(AtlasSpecDest)

    ! Now make Dest have that size
    call ALLOC(AtlasSpecDest, SIZE=SIZE(AtlasSpecSrc))
    
    ! Finally, just copy the data from Src into Dest
    AtlasSpecDest = AtlasSpecSrc
    return
  end subroutine AtlasSpecClone
  
  subroutine AtlasSpecRealloc(AtlasSpec, SIZE)
    ! Reallocate AtlasSpec to the new size.  If the new size is larger
    ! than the current size, then the new spec contains the original
    ! data at the beginning of it's lists.  
    type (BC_Atlas_Spec), intent(INOUT) :: AtlasSpec
    integer( int_kind),              intent(IN   ) :: SIZE

    ! Local variables
    type(BC_Atlas_Spec) :: Orig_Spec

    ! Make a copy of the original
    call INITIALIZE(Orig_Spec)
    call ALLOC(Orig_Spec, SIZE=AtlasSpecSize(AtlasSpec))
    Orig_Spec = AtlasSpec

    ! Now get rid of the input
    call FREE(AtlasSpec)

    ! Re-Create in the new size
    call ALLOC(AtlasSpec, SIZE=SIZE)
    
    ! Set the fields to BC_INVALID
    AtlasSpec%Lengths = BC_INVALID_LENGTH
    AtlasSpec%Offset = BC_INVALID_OFFSET

    ! Copy original data into new space
    AtlasSpec = Orig_Spec

    ! Free up the temporary
    call FREE(Orig_Spec)

    return
  end subroutine AtlasSpecRealloc

  subroutine AtlasSpecAppend(AtlasSpec, AtlasIndex, LENGTH, OFFSET)
    ! Append the data to the Spec.  Assumes that there is enough room.
    ! Does not do any reallocation of memory.
    implicit none
    type(BC_Atlas_Spec), intent(INOUT) :: AtlasSpec
    integer( int_kind),             intent(IN   ) :: AtlasIndex
    integer( int_kind),             intent(IN   ) :: LENGTH, OFFSET

    AtlasSpec%LENGTHS(AtlasIndex) = LENGTH
    AtlasSpec%OFFSET(AtlasIndex) = OFFSET
    return
  end subroutine AtlasSpecAppend


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routines to create and destroy Atlas_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine InitializeAtlasDATA(AtlasDATA)
    implicit none
    type (BC_Atlas_Data), intent(INOUT) :: AtlasDATA
    NULLIFY(AtlasData%Face)
    NULLIFY(AtlasData%Cell)
    NULLIFY(AtlasData%Values)
    NULLIFY(AtlasData%ValueIndex)
    NULLIFY(AtlasData%UseFunction)
    NULLIFY(AtlasData%Positions)
    return
  end subroutine InitializeAtlasDATA

  subroutine AllocAtlasDATA(AtlasDATA, SIZE, DIMENSIONALITY, DOF)
    implicit none
    type (BC_Atlas_Data), intent(INOUT) :: AtlasDATA
    integer( int_kind),              intent(IN   ) :: SIZE, DIMENSIONALITY, DOF

    ALLOCATE(AtlasData%Face(SIZE))
    ALLOCATE(AtlasData%Cell(SIZE))
    ALLOCATE(AtlasData%Values(DOF, SIZE))
    ALLOCATE(AtlasData%ValueIndex(SIZE))
    ALLOCATE(AtlasData%UseFunction(SIZE))
    ALLOCATE(AtlasData%Positions(Dimensionality, SIZE))
    AtlasData%Face = BC_INVALID_FACE
    AtlasData%Cell = BC_INVALID_CELL
    AtlasData%ValueIndex = BC_INVALID_CELL
    AtlasData%UseFunction = .FALSE.

    return
  end subroutine AllocAtlasDATA

  subroutine FreeAtlasData(AtlasData)
    implicit none
    type (BC_Atlas_Data), intent(INOUT) :: AtlasData
    
    if (ASSOCIATED(AtlasData%Face))    DEALLOCATE(AtlasData%Face)
    if (ASSOCIATED(AtlasData%Cell))    DEALLOCATE(AtlasData%Cell)
    if (ASSOCIATED(AtlasData%Values))    DEALLOCATE(AtlasData%Values)
    if (ASSOCIATED(AtlasData%ValueIndex))    DEALLOCATE(AtlasData%ValueIndex)
    if (ASSOCIATED(AtlasData%UseFunction))    DEALLOCATE(AtlasData%UseFunction)
    if (ASSOCIATED(AtlasData%Positions)) DEALLOCATE(AtlasData%Positions)
    return
  end subroutine FreeAtlasData

  function AtlasDataSize(AtlasData)
    implicit none
    type(BC_Atlas_Data), intent(IN) :: AtlasData
    integer( int_kind) :: AtlasDataSize
    AtlasDataSize = SIZE(AtlasData%Values,2)
    return
  END function AtlasDataSize
    
  function AtlasDataDimensionality(AtlasData)
    implicit none
    type(BC_Atlas_Data), intent(IN) :: AtlasData
    integer( int_kind) :: AtlasDataDimensionality
    AtlasDataDimensionality = SIZE(AtlasData%Positions,1)
    return
  END function AtlasDataDimensionality
    
  function AtlasDataGetDOF(AtlasData) RESULT(DOF)
    implicit none
    type(BC_Atlas_Data), intent(IN) :: AtlasData
    integer :: DOF
    DOF = SIZE(AtlasData%Values,1)
    return
  end function AtlasDataGetDOF

  subroutine AtlasDataCopy(AtlasDataDest, AtlasDataSrc)
    ! Assign AtlasDataDest contents of AtlasDataSrc.  Okay if AtlasDataDest is larger.  This
    ! copies AtlasDataSrc into low portion of AtlasDataDest.
    ! This does not do any memory allocation.  Assumes that is already done.
    implicit none
    type(BC_Atlas_Data), intent(INOUT) :: AtlasDataDest
    type(BC_Atlas_Data), intent(IN   ) :: AtlasDataSrc

    ! Local variables
    integer :: DestChartDataSpace, SrcChartDataSpace
    integer :: SrcDim, DestDim, SrcDOF, DestDOF

    SrcDim  = Dimensionality(AtlasDataSrc)
    DestDim = Dimensionality(AtlasDataDest)
    SrcDOF  = BC_Get_DOF(AtlasDataSrc)
    DestDOF = BC_Get_DOF(AtlasDataDest)
    SrcChartDataSpace  = SIZE(AtlasDataSrc)
    DestChartDataSpace = SIZE(AtlasDataDest)

#ifdef BC_TESTS_ON
    ! Check that both have same dimensionality
    if (SrcDim /= DestDim) then
       print *, 'Dimensionality does not match in AtlasCopy'
    end if

    ! Check that both have same DOF
    if (SrcDOF /= DestDOF) then
       print *, 'Degrees of Freedom do not match in AtlasCopy.'
    end if

    ! Check that dest has enough chart space
    if (DestChartDataSpace < SrcChartDataSpace) then
       print *, 'ChartDataSpace not large enough in LHS of AtlasDataCopy'
    end if

#endif

    ! Assume that dest is already initialized as needed.
    AtlasDataDest%Face(1:SrcChartDataSpace) = AtlasDataSrc%Face(1:SrcChartDataSpace)
    AtlasDataDest%Cell(1:SrcChartDataSpace) = AtlasDataSrc%Cell(1:SrcChartDataSpace)
    AtlasDataDest%Values(:,1:SrcChartDataSpace) = AtlasDataSrc%Values(:,1:SrcChartDataSpace)
    AtlasDataDest%ValueIndex(1:SrcChartDataSpace) = AtlasDataSrc%ValueIndex(1:SrcChartDataSpace)
    AtlasDataDest%UseFunction(1:SrcChartDataSpace) = AtlasDataSrc%UseFunction(1:SrcChartDataSpace)
    AtlasDataDest%Positions(:, 1:SrcChartDataSpace) = AtlasDataSrc%Positions(:, 1:SrcChartDataSpace)
    return
    
  end subroutine AtlasDataCopy
  
  subroutine AtlasDataClone(AtlasDataDest, AtlasDataSrc)
    ! Given an empty AtlasDataDest, and a non-empty AtlasDataSrc,
    ! CLONE allocateds Dest to be the same size as Src, then copies
    ! all Src data into it.
    ! The cloned AtlasData has all of its own memory, and is fully
    ! indepenent of the src.  
    implicit none
    type(BC_Atlas_Data), intent(  OUT) :: AtlasDataDest
    type(BC_Atlas_Data), intent(IN   ) :: AtlasDataSrc

    ! Before we can do anything with Dest, we need to initialize it.
    call INITIALIZE(AtlasDataDest)

    ! Now make Dest of the size of the Src
    call ALLOC(AtlasDataDest, SIZE=SIZE(AtlasDataSrc) , &
                              DIMENSIONALITY=DIMENSIONALITY(AtlasDataSrc),&
                              DOF = BC_Get_DOF(AtlasDataSrc))
    
    ! Finally, just copy the data from Src into Dest
    AtlasDataDest = AtlasDataSrc
    return
  end subroutine AtlasDataClone
  
  subroutine AtlasDataRealloc(AtlasData, SIZE)
    ! Reallocate AtlasData to the new size.  If the new size is larger
    ! than the current size, then the new Data contains the original
    ! data at the beginning of it's lists.  
    type (BC_Atlas_Data), intent(INOUT) :: AtlasData
    integer( int_kind),              intent(IN   ) :: SIZE

    ! Local variables
    type(BC_Atlas_Data) :: Orig_Data

    ! Make a copy of the original
    call INITIALIZE(Orig_Data)
    call ALLOC(Orig_Data, SIZE=AtlasDataSize(AtlasData), &
                          DIMENSIONALITY = DIMENSIONALITY(AtlasData), &
                          DOF = BC_Get_DOF(AtlasData))
    Orig_Data = AtlasData

    ! Now get rid of the input
    call FREE(AtlasData)

    ! Re-Create in the new size
    call ALLOC(AtlasData, SIZE=SIZE, &
                          DIMENSIONALITY = DIMENSIONALITY(Orig_Data),&
                          DOF = BC_Get_DOF(Orig_Data))
    
    ! Set the fields to BC_INVALID
    AtlasData%Face = BC_INVALID_FACE
    AtlasData%Cell = BC_INVALID_CELL
    AtlasData%ValueIndex = BC_INVALID_CELL
    AtlasData%UseFunction = .FALSE.

    ! Copy original data into new space
    AtlasData = Orig_Data

    ! Free up the temporary
    call FREE(Orig_Data)

    return
  end subroutine AtlasDataRealloc

  subroutine AtlasDataAppend(AtlasData, OFFSET, VALUES, POSITIONS, USEFUNCTION, VALUEINDEX, CELL, FACE)
    ! Append the data to the Data.  Assumes that there is enough room.
    ! Does not do any reallocation of memory.
    implicit none
    type(BC_Atlas_Data), intent(INOUT) :: AtlasData
    integer( int_kind),             intent(IN   ) :: OFFSET
    real(real_kind), dimension(:,:),intent(IN   ) :: Values
    real(real_kind), dimension(:,:),intent(IN   ) :: Positions
    integer( int_kind),dimension(:),intent(IN   ) :: VALUEINDEX
    logical( log_kind),dimension(:),intent(IN   ) :: USEFUNCTION
    integer( int_kind),             intent(IN   ) :: CELL
    integer( int_kind), OPTIONAL,   intent(IN   ) :: FACE

    ! Local variables
    integer :: AppendSize
    AppendSize = SIZE(Values,2)

    ! Check that dimensionalities are same.
    if (DIMENSIONALITY(AtlasData) /= SIZE(POSITIONS,1)) then
       print *, 'ERROR: in AtlasDataAppend, dimensionalities different'
    end if
    
    ! Check that DOF are the same
    if (BC_Get_DOF(AtlasData) /= SIZE(Values,1)) then
       print *, 'ERROR: in AtlasDataAppend, Degrees of Freedom different.'
    endif

    AtlasData%Values(:,Offset:Offset+AppendSize-1) = Values(:,:)
    AtlasData%Positions(:,Offset:Offset+AppendSize-1) = Positions(:,:)
    AtlasData%CELL(Offset:Offset+AppendSize-1) = CELL
    AtlasData%ValueIndex(Offset:Offset+AppendSize-1) = VALUEINDEX
    AtlasData%UseFunction(Offset:Offset+AppendSize-1) = USEFUNCTION
    if (PRESENT(FACE)) then
       AtlasData%FACE(Offset:Offset+AppendSize-1) = FACE
    end if
    return
  end subroutine AtlasDataAppend

  subroutine AtlasDataPermute(AtlasData, Rank)
    ! Permute the atlas data according to the input vector RANK.
    ! 
    use pgslib_module, ONLY: PGSLib_GS_Trace, PGSLib_Permute, &
                             PGSLib_DEALLOCATE_TRACE
    implicit none
    type(BC_Atlas_Data),   TARGET,            intent(INOUT) :: AtlasData
    integer( int_kind), dimension(SIZE(AtlasData%Face)), intent(IN   ) :: Rank

    ! Local variables
    integer :: DataSize, d
    integer( int_kind),  dimension(SIZE(AtlasData%Face)  ) :: Cells_Temp
    integer( int_kind),  dimension(SIZE(AtlasData%Face)  ) :: Faces_Temp
    integer( int_kind),  dimension(SIZE(AtlasData%Face)  ) :: ValueIndex_Temp
    logical( log_kind),  dimension(SIZE(AtlasData%Face)  ) :: UseFunction_Temp
    real(real_kind),     dimension(SIZE(AtlasData%Values,1), &
                                   SIZE(AtlasData%Values,2)):: Values_Temp

    real(real_kind),     dimension(SIZE(AtlasData%Positions,1), &
                                   SIZE(AtlasData%Positions,2)) :: Positions_Temp
    type (PGSLib_GS_Trace), POINTER                     :: Permute_Trace

    integer( int_kind), POINTER, dimension(:  ) :: Cells
    integer( int_kind), POINTER, dimension(:  ) :: Faces
    integer( int_kind), POINTER, dimension(:  ) :: ValueIndex
    logical( log_kind), POINTER, dimension(:  ) :: UseFunction
    real(real_kind),    POINTER, dimension(:,:) :: Values
    real(real_kind),    POINTER, dimension(:,:) :: Positions
    
    DataSize = SIZE(AtlasData)

    Cells => AtlasData%Cell
    Faces => AtlasData%Face
    ValueIndex => AtlasData%ValueIndex
    UseFunction => AtlasData%UseFunction
    Values => AtlasData%Values
    Positions => AtlasData%Positions

    ! Permute data into temporary space, then overwrite
    NULLIFY(Permute_Trace)
    call PGSLib_Permute(DEST = Cells_Temp, &
                        SOURCE = Cells,    &
                        INDEX  = Rank,     &
                        TRACE  = Permute_Trace)
    Cells     = Cells_Temp
    call PGSLib_Permute(DEST = Faces_Temp, &
                        SOURCE = Faces,    &
                        INDEX  = Rank,     &
                        TRACE  = Permute_Trace)
    Faces     = Faces_Temp
    do d = 1, BC_Get_DOF(AtlasData)
       call PGSLib_Permute(DEST = Values_Temp(d,:), &
                           SOURCE = Values(d,:),    &
                           INDEX  = Rank,           &
                           TRACE  = Permute_Trace)
       Values(d,:)     = Values_Temp(d,:)
    end do
    call PGSLib_Permute(DEST = ValueIndex_Temp, &
                        SOURCE = ValueIndex,    &
                        INDEX  = Rank,     &
                        TRACE  = Permute_Trace)
    ValueIndex     = ValueIndex_Temp
    call PGSLib_Permute(DEST = UseFunction_Temp, &
                        SOURCE = UseFunction,    &
                        INDEX  = Rank,     &
                        TRACE  = Permute_Trace)
    UseFunction    = UseFunction_Temp
    do d = 1, DIMENSIONALITY(AtlasData)
       call PGSLib_Permute(DEST = Positions_Temp(d,:), &
                           SOURCE = Positions(d,:),    &
                           INDEX  = Rank,     &
                           TRACE  = Permute_Trace)
    end do

    Positions = Positions_Temp
    
    ! Done with trace
    call PGSLib_DEALLOCATE_TRACE(Permute_Trace)

    return
  end subroutine AtlasDataPermute

  subroutine AtlasDataRedistribute(Dest, Source)
    ! Redistribute the source atlas into the dest atlas
    ! Also set ChartDataSize to number of valid entries
    use pgslib_module,      ONLY: PGSLib_GS_Trace, PGSLib_Redistribute, &
                                  PGSLib_DEALLOCATE_TRACE
    implicit none
    type(BC_Atlas_Data),   TARGET,            intent(INOUT) :: Dest
    type(BC_Atlas_Data),   TARGET,            intent(IN   ) :: Source

    ! Local variables
    integer :: d
    type (PGSLib_GS_Trace), POINTER                     :: Redistribute_Trace

    integer( int_kind), POINTER, dimension(:  ) :: Cells_Dest, Cells_Source
    integer( int_kind), POINTER, dimension(:  ) :: Faces_Dest, Faces_Source
    integer( int_kind), POINTER, dimension(:  ) :: ValueIndex_Dest, ValueIndex_Source
    logical( log_kind), POINTER, dimension(:  ) :: UseFunction_Dest, UseFunction_Source
    real(real_kind),    POINTER, dimension(:,:) :: Values_Dest, Values_Source
    real(real_kind),    POINTER, dimension(:,:) :: Positions_Dest, Positions_Source
    
    Cells_Dest      => Dest%Cell
    Faces_Dest      => Dest%Face
    ValueIndex_Dest => Dest%ValueIndex
    UseFunction_Dest => Dest%UseFunction
    Values_Dest     => Dest%Values
    Positions_Dest  => Dest%Positions

    Cells_Source      => Source%Cell
    Faces_Source      => Source%Face
    ValueIndex_Source => Source%ValueIndex
    UseFunction_Source => Source%UseFunction
    Values_Source     => Source%Values
    Positions_Source  => Source%Positions

    ! Redistribute from source to dest
    NULLIFY(Redistribute_Trace)

    call PGSLib_Redistribute(DEST = Cells_Dest, SOURCE = Cells_Source, TRACE = Redistribute_Trace)
    call PGSLib_Redistribute(DEST = Faces_Dest, SOURCE = Faces_Source, TRACE = Redistribute_Trace)
    call PGSLib_Redistribute(DEST = ValueIndex_Dest, SOURCE = ValueIndex_Source, TRACE = Redistribute_Trace)
    call PGSLib_Redistribute(DEST = UseFunction_Dest, SOURCE = UseFunction_Source, TRACE = Redistribute_Trace)
    do d = 1, BC_Get_DOF(Dest)
       call PGSLib_Redistribute(DEST   = Values_Dest(d,:),   &
                                SOURCE = Values_Source(d,:), &
                                TRACE  = Redistribute_Trace)
    end do
    do d = 1, DIMENSIONALITY(Dest)
       call PGSLib_Redistribute(DEST = Positions_Dest(D,:),     &
                                SOURCE = Positions_Source(d,:), &
                                TRACE = Redistribute_Trace)
    end do

    ! Done with trace
    call PGSLib_DEALLOCATE_TRACE(Redistribute_Trace)

    return
  end subroutine AtlasDataRedistribute

END Module BC_ATLASES_INTERNAL

