!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_ATLASES_DATA_TYPES
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition atlas
  !   data types.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  use bc_atlases_internal
  use bc_enum_types
  Implicit None
  Private

  ! Stuff inherited from bc_atlases_internal
  PUBLIC :: BC_Atlas_Data
  PUBLIC :: BC_Atlas_Spec

  ! Stuff defined or extended here.
  PUBLIC :: BC_Atlas
  PUBLIC :: BC_Get_Face
  PUBLIC :: BC_Get_Cell
  PUBLIC :: BC_Get_Offset
  PUBLIC :: BC_Get_Length
  PUBLIC :: BC_Get_Values
  PUBLIC :: BC_Get_ValueIndex
  PUBLIC :: BC_Get_UseFunction
  PUBLIC :: BC_Get_Positions
  PUBLIC :: BC_Get_Data
  PUBLIC :: BC_Set_Data_Size
  PUBLIC :: BC_NumberOfCharts
  PUBLIC :: BC_Get_Name
  PUBLIC :: BC_Set_Name
  PUBLIC :: GET_SCOPE
  PUBLIC :: SET_SCOPE
  PUBLIC :: BC_Get_DOF
  PUBLIC :: BC_Set_DOF
  PUBLIC :: DATA_SIZE
  PUBLIC :: SIZE
  PUBLIC :: PERMUTE, CLONE, REDISTRIBUTE
  PUBLIC :: Dimensionality
  PUBLIC :: INITIALIZE, FREE, ALLOC, REALLOC, APPEND
  PUBLIC :: ASSIGNMENT(=)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! An Atlas
  type BC_Atlas
     PRIVATE
     ! These fields cannot change unless the atlas is re-allocated
     type(BC_Atlas_Max_Extents) :: MaxSizes
     ! These fields change depending on data and use of the atlas
     type(BC_Atlas_Current) :: Current

     type(BC_Atlas_Spec) :: Spec
     type(BC_Atlas_Data) :: Data

     CHARACTER(BC_STRING_LEN) :: Name = 'NO NAME'
  end type BC_ATLAS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE InitializeAtlas
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE AllocAtlas
  END INTERFACE

  INTERFACE FREE
     MODULE PROCEDURE FreeAtlas
  END INTERFACE

  INTERFACE REALLOC
     MODULE PROCEDURE ReAllocAtlas
  END INTERFACE

  INTERFACE BC_NumberOfCharts
     MODULE PROCEDURE AtlasNumberOfCharts
  END INTERFACE
  
  INTERFACE BC_MaxNumberOfCharts
     MODULE PROCEDURE AtlasMaxNumberOfCharts
  END INTERFACE
  
  INTERFACE BC_MaxAtlasDataSize
     MODULE PROCEDURE AtlasMaxAtlasDataSize
  END INTERFACE
  
  INTERFACE DIMENSIONALITY
     MODULE PROCEDURE AtlasDimensionality
  END INTERFACE

  INTERFACE BC_Get_Face
     MODULE PROCEDURE AtlasFaceFromIndex
     MODULE PROCEDURE AtlasFaceFromAtlas
  END INTERFACE
  
  INTERFACE BC_Get_Cell
     MODULE PROCEDURE AtlasCellFromIndex
     MODULE PROCEDURE AtlasCellFromAtlas
  END INTERFACE

  INTERFACE BC_Get_Offset
     MODULE PROCEDURE AtlasOffsetFromIndex
     MODULE PROCEDURE AtlasOffsetFromAtlas
  END INTERFACE

  INTERFACE BC_Get_Length
     MODULE PROCEDURE AtlasLengthFromIndex
     MODULE PROCEDURE AtlasLengthFromAtlas
  END INTERFACE

  INTERFACE BC_Get_Values
     MODULE PROCEDURE AtlasValuesFromIndex
     MODULE PROCEDURE AtlasValuesFromAtlas
  END INTERFACE
  
  INTERFACE BC_Get_ValueIndex
     MODULE PROCEDURE AtlasValueIndexFromIndex
     MODULE PROCEDURE AtlasValueIndexFromAtlas
  END INTERFACE
  
  INTERFACE BC_Get_UseFunction
     MODULE PROCEDURE AtlasUseFunctionFromIndex
     MODULE PROCEDURE AtlasUseFunctionFromAtlas
  END INTERFACE
  
  INTERFACE BC_Get_Positions
     MODULE PROCEDURE AtlasPositionsFromIndex
     MODULE PROCEDURE AtlasPositionsFromAtlas
  END INTERFACE
  
  INTERFACE BC_Get_Data
     MODULE PROCEDURE AtlasData
  END INTERFACE
  
  INTERFACE Get_SCOPE
     MODULE PROCEDURE AtlasGetScope
  END INTERFACE

  INTERFACE Set_SCOPE
     MODULE PROCEDURE AtlasSetScope
  END INTERFACE

  INTERFACE BC_Get_DOF
     MODULE PROCEDURE AtlasGetDOF
  END INTERFACE
  INTERFACE BC_Set_DOF
     MODULE PROCEDURE AtlasSetDOF
  END INTERFACE

  INTERFACE SIZE
     MODULE PROCEDURE AtlasAtlasDataSize
  END INTERFACE

  INTERFACE DATA_SIZE
     MODULE PROCEDURE AtlasAtlasDataSize
  END INTERFACE

  INTERFACE BC_Set_Data_Size
     MODULE PROCEDURE AtlasSetAtlasDataSize
  END INTERFACE
  
  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE AtlasCopy
  END INTERFACE

  INTERFACE APPEND
     MODULE PROCEDURE AtlasAppend
  END INTERFACE

  INTERFACE CLONE
     MODULE PROCEDURE AtlasClone
  END INTERFACE

  INTERFACE BC_Get_Name
     MODULE PROCEDURE AtlasGetName
  END INTERFACE
  INTERFACE BC_Set_Name
     MODULE PROCEDURE AtlasSetName
  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routines to create and destroy Atlas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine InitializeAtlas(Atlas)
    ! Set the size of an atlas to -1, do not allocate any memory.
    type(BC_Atlas), intent(INOUT) :: Atlas

    call INITIALIZE(Atlas%MaxSizes)
    call INITIALIZE(Atlas%Current)
    call INITIALIZE(Atlas%Spec)
    call INITIALIZE(Atlas%Data)
    call BC_Set_Name(Atlas, BC_NO_NAME)

  end subroutine InitializeAtlas

  subroutine AllocAtlas(Atlas, MaxNumberOfCharts, MaxAtlasDataSize, Dimensionality, DOF)
    ! Alloc memory for chart data and specifiers
    type(BC_Atlas), intent(INOUT) :: Atlas
    integer, intent(IN) :: MaxNumberOfCharts, MaxAtlasDataSize
    integer, intent(IN) :: DIMENSIONALITY, DOF

    call ALLOC(Atlas%Spec, SIZE_=MaxNumberOfCharts)
    call ALLOC(Atlas%Data, SIZE_=MaxAtlasDataSize, DIMENSIONALITY=DIMENSIONALITY, DOF = DOF)

    Atlas%MaxSizes%MaxNumberOfCharts = MaxNumberOfCharts
    Atlas%MaxSizes%MaxAtlasDataSize  = MaxAtlasDataSize
    Atlas%MaxSizes%Dimensionality    = DIMENSIONALITY
    call BC_Set_DOF(Atlas%MaxSizes, DOF)

    Atlas%Current%AtlasDataSize     = 0
    Atlas%Current%NumberOfCharts    = 0
    Atlas%Current%Scope             = BC_Cells_No_Scope

  end subroutine AllocAtlas

  subroutine FreeAtlas(Atlas)
    ! Free memory for chart data and specifiers
    type(BC_Atlas), intent(INOUT) :: Atlas
    call FREE(Atlas%Spec)
    call FREE(Atlas%Data)
    call INITIALIZE(Atlas)
  end subroutine FreeAtlas

  subroutine AtlasCopy(AtlasDest, AtlasSrc)
    ! Assign AtlasDest contents of AtlasSrc.  Okay if AtlasDest is larger.  This
    ! copies AtlasSrc into low portion of AtlasDest.
    ! This does not do any memory allocation.  Assumes that is already done.
    ! Since memory allocation has been done, and associated fields filed in,
    ! AtlasDest must be INOUT.
    type(BC_Atlas), intent(INOUT) :: AtlasDest
    type(BC_Atlas), intent(IN) :: AtlasSrc

    ! Copy the specs
    AtlasDest%Spec = AtlasSrc%Spec

    ! Copy the data
    AtlasDest%Data = AtlasSrc%Data

    ! Fill in the fields which may change 
    AtlasDest%Current =  AtlasSrc%Current

    AtlasDest%Name = BC_Get_Name(AtlasSrc)

  end subroutine AtlasCopy

  subroutine AtlasClone(AtlasDest, AtlasSrc)
    ! Given an empty AtlasDest, and a non-empty AtlasSrc,
    ! CLONE allocateds Dest to be the same size as Src, then copies
    ! all Src data into it.
    ! The cloned Atlas has all of its own memory, and is fully
    ! indepenent of the src.  
    type(BC_Atlas), intent(OUT) :: AtlasDest
    type(BC_Atlas), intent(IN) :: AtlasSrc

    ! First we have to initialize the dest
    call INITIALIZE(AtlasDest)

    ! Now make Dest the same size as Src
    call ALLOC(AtlasDest, MaxNumberOfCharts   = BC_MaxNumberOfCharts(AtlasSrc), &
                          MaxAtlasDataSize    = BC_MaxAtlasDataSize(AtlasSrc),  &
                          Dimensionality      = DIMENSIONALITY(AtlasSrc),       &
                          DOF                 = BC_Get_DOF(AtlasSrc))
    
    ! Finally, just copy the data from Src into Dest
    AtlasDest = AtlasSrc

  end subroutine AtlasClone

  subroutine ReAllocAtlas(Atlas, MaxNumberOfCharts, MaxAtlasDataSize)
    ! Resize Atlas to hold more charts, and more chart data space.
    ! Returned Atlas has original data at head of all lists
    type(BC_Atlas), intent(INOUT) :: Atlas
    integer, intent(IN) :: MaxNumberOfCharts, MaxAtlasDataSize

    ! Local variables
    type(BC_Atlas) :: OrigAtlas

    ! Make space to hold the original data
    call INITIALIZE(OrigAtlas)
    call ALLOC(OrigAtlas, MaxNumberOfCharts   = BC_NumberOfCharts(Atlas), &
                          MaxAtlasDataSize    = SIZE(Atlas),              &
                          DIMENSIONALITY      = DIMENSIONALITY(Atlas),    &
                          DOF                 = BC_Get_DOF(Atlas))

    ! Copy original data into hold area
    OrigAtlas = Atlas

    ! Destroy input.
    call FREE(Atlas)
    
    ! Create new space
    call ALLOC(Atlas, MaxNumberOfCharts = MaxNumberOfCharts, &
                      MaxAtlasDataSize  = MaxAtlasDataSize,   &
                      DIMENSIONALITY    = DIMENSIONALITY(OrigAtlas), &
                      DOF               = BC_Get_DOF(OrigAtlas))

    ! Copy data to new structure
    Atlas = OrigAtlas

    ! Destroy the temp
    call FREE(OrigAtlas)

  end subroutine ReAllocAtlas
  
  subroutine AtlasAppend(Atlas, Values, Positions, ValueIndex, UseFunction, Cell, Face)
    ! Append a chart of data onto an Atlas
    type(BC_Atlas),       intent(INOUT) :: Atlas
    real(r8), dimension(:,:), intent(IN) :: Values
    real(r8), dimension(:,:), intent(IN) :: Positions
    integer,dimension(:), intent(IN) :: ValueIndex
    logical,dimension(:), intent(IN) :: UseFunction
    integer, intent(IN) :: Cell
    integer, OPTIONAL, intent(IN) :: Face

    ! Local variables
    integer :: Offset, Length, AtlasIndex, NewAtlasDataSize
    integer :: NewMaxAtlasDataSize, NewMaxNumberOfCharts
    
    ! The offset is just one more than the current chart data size
    Offset = SIZE(Atlas) + 1
    
    ! The length is the amount of data to append
    Length = SIZE(Values,2)

    ! The new chart size is the old plus the additional
    NewAtlasDataSize = SIZE(Atlas) + Length

    ! If we don't have enough chart space, then we need to make more
    if (NewAtlasDataSize > BC_MaxAtlasDataSize(Atlas)) then
       ! Grow by doubling size
       NewMaxAtlasDataSize = 2 * BC_MaxAtlasDataSize(Atlas)
       call REALLOC(Atlas%Data, SIZE_=NewMaxAtlasDataSize)
       Atlas%MaxSizes%MaxAtlasDataSize = NewMaxAtlasDataSize
    end if

    ! The AtlasIndex is just the current number of charts + 1
    AtlasIndex = BC_NumberOfCharts(Atlas) + 1

    ! If we don't have enough chart space we need to make more
    if (AtlasIndex > BC_MaxNumberOfCharts(Atlas)) then
       ! Grow by doubling
       NewMaxNumberOfCharts = 2*BC_MaxNumberOfCharts(Atlas)
       call REALLOC(Atlas%Spec, SIZE_=NewMaxNumberOfCharts)
       Atlas%MaxSizes%MaxNumberOfCharts = 2* NewMaxNumberOfCharts
    end if
    
    ! Append the spec info
    call APPEND(Atlas%Spec, AtlasIndex, LENGTH=Length, OFFSET=OFFSET)
    Atlas%Current%NumberOfCharts = AtlasIndex
    ! Append the data
    call APPEND(Atlas%Data, OFFSET=Offset,              &
                            VALUES=Values,              &
                            POSITIONS=Positions,        &
                            VALUEINDEX=ValueIndex,      &
                            USEFUNCTION=UseFunction,    &
                            CELL=Cell,                  &
                            FACE=Face)
    Atlas%Current%AtlasDataSize = NewAtlasDataSize
  end subroutine AtlasAppend
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routines to access Atlas fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function AtlasData(Atlas)
    type(BC_Atlas), TARGET, intent(IN):: Atlas
    type(BC_Atlas_Data), POINTER :: AtlasData
    AtlasData => Atlas%Data
  end function AtlasData

  function AtlasNumberOfCharts(Atlas)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: AtlasNumberOfCharts
    AtlasNumberOfCharts = Atlas%Current%NumberOfCharts
  end function AtlasNumberOfCharts

  function AtlasMaxNumberOfCharts(Atlas)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: AtlasMaxNumberOfCharts
    AtlasMaxNumberOfCharts = Atlas%MaxSizes%MaxNumberOfCharts
  end function AtlasMaxNumberOfCharts

  function AtlasDimensionality(Atlas)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: AtlasDimensionality
    AtlasDimensionality = Atlas%MaxSizes%Dimensionality
  end function AtlasDimensionality

  subroutine AtlasSetDOF(Atlas, DOF)
    type(BC_Atlas), intent(INOUT) :: Atlas
    integer, intent(IN) :: DOF
    call BC_Set_DOF(Atlas%MaxSizes, DOF)
  end subroutine AtlasSetDOF

  function AtlasGetDOF(Atlas) RESULT(DOF)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: DOF
    DOF = BC_Get_DOF(Atlas%MaxSizes)
  end function AtlasGetDOF

  function AtlasFaceFromIndex(Atlas, AtlasIndex) RESULT(Face)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    integer :: Face
    Face = Atlas%Data%Face(BC_Get_Offset(Atlas, AtlasIndex))
  end function AtlasFaceFromIndex

  function AtlasFaceFromAtlas(Atlas) RESULT(Faces)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, dimension(:), POINTER:: Faces
    Faces => Atlas%Data%Face
  end function AtlasFaceFromAtlas

  function AtlasCellFromIndex(Atlas, AtlasIndex) RESULT(Cell)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    integer :: Cell
    Cell = Atlas%Data%Cell(BC_Get_Offset(Atlas, AtlasIndex))
  end function AtlasCellFromIndex

  function AtlasCellFromAtlas(Atlas) RESULT(Cells)
    type(BC_Atlas), intent(IN)   :: Atlas
    integer, dimension(:), POINTER:: Cells
    Cells => Atlas%Data%Cell
  end function AtlasCellFromAtlas

  function AtlasOffsetFromIndex(Atlas, AtlasIndex) RESULT(Offset)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    integer :: Offset
    Offset = Atlas%Spec%Offset(AtlasIndex)
  end function AtlasOffsetFromIndex

  function AtlasOffsetFromAtlas(Atlas) RESULT(Offsets)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, dimension(:), POINTER:: Offsets
    Offsets => Atlas%Spec%Offset
  end function AtlasOffsetFromAtlas

  function AtlasLengthFromIndex(Atlas, AtlasIndex) RESULT(Length)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    integer :: Length
    Length = Atlas%Spec%Lengths(AtlasIndex)
  end function AtlasLengthFromIndex

  function AtlasLengthFromAtlas(Atlas) RESULT(Lengths)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, dimension(:), POINTER:: Lengths
    Lengths => Atlas%Spec%Lengths
  end function AtlasLengthFromAtlas

  function AtlasMaxAtlasDataSize(Atlas)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: AtlasMaxAtlasDataSize
    AtlasMaxAtlasDataSize = SIZE(Atlas%Data)
  end function AtlasMaxAtlasDataSize

  function AtlasAtlasDataSize(Atlas)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: AtlasAtlasDataSize
    AtlasAtlasDataSize = Atlas%Current%AtlasDataSize
  end function AtlasAtlasDataSize

  subroutine AtlasSetAtlasDataSize(Atlas, SIZE)
    type(BC_Atlas), intent(INOUT) :: Atlas
    integer :: SIZE
    Atlas%Current%AtlasDataSize = SIZE
  end subroutine AtlasSetAtlasDataSize

  function AtlasValuesFromIndex(Atlas, AtlasIndex) RESULT(Values)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    real(r8), dimension(:,:), POINTER :: Values

    ! Local variables
    integer :: Offset, Length

    Offset = BC_Get_Offset(Atlas, AtlasIndex)
    Length = BC_Get_Length(Atlas, AtlasIndex)
    
    Values => Atlas%Data%Values(:,Offset: Offset+Length-1)

  end function AtlasValuesFromIndex

  function AtlasValuesFromAtlas(Atlas) RESULT(Values)
    type(BC_Atlas), intent(IN) :: Atlas
    real(r8), dimension(:,:), POINTER :: Values

    Values => Atlas%Data%Values

  end function AtlasValuesFromAtlas

  function AtlasValueIndexFromIndex(Atlas, AtlasIndex) RESULT(ValueIndex)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    integer, dimension(:), POINTER:: ValueIndex

    ! Local variables
    integer :: Offset, Length

    Offset = BC_Get_Offset(Atlas, AtlasIndex)
    Length = BC_Get_Length(Atlas, AtlasIndex)
    
    ValueIndex => Atlas%Data%ValueIndex(Offset: Offset+Length-1)

  end function AtlasValueIndexFromIndex

  function AtlasValueIndexFromAtlas(Atlas) RESULT(ValueIndex)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, dimension(:), POINTER:: ValueIndex

    ValueIndex => Atlas%Data%ValueIndex

  end function AtlasValueIndexFromAtlas

  function AtlasUseFunctionFromIndex(Atlas, AtlasIndex) RESULT(UseFunction)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    logical, dimension(:), POINTER:: UseFunction

    ! Local variables
    integer :: Offset, Length

    Offset = BC_Get_Offset(Atlas, AtlasIndex)
    Length = BC_Get_Length(Atlas, AtlasIndex)
    
    UseFunction => Atlas%Data%UseFunction(Offset: Offset+Length-1)

  end function AtlasUseFunctionFromIndex

  function AtlasUseFunctionFromAtlas(Atlas) RESULT(UseFunction)
    type(BC_Atlas), intent(IN) :: Atlas
    logical, dimension(:), POINTER:: UseFunction
    UseFunction => Atlas%Data%UseFunction
  end function AtlasUseFunctionFromAtlas

  function AtlasPositionsFromIndex(Atlas, AtlasIndex) RESULT(Positions)
    type(BC_Atlas), intent(IN) :: Atlas
    integer, intent(IN) :: AtlasIndex
    real(r8),dimension(:,:), POINTER  :: Positions

    ! Local variables
    integer :: Offset, Length

    Offset = BC_Get_Offset(Atlas, AtlasIndex)
    Length = BC_Get_Length(Atlas, AtlasIndex)
    
    Positions => Atlas%Data%Positions(:, Offset: Offset+Length-1)

  end function AtlasPositionsFromIndex

  function AtlasPositionsFromAtlas(Atlas) RESULT(Positions)
    type(BC_Atlas), intent(IN) :: Atlas
    real(r8),dimension(:,:), POINTER  :: Positions
    Positions => Atlas%Data%Positions
  end function AtlasPositionsFromAtlas

  function AtlasGetScope(Atlas) RESULT(SCOPE)
    type(BC_Atlas), intent(IN) :: Atlas
    integer :: SCOPE
    SCOPE = Atlas%Current%Scope
  end function AtlasGetScope

  subroutine AtlasSetScope(Atlas, SCOPE)
    type(BC_Atlas), intent(INOUT) :: Atlas
    integer, intent(IN) :: SCOPE
    Atlas%Current%Scope = SCOPE
  end subroutine AtlasSetScope
    

  subroutine AtlasSetName(Atlas, Name)
    ! Set the name field in the Atlas
    type(BC_Atlas), intent(INOUT) :: Atlas
    character(*), INTENT(IN) :: Name
    Atlas%Name = Name
  end subroutine AtlasSetName
  
  function AtlasGetName(Atlas) RESULT(Name)
    ! Get the name field in the Atlas
    type(BC_Atlas), intent(IN) :: Atlas
    character(BC_STRING_LEN) :: Name
    Name = TRIM(Atlas%Name)
  end function AtlasGetName
    

END Module BC_ATLASES_DATA_TYPES

