!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_Regions
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition regions
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@lanl.gov)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  use bc_enum_types
  Implicit None
  Private

  PUBLIC :: BC_Region
  PUBLIC :: ASSIGNMENT(=), &
            ALLOC,         &
            APPEND,        &
            FREE,          &
            INITIALIZE,    &
            INSERT,        &
            REALLOC,       &
            SIZE,          &
            DIMENSIONALITY,&
            COLLATE

  PUBLIC :: ORDER
  PUBLIC :: PERMUTE
  PUBLIC :: CANONICALIZE

  PUBLIC :: BC_Region_Dimension, &
            BC_Region_Faces,     &
            BC_Region_OwnerCells,&
            BC_Region_Positions, &
            BC_Region_Values,    &
            BC_Region_UseFunction,&
            BC_Get_Name,         &
            BC_Set_Name,         &
            BC_Get_DOF,          &
            BC_Set_DOF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Type to define a region
  type BC_Region
     integer                           :: NumberOfEntries
     integer                           :: DegreesOfFreedom
     integer,  pointer, dimension(:)   :: OwnerCell => null()
     integer,  pointer, dimension(:)   :: Face => null()
     real(r8), pointer, dimension(:,:) :: Value => null()
     real(r8), pointer, dimension(:,:) :: Position => null()
     logical,  pointer, dimension(:)   :: UseFunction => null()
     character(BC_STRING_LEN)          :: NAME = BC_NO_NAME
  end type BC_Region
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE ALLOC
     MODULE PROCEDURE BC_Region_Alloc
  END INTERFACE
  
  INTERFACE APPEND
     MODULE PROCEDURE BC_Region_Append
  END INTERFACE
  
  INTERFACE FREE
     MODULE PROCEDURE BC_Region_Free
  END INTERFACE
  
  INTERFACE INITIALIZE
     MODULE PROCEDURE BC_Region_Initialize
  END INTERFACE
  
  INTERFACE INSERT
     MODULE PROCEDURE BC_Region_INSERT
  END INTERFACE

  INTERFACE REALLOC
     MODULE PROCEDURE BC_Region_REAlloc
  END INTERFACE

  INTERFACE SIZE
     MODULE PROCEDURE BC_Region_Size
  END INTERFACE

  INTERFACE DIMENSIONALITY
     MODULE PROCEDURE BC_Region_Dimension
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE BC_Region_Copy
  END INTERFACE
  
  INTERFACE CANONICALIZE
     MODULE PROCEDURE BC_Region_Canonical
  END INTERFACE

  INTERFACE ORDER
     MODULE PROCEDURE BC_Region_Order
  END INTERFACE

  INTERFACE COLLATE
     MODULE PROCEDURE BC_Region_Collate
  END INTERFACE

  INTERFACE PERMUTE
     MODULE PROCEDURE RegionPermute
  END INTERFACE

  INTERFACE BC_Region_UseFunction
     MODULE PROCEDURE UseFunction
  END INTERFACE

  INTERFACE BC_Get_Name
     MODULE PROCEDURE RegionGetName
  END INTERFACE
  INTERFACE BC_Set_Name
     MODULE PROCEDURE RegionsetName
  END INTERFACE

  INTERFACE BC_Get_DOF
     MODULE PROCEDURE RegionGetDegreesOfFreedom
  END INTERFACE
  INTERFACE BC_Set_DOF
     MODULE PROCEDURE RegionSetDegreesOfFreedom
  END INTERFACE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  !====================================================================
  !
  ! Routines For Maintaining and Accessing BC_Region
  !
  !====================================================================

  function BC_Region_Size(Region)
    type(BC_Region), intent(IN) :: Region
    integer :: BC_Region_Size
    BC_Region_Size = Region%NumberOfEntries
  end function BC_Region_Size
  
  function BC_Region_Dimension(Region)
    type(BC_Region), intent(IN) :: Region
    integer :: BC_Region_Dimension
    BC_Region_Dimension = SIZE(BC_Region_Positions(Region),1)
  end function BC_Region_Dimension
  
  function BC_Region_OwnerCells(Region)
    type(BC_Region), intent(IN) :: Region
    integer, POINTER, dimension(:) :: BC_Region_OwnerCells
    BC_Region_OwnerCells => Region%OwnerCell
  end function BC_Region_OwnerCells

  function BC_Region_Faces(Region)
    type(BC_Region), intent(IN) :: Region
    integer, POINTER, dimension(:) :: BC_Region_Faces
    BC_Region_Faces => Region%Face
  end function BC_Region_Faces

  function BC_Region_Values(Region)
    type(BC_Region), intent(IN) :: Region
    real(r8), POINTER, dimension(:,:) :: BC_Region_Values
    BC_Region_Values => Region%Value
  end function BC_Region_Values

  function UseFunction(Region) 
    type(BC_Region), intent(IN) :: Region
    logical, POINTER, dimension(:) :: UseFunction
    UseFunction => Region%UseFunction
  end function UseFunction

  function BC_Region_Positions(Region)
    type(BC_Region), intent(IN)  :: Region
    real(r8), POINTER, dimension(:,:) :: BC_Region_Positions
    BC_Region_Positions => Region%Position
  end function BC_Region_Positions

  function RegionGetName(Region) RESULT(NAME)
    type(BC_Region), intent(IN) :: Region
    character(BC_STRING_LEN) :: NAME
    NAME = TRIM(Region%Name)
  END function RegionGetName

  subroutine RegionSetName(Region, NAME)
    type(BC_Region), intent(INOUT)  :: Region
    character(*), intent(IN) :: NAME
    Region%Name = TRIM(NAME)
  END subroutine RegionSetName

  function RegionGetDegreesOfFreedom(Region) RESULT(DOF)
    type(BC_Region), intent(IN) :: Region
    integer :: DOF
    DOF = Region%DegreesOfFreedom
  END function RegionGetDegreesOfFreedom

  subroutine RegionSetDegreesOfFreedom(Region, DOF)
    type(BC_Region), intent(INOUT) :: Region
    integer, intent(IN) :: DOF
    Region%DegreesOfFreedom = DOF
  END subroutine RegionSetDegreesOfFreedom

  subroutine BC_Region_Initialize(Region)
    ! Set the size of the region to -1.
    ! This does not deallocate any memory, so it is appropraite to call
    ! this routine on new structure.  It is an error to call this routine on
    ! a Region which has allocated memory.  Call FreeBC_Region instead.
    type(BC_Region), intent(INOUT) :: Region

    Region%NumberOfEntries = BC_INVALID_COUNT
    Region%DegreesOfFreedom = BC_INVALID_COUNT

    NULLIFY(Region%OwnerCell)
    NULLIFY(Region%Face)
    NULLIFY(Region%Value)
    NULLIFY(Region%UseFunction)
    NULLIFY(Region%Position)
    Call BC_Set_Name(Region, BC_NO_NAME)

  end subroutine BC_Region_Initialize
  
  subroutine BC_Region_Free(Region)
    ! Free the memory used by a BC_Region, and initialize the
    ! region to be empty.
    type(BC_Region), intent(INOUT) :: Region

    if (ASSOCIATED(Region%OwnerCell)) DEALLOCATE(Region%OwnerCell)
    if (ASSOCIATED(Region%Face))      DEALLOCATE(Region%Face)
    if (ASSOCIATED(Region%Value))     DEALLOCATE(Region%Value)
    if (ASSOCIATED(Region%UseFunction)) DEALLOCATE(Region%UseFunction)
    if (ASSOCIATED(Region%Position))  DEALLOCATE(Region%Position)

    call BC_Region_Initialize(Region)

  end subroutine BC_Region_Free

  subroutine BC_Region_Alloc(Region, SIZE, DIMENSIONALITY, DOF)
    ! ALLOCATE memory internal to a BC_Region
    type(BC_Region), intent(INOUT) :: Region
    integer, intent(IN) :: SIZE
    integer, intent(IN) :: DIMENSIONALITY
    integer, intent(IN) :: DOF

    Region%NumberOfEntries = SIZE
    Region%DegreesOfFreedom = DOF

    ALLOCATE(Region%OwnerCell(SIZE))
    Region%OwnerCell = BC_INVALID_CELL
    ALLOCATE(Region%Face(SIZE))
    ALLOCATE(Region%Value(DOF, SIZE))
    ALLOCATE(Region%UseFunction(SIZE))
    ALLOCATE(Region%Position(DIMENSIONALITY, SIZE))

  end subroutine BC_Region_Alloc

  subroutine BC_Region_Copy(Dest_Region, Src_Region)
    ! Copy Dest_Region into Src_Region
    ! The Dest_Region must have data fields previously allocated.
    ! The Dest_Region must be at least as big as the Src_Region.
    type(BC_Region), intent(INOUT) :: Dest_Region
    type(BC_Region), intent(IN   ) :: Src_Region

    integer :: DestSize, SrcSize
    DestSize = SIZE(Dest_Region)
    SrcSize  = SIZE(Src_Region)
#ifdef BC_TESTS_ON
    if (DestSize < SrcSize) then
       print *, "ERROR: BC_Region_Copy, Dest region is smaller than Src region."
    end if
    if (DIMENSIONALITY(Dest_Region) /= DIMENSIONALITY(Src_Region)) then
       print *, "ERROR: BC_Region_Copy, Dest region and Src Region do not have same DIMENSIONALITY."
    end if
    if (BC_Get_DOF(Dest_Region) /= BC_Get_DOF(Src_Region)) then
       print *, "ERROR: BC_Region_Copy, Dest region and Src Region do not have same Degrees Of Freedom."
    end if
#endif

    Dest_Region%OwnerCell(1:SrcSize)   = BC_Region_OwnerCells(Src_Region)
    Dest_Region%Face(1:SrcSize)        = BC_Region_Faces(Src_Region)
    Dest_Region%Value(:,1:SrcSize)     = BC_Region_Values(Src_Region)
    Dest_Region%UseFunction(1:SrcSize) = BC_Region_UseFunction(Src_Region)
    Dest_Region%Position(:,1:SrcSize)  = BC_Region_Positions(Src_Region)

    Call BC_Set_Name(Dest_Region, BC_Get_Name(Src_Region))
  end subroutine BC_Region_Copy
  
  subroutine BC_Region_Insert(Region, CellList, FaceList, ValueList, UseFunctionList, PositionList)
    ! Insert data lists onto a BC_Region
    ! If Region has slots for REGION_SIZE entries, the lists must have at least REGION_SIZE
    ! items.  This routine inserts the first REGION_SIZE items from each list into Region

    ! Parameter List
    type (BC_Region), intent(INOUT) :: Region
    integer,  intent(IN) :: CellList(:)
    integer,  intent(IN) :: FaceList(:)
    real(r8), intent(IN) :: ValueList(:,:)
    logical,  intent(IN) :: UseFunctionList(:)
    real(r8), intent(IN) :: PositionList(:,:)

    ! Local variables
    integer :: region_size, insert_size

    region_size = BC_Region_Size(Region)
    insert_size = SIZE(CellList)

#ifdef BC_TESTS_ON
    ! Test that the lists all consistent
    if (.NOT. ListAllSameSize(CellList, FaceList, ValueList, UseFunctionList, PositionList) ) then
       print *, "ERROR: BC_Region_Insert lists not all same size"
    endif

    ! Test that the Positions have the proper dimensions
    if (.NOT. ConsistentDimensions(Region, PositionList)) then
       print *, "ERROR: BC_Region_Insert PositionList different DIMENSIONALITY from Region"
    endif

    ! Test that the list sizes are large enough
    if (region_size > insert_size) then
       print *, "ERROR: BC_Region_Insert Lists not large enough"
    end if
#endif

    ! Put the first insert_size items from the list into Region
    Region%OwnerCell(1:insert_size)  = CellList(1:insert_size)
    Region%Face(1:insert_size)       = FaceList(1:insert_size)
    Region%Value(:,1:insert_size)      = ValueList(:,1:insert_size)
    Region%UseFunction(1:insert_size)= UseFunctionList(1:insert_size)
    Region%Position(:,1:insert_size) = PositionList(:,1:insert_size)

  end subroutine BC_Region_Insert

  subroutine BC_Region_REAlloc(Region, SIZE)
    ! RE-Allocate a BC_Region.  The returned region has the data
    ! of the input region as the first items in its lists.  The
    ! tails of the lists are empty.
    ! Parameter list
    type(BC_Region), intent(INOUT) :: Region
    integer, intent(IN) :: SIZE

    ! Local variables
    integer :: old_size, region_dim, region_DOF
    type(BC_Region) :: Temp_Region

    old_size   = BC_Region_Size(Region)
    region_dim = BC_Region_Dimension(Region)
    region_DOF = BC_Get_DOF(Region)
    
    ! Make a temporary region to hold the current region
    call INITIALIZE(Temp_Region)
    call ALLOC(Temp_Region, SIZE=old_size, DIMENSIONALITY = region_dim, DOF = region_DOF)
    Temp_Region = Region

    ! Now enlarge Region to have room for the new lists
    call FREE(Region)
    call ALLOC(Region, SIZE=SIZE, DIMENSIONALITY = region_dim, DOF = region_DOF)

    ! And Copy temp into new Region
    Region = Temp_Region

    ! Done with Temp_Region, so release it
    call FREE(Temp_Region)

  end subroutine BC_Region_REAlloc

  subroutine BC_Region_Append(Region, CellList, FaceList, ValueList, UseFunctionList, PositionList)
    ! Append data lists onto a BC_Region.  The returned list is re-sized
    ! to hold all the original data and the new data.  The returned list
    ! has the original data at the head of each internal list, and the new
    ! list data at the tail of each list.
    ! Parameter List
    type (BC_Region), intent(INOUT) :: Region
    integer,  intent(IN) :: CellList(:)
    integer,  intent(IN) :: FaceList(:)
    real(r8), intent(IN) :: ValueList(:,:)
    logical,  intent(IN) :: UseFunctionList(:)
    real(r8), intent(IN) :: PositionList(:,:)

    ! Local variables
    integer :: old_size, increment, new_size, region_dim
    integer :: start, end

    ! Start of subroutine
    old_size    = BC_Region_Size(Region)
    region_dim  = BC_Region_Dimension(Region)
    increment   = SIZE(CellList)
    
    ! Test for consitency
#ifdef ADD_DEBUG_CODE_HERE
    ! Test that the lists all consistent
    if (.NOT. ListAllSameSize(CellList, FaceList, ValueList, UseFunctionList, PositionList) ) then
       print *, "ERROR: BC_Region_Append lists not all same size"
    endif

    ! Test that the Positions have the proper dimensions
    if (.NOT. ConsistentDimensions(Region, PositionList)) then
       print *, "ERROR: BC_Region_Append PositionList different DIMENSIONALITY from Region"
    endif
#endif

    ! The size of the new region is large enough for the existing region and the
    ! lists to be appended
    new_size = old_size + increment
    
    ! Reallocate Region to be large enough to hold the new lists
    ! This returns a new region sized to hold the old and new lists,
    ! with the original data at the head of the lists
    call REALLOC(Region, SIZE=new_size)

    ! Now append the new data to Region
    start = old_size + 1
    end   = new_size
    Region%OwnerCell(start:end)  = CellList
    Region%Face(start:end)       = FaceList
    Region%Value(:,start:end)      = ValueList
    Region%UseFunction(start:end)= UseFunctionList
    Region%Position(:,start:end) = PositionList

  end subroutine BC_Region_Append

  function ListAllSameSize(CellList, FaceList, ValueList, UseFunctionList, PositionList)
    ! Test that all the given lists have the same size.
    ! This returns TRUE if that is the case, otherwise FALSE.
    ! Parameter List
    logical :: ListAllSameSize
    integer,  intent(IN) :: CellList(:)
    integer,  intent(IN) :: FaceList(:)
    real(r8), intent(IN) :: ValueList(:,:)
    logical,  intent(IN) :: UseFunctionList(:)
    real(r8), intent(IN) :: PositionList(:,:)

    ! Local variables
    integer :: CellList_Size, FaceList_Size
    integer :: ValueList_Size, UseFunctionList_Size, PositionList_Size

    ! Start of subroutine
    CellList_Size        = SIZE(CellList)
    FaceList_Size        = SIZE(FaceList)
    ValueList_Size       = SIZE(ValueList,2)
    UseFunctionList_Size = SIZE(UseFunctionList)
    PositionList_Size    = SIZE(PositionList,2)
    
    ! Note that this test by associativity this test is both necessary and sufficient
    ListAllSameSize =       (CellList_Size == FaceList_Size)        &
                      .AND. (CellList_Size == ValueList_Size)       &
                      .AND. (CellList_Size == UseFunctionList_Size) &
                      .AND. (CellList_Size == PositionList_Size) 

  end function ListAllSameSize

  function ConsistentDimensions(Region, PositionList)
    ! Returns true if the Position field of Region and the PositionList
    ! have the same first dimension
    ! Parameter List
    logical :: ConsistentDimensions
    type(BC_Region), intent(INOUT) :: Region
    real(r8), dimension(:,:), intent(IN) :: PositionList
    ConsistentDimensions = BC_Region_Dimension(Region) == SIZE(PositionList,1)
  end function ConsistentDimensions

  subroutine BC_Region_Canonical(Region)
    !-----------------------------------------------------------------------------
    ! Purpose:
    !   Re-order the region so that it is canonical form.
    !   At the moment that means that the entries are in (cell,face) order.
    !-----------------------------------------------------------------------------
    ! Arguments
    type(BC_Region), intent(INOUT) :: Region
    call ORDER(Region)
  end subroutine BC_Region_Canonical
  

  subroutine BC_Region_Order(Region)
    !-----------------------------------------------------------------------------
    ! Purpose:
    !   Order the region so that it is in increasing cell number,
    !   and within a cell, in increasing face number.
    !
    !-----------------------------------------------------------------------------
    use parameter_module
    use pgslib_module, only : PGSLib_GRADE_UP, PGSLib_Global_EOSHIFT,  &
                       PGSLIB_PARITY_PREFIX

    ! Arguments
    type(BC_Region), intent(INOUT) :: Region

    ! Local variables
    integer, dimension(:), pointer :: Cells
    integer, dimension(:), pointer :: Rank
    logical, dimension(:), pointer :: Mask, Segment

    ! Before we stuff this into an Atlas, want to get data into (Face,Cell) order.

    ALLOCATE(Rank(SIZE(Region)))
    ALLOCATE(Mask(SIZE(Region)))
    ALLOCATE(Segment(SIZE(Region)))

    ! First sort by cell number, then segment sort by face number
    Cells => BC_Region_OwnerCells(Region)
    Rank = PGSLib_GRADE_UP(Cells)
    call PERMUTE(Region, Rank)

    ! Construct segments at cell boundaries.
    Mask = Cells /= PGSlib_GLOBAL_EOSHIFT(Cells, SHIFT=-1, BOUNDARY = BC_INVALID_CELL)
    Segment = PGSLIB_PARITY_PREFIX(Segment)

    ! Now sort the segments
    Rank = PGSLib_GRADE_UP(Cells, Segment)
    call PERMUTE(Region, Rank)

    ! Done with everything
    DEALLOCATE(Segment, Mask, Rank)

  end subroutine BC_Region_Order

  subroutine BC_Region_Collate(Collated_Region, Local_Region)
    ! Collate Local_Region into Collated_Region.
    ! This routine assumes that Collated_Region has been appropriately allocated,
    ! so Collated_Region is an INOUT argument.

    use parallel_util_module
    use pgslib_module,       ONLY: PGSLib_Global_SUM, PGSlib_Collate
    type(BC_Region), intent(INOUT) :: Collated_Region
    type(BC_Region), intent(IN   ) :: Local_Region

    ! Local variables
    integer :: Collated_Size, d, DOF
    integer,  pointer, dimension(:)   :: CollatedCell, LocalCell
    integer,  pointer, dimension(:)   :: CollatedFace, LocalFace
    real(r8), pointer, dimension(:,:) :: CollatedValue, LocalValue
    logical,  pointer, dimension(:)   :: CollatedUseF, LocalUseF
    real(r8), pointer, dimension(:,:) :: CollatedPosition, LocalPosition

    Collated_Size = PGSLib_Global_SUM(SIZE(Local_Region))
    if (.NOT. Is_IO_PE()) Collated_Size = 0
    DOF = BC_Get_DOF(Local_Region)

    ALLOCATE(CollatedCell(Collated_Size))
    ALLOCATE(CollatedFace(Collated_Size))
    ALLOCATE(CollatedValue(DOF,Collated_Size))
    ALLOCATE(CollatedUseF(Collated_Size))
    ALLOCATE(CollatedPosition(DIMENSIONALITY(Local_Region), Collated_Size))

    ! Collate the lists
    LocalCell => BC_Region_OwnerCells(Local_Region)
    call PGSLib_Collate(CollatedCell, LocalCell)
    LocalFace => BC_Region_Faces(Local_Region)
    call PGSLib_Collate(CollatedFace, LocalFace)
    LocalValue => BC_Region_Values(Local_Region)
    do d = 1, DOF
       call PGSLib_Collate(CollatedValue(d,:), LocalValue(d,:))
    end do
    LocalUseF => BC_Region_UseFunction(Local_Region)
    call PGSLib_Collate(CollatedUseF, LocalUseF)
    LocalPosition => BC_Region_Positions(Local_Region)
    do d = 1, DIMENSIONALITY(Local_Region)
       call PGSLib_Collate(CollatedPosition(d,:), LocalPosition(d,:))
    end do

    ! Now insert the collated lists into the collated region
    if (Is_IO_PE()) then
       call INSERT(Collated_Region, CellList        = CollatedCell, &
                                    FaceList        = CollatedFace, &
                                    ValueList       = CollatedValue, &
                                    UseFunctionList = CollatedUseF, &
                                    PositionList    = CollatedPosition)
    end if

    call BC_Set_Name(Collated_Region, BC_Get_Name(Local_Region))
    ! Clean up 
    DEALLOCATE(CollatedPosition)
    DEALLOCATE(CollatedValue)
    DEALLOCATE(CollatedUseF)
    DEALLOCATE(CollatedFace)
    DEALLOCATE(CollatedCell)

  END subroutine BC_Region_Collate    

  subroutine RegionPermute(Region, Rank)
    ! Permute a region by permute vector rank.  By default permuting means
    ! to re-order the data.
    use pgslib_module, ONLY: PGSLib_GS_Trace, PGSLib_Permute, PGSLib_DEALLOCATE_TRACE
    type(BC_Region), intent(INOUT) :: Region
    integer, dimension(:), intent(IN) :: Rank

    ! Local variables
    integer :: d
    integer,  pointer, dimension(:)   :: Cell, Face
    real(r8), pointer, dimension(:,:) :: Value
    logical,  pointer, dimension(:)   :: UseF
    real(r8), pointer, dimension(:,:) :: Position

    integer,  dimension(SIZE(Region%Face)) :: TempCell, TempFace
    real(r8), dimension(SIZE(Region%Value,1),SIZE(Region%Value,2)) :: TempValue
    logical,  dimension(SIZE(Region%Value,2)) :: TempLog
    real(r8), dimension(SIZE(Region%Position,1), SIZE(Region%Position,2)) :: TempPosition
    
    type (PGSLib_GS_Trace), POINTER :: Permute_Trace

    Cell     => BC_Region_OwnerCells(Region)
    Face     => BC_Region_Faces(Region)
    Value    => BC_Region_Values(Region)
    UseF     => BC_Region_UseFunction(Region)
    Position => BC_Region_Positions(Region)


    ! Permute data into temporary space, then overwrite
    NULLIFY(Permute_Trace)
    call PGSLib_PERMUTE(DEST   = TempCell, &
                        SOURCE = Cell,     &
                        INDEX  = Rank,     &
                        TRACE  = Permute_Trace)
    Cell = TempCell

    call PGSLib_PERMUTE(DEST   = TempFace, &
                        SOURCE = Face,     &
                        INDEX  = Rank,     &
                        TRACE  = Permute_Trace)
    Face = TempFace

    do d = 1, BC_Get_DOF(Region)
       call PGSLib_PERMUTE(DEST   = TempValue(d,:), &
                           SOURCE = Value(d,:),     &
                           INDEX  = Rank,           &
                           TRACE  = Permute_Trace)
       Value(d,:) = TempValue(d,:)
    end do

    call PGSLib_PERMUTE(DEST   = TempLog,   &
                        SOURCE = UseF,      &
                        INDEX  = Rank,      &
                        TRACE  = Permute_Trace)
    UseF = TempLog

    do d = 1, DIMENSIONALITY(Region)
       call PGSLib_PERMUTE(DEST   = TempPosition(d,:), &
                           SOURCE = Position(d,:),     &
                           INDEX  = Rank,              &
                           TRACE  = Permute_Trace)
       Position(d,:) = TempPosition(d,:)
    end do

    !Done with trace
    call PGSLib_DEALLOCATE_TRACE(Permute_Trace)
    
  end subroutine RegionPermute

END Module BC_Regions

