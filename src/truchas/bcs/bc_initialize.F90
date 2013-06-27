Module BC_Initialize

  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide some stuff for initializing BC's.  This is for
  !   "New BCs".  These are support routines.  Other modules
  !   have the variable specific stuff.
  !   These routines require telluride/trucahs specific modules.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  use bc_data_types
  use parameter_module
  implicit none
  Private

  PUBLIC :: APPEND
  PUBLIC :: BC_Atlas_From_Region

  INTERFACE APPEND
     MODULE PROCEDURE Append_BC_Region_From_Mask
  END INTERFACE
  

CONTAINS
  subroutine Append_BC_Region_From_Mask(Region, BC_Mask, BC_Values, BC_UseFunction, POSITIONS)
    ! Append the boundary values flagged in BC_Mask to Region
    use pgslib_module, ONLY: PGSLib_SUM_PREFIX
    use mesh_module

    type(BC_Region), intent(INOUT) :: Region
    logical, intent(IN), dimension(nfc,ncells) :: BC_Mask
    ! BC_Values is (DOF, nfc, ncells)
    real(r8), intent(IN), dimension(:,:,:) :: BC_Values
    logical,  intent(IN), dimension(nfc,ncells) :: BC_UseFunction
    real(r8), intent(IN), dimension(ndim, nfc,ncells) :: Positions

    integer :: ItemsInList, C, F, D, ListCounter, DOF
    integer, dimension(ncells) :: Global_Cell_Number
    integer :: Global_Offset
    integer, pointer, dimension(:) :: CellList
    integer, pointer, dimension(:) :: FaceList
    real(r8), pointer, dimension(:,:) :: ValueList
    logical,  pointer, dimension(:) :: UseFunction
    real(r8), pointer, dimension(:,:) :: PositionList

    ! How many BC points do we get from this side?
    ItemsInList = COUNT(BC_Mask)

    ! How many degrees of freedom?
    DOF = BC_Get_DOF(Region)

    ! Make space for the lists.
    ALLOCATE(CellList(ItemsInList), FaceList(ItemsInList), &
             ValueList(DOF, ItemsInList), UseFunction(ItemsInList), PositionList(ndim, ItemsInList))

    ! Need the global cell number to put into region.
    ! That is just local cell number + global offset
    Global_Cell_Number = PGSLib_SUM_PREFIX( (/ (1, c = 1, ncells) /) )
    Global_Offset      = Global_Cell_Number(1) - 1

    ! Pack the data from the Mesh into the lists
    ListCounter = 1
    do C = 1, ncells
       do F = 1, nfc
          if (BC_Mask(F, C)) then
             CellList(ListCounter) = C + Global_Offset
             FaceList(ListCounter) = F
             ValueList(:,ListCounter) = BC_Values(:,F,C)
             UseFunction(ListCounter) = BC_UseFunction(F,C)
             ! Position set to face centroid, but IS THIS CORRECT???
             do d = 1, ndim
                PositionList(d, ListCounter) = Positions(d, F, C)
             end do
             ListCounter = ListCounter + 1
          end if
       end do
    end do

    ! Append the lists onto the region
    call APPEND(REGION       = Region,     &
                CELLLIST     = CellList,   &
                FACELIST     = FaceList,   &
                VALUELIST    = ValueList,  &
                UseFunctionList = UseFunction, &
                POSITIONLIST = Positionlist)

    ! Coda
    ! Done with the lists, clean up and go home.
    DEALLOCATE(CellList, FaceList, ValueList, UseFunction, PositionList)

  end subroutine Append_BC_Region_From_Mask

  subroutine BC_Atlas_From_Region(Atlas, Region)
    !-----------------------------------------------------------------------------
    ! Purpose:
    !   Initialize whatever is needed for the BC Atlas.
    !   given a Region.
    !   Eventually this routine will be vanish and a general
    !   BC initialization will be done.
    !
    !-----------------------------------------------------------------------------
    ! Arguments
    type (BC_Atlas), intent(INOUT)  :: Atlas
    type (BC_Region), intent(IN   ) :: Region

    ! Local variables
    integer :: c, d
    
    integer,  dimension(:), pointer :: Faces, Cells
    real(r8), dimension(:,:), pointer :: Values
    logical,  pointer, dimension(:) :: UseFunction
    real(r8), dimension(:,:), pointer :: Positions
    integer,  dimension(1) :: ChartValueIndex
    real(r8), dimension(ndim,1) :: ChartPositions

    ! Build the Atlas from the Region data

    call INITIALIZE(Atlas)
    
    call BC_Set_Name(Atlas, BC_Get_Name(Region))
    ! Pick a size that should be large enough, okay if it isn't.
    ! In this simple case each point in the region generates a single
    ! chart of size 1.
    call ALLOC(Atlas, MaxNumberOfCharts   = SIZE(Region), &
                      MaxAtlasDataSize    = SIZE(Region), &
                      DIMENSIONALITY = DIMENSIONALITY(Region),&
                      DOF = BC_Get_DOF(Region))

    ! Cell numbers are global

    call SET_SCOPE(Atlas, SCOPE = BC_Cells_Scope_Global)

    ! To set the atlas we will loop over all the points in the region

    Faces       => BC_Region_Faces(Region)
    Cells       => BC_Region_OwnerCells(Region)
    Values      => BC_Region_Values(Region)
    UseFunction => BC_Region_UseFunction(Region)
    Positions   => BC_Region_Positions(Region)

    do c = 1, SIZE(Faces)
       ! all charts have size 1, so this is easy
!       ChartValue = Values(c)
       ! Value Index is not functioning yet.  The problem is we might want that to be local,
       ! but don't know how to assure both it at Cells are local if they aren't the same.
       ! For now assume they are the same, and don't use ValueIndex for anything.
       ! Put BC_INVALID_CELL in it to catch misuses.
       ! ChartValueIndex = Cells(c)
       ChartValueIndex = BC_INVALID_CELL
       do d = 1, ndim
          ChartPositions(d,1) = Positions(d,c)
       end do
       call APPEND(Atlas,                            &
                   VALUES      = Values(:,c:c),      &
                   USEFUNCTION = (/UseFunction(c)/), &
                   VALUEINDEX  = ChartValueIndex,   &
                   POSITIONS   =    ChartPositions, &
                   CELL        =    Cells(c),       &
                   FACE        =    Faces(c))
    end do
    
  end subroutine BC_Atlas_From_Region

END Module BC_Initialize


    
    
    
