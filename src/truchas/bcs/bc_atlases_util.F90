!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_ATLASES_Util
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the routines to manipulate Atlases.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  use bc_atlases_data_types
  use bc_enum_types
  use truchas_logging_services
  Implicit None
  Private

  PUBLIC :: CANONICALIZE
  PUBLIC :: COLLATE
  PUBLIC :: ORDER
  PUBLIC :: PERMUTE
  PUBLIC :: RENUMBER_CELLS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE CANONICALIZE
     MODULE PROCEDURE AtlasCanonical
  END INTERFACE
  
  INTERFACE ORDER
     MODULE PROCEDURE AtlasOrder
  END INTERFACE

  INTERFACE LOCALIZE
     MODULE PROCEDURE AtlasLocalize
  END INTERFACE

  INTERFACE PERMUTE
     MODULE PROCEDURE AtlasPermute
  END INTERFACE

  INTERFACE REDISTRIBUTE
     MODULE PROCEDURE AtlasRedistribute
  END INTERFACE

  INTERFACE COLLATE
     MODULE PROCEDURE AtlasCollate
  END INTERFACE

  INTERFACE RENUMBER_CELLS
     MODULE PROCEDURE AtlasRenumberCells
  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routines to canonicalize an atlas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AtlasCanonical(Atlas)
    ! Permute, sort, re-arrange Atlas so it is in canonical order.
    type(BC_Atlas), intent(INOUT), target :: Atlas
    
    ! An atlas is canonical if 
    ! 1) All data for a chart is on the same processor as 
    !    the chart reference point.
    ! 2) Locally on each processor the charts are ordered by
    !    increasing cell number, and within a cell by face.
    ! 3) Cell numbers refer to local cell numbers, not global cell
    !    numbers.

    ! First we need to permute the Atlas so that the data
    ! is local to the cell it applies to.

    call LOCALIZE(Atlas)
    call ORDER(Atlas)
    call RENUMBER_CELLS(Atlas, SCOPE = BC_Cells_Scope_Local)

  end subroutine AtlasCanonical

  subroutine AtlasLocalize(Atlas)
    ! Permute the Atlas (just the AtlasData, actually) so that
    ! the data is local to the processor it applies to.  More specifically, 
    ! permute the atlas so that all the data for each chart is local
    ! to the processor which owns that chart.  Chart ownership
    ! is determined by the cell identifier for that chart.
    use legacy_mesh_api, only: ncells
    use parallel_communication, only: this_PE
    use pgslib_module, ONLY: PGSLib_Gather, PGSLib_GRADE_UP, PGSLib_Scatter_SUM                                  
    type(BC_Atlas), intent(INOUT), target :: Atlas

    ! Local variables
    integer :: i, NewDataSize
    integer, allocatable, dimension(:) :: ChartProcNum
    integer, allocatable, dimension(:) :: Rank
    integer, dimension(1) :: ItemsThisProc
!!$    integer, allocatable, dimension(:) :: ItemsPerProc
!!$    logical, allocatable, dimension(:) :: ProcNumSeg
!!$    logical, allocatable, dimension(:) :: ProcNumMask
    integer, dimension(ncells)         :: CellProcNum
    integer, POINTER, dimension(:)     :: Cells
    type(BC_Atlas)                    :: Temp_Atlas
    integer                            :: status

    allocate (ChartProcNum(SIZE(Atlas)), STAT=status)
    if (status /= 0) call TLS_panic ('AtlasLocalize: allocate failed: ChartProcNum')
    allocate (Rank(SIZE(Atlas)), STAT=status)
    if (status /= 0) call TLS_panic ('AtlasLocalize: allocate failed: Rank')
!!$    allocate (ItemsPerProc(SIZE(Atlas)), STAT=status)
!!$    if (status /= 0) call PUNT((/'allocate failed: ItemsPerProc'/), 'AtlasLocalize')
!!$    allocate (ProcNumMask(SIZE(Atlas)), STAT=status)
!!$    if (status /= 0) call PUNT((/'allocate failed: ProcNumMask'/), 'AtlasLocalize')
!!$    allocate (ProcNumSeg(SIZE(Atlas)), STAT=status)
!!$    if (status /= 0) call PUNT((/'allocate failed: ProcNumSeg'/), 'AtlasLocalize')

    ! First we have to find the processor number for each chart
    ! At this stage Cells refers to global cell number, since Atlas is
    ! not yet canonical.
    Cells => BC_Get_Cell(Atlas)
    ! Find the processor number for each cell.
    CellProcNum = this_PE
    ! Gather that processor number to the charts
    call PGSLib_GATHER(DEST   = ChartProcNum, &
                       SOURCE = CellProcNum,  &
                       INDEX  = Cells)

    ! Now sort by processor number, so that all charts with same processor number
    ! get together.
    Rank = PGSLib_Grade_UP(ChartProcNum)

    ! And permute atlas to get it into this order.  No point in re-indexing 
    ! at this moment, since lots of finagling left to do.
    call PERMUTE(Atlas, Rank, REINDEX = .FALSE.)
    
    ! Finally, we need to redistribute the atlas data so that it
    ! is on the appropriate processor.

!!$    ! We need the ChartProcNum, but it isn't in the appropriate order.
!!$    ! We could permute it, or re-get it.  Choose latter.
!!$    call PGSLib_GATHER(DEST   = ChartProcNum, &
!!$                       SOURCE = CellProcNum,  &
!!$                       INDEX  = Cells)
!!$    
!!$    ! Now need to count the number of data items on each processor.
!!$    ProcNumMask = ChartProcNum /= PGSLib_Global_EOSHIFT(ChartProcNum, SHIFT = -1, BOUNDARY = -1)
!!$    ProcNumSeg  = PGSLib_PARITY_PREFIX(ProcNumMask)
!!$    ItemsPerProc = PGSLib_SUM_SUFFIX((/ (1, i=1,SIZE(ItemsPerProc)) /), SEGMENT = ProcNumSeg)
!!$    call PGSLib_GLOBAL_PACK(DEST = ItemsThisProc, SOURCE = ItemsPerProc, MASK=ProcNumMask)
    
    
    ! Now need to count the number of data items on each processor.
    ! Do that by scattering 1 for each item to the appropriate processor.
    ItemsThisProc = 0
    call pgslib_scatter_sum(DEST = ItemsThisProc, &
                            SOURCE = (/ (1, i=1, SIZE(ChartProcNum)) /), &
                            INDEX = ChartProcNum)

    NewDataSize = ItemsThisProc(1)

    ! Re-Distribute moves data from one atlas to another.  To use it we need
    ! an atlas of the proper distribution.  However, we want to return the input
    ! atlas.  So, first copy the input atlas to a temporary, then delete and allocate it
    ! newly with the new size.
    call CLONE(Temp_Atlas, Atlas)

    call FREE(Atlas)
    ! We don't care how many charts we allow for, since we will re-index the charts
    ! soon.
    call ALLOC(Atlas, MaxNumberOfCharts   = 0,            &
                      MaxAtlasDataSize    = NewDataSize,  &
                      DIMENSIONALITY = DIMENSIONALITY(Temp_Atlas),&
                      DOF = BC_Get_DOF(Temp_Atlas))
    call BC_Set_Name(Atlas, BC_Get_Name(Temp_Atlas))

    ! Finally, redistribute Temp_Atlas into Atlas
    call REDISTRIBUTE(DEST = Atlas, SOURCE = Temp_Atlas)

    ! And get rid of Temp_Atlas
    call FREE(Temp_Atlas)

!!$    deallocate (ProcNumSeg)
!!$    deallocate (ProcNumMask)
!!$    deallocate (ItemsPerProc)
    deallocate (Rank)
    deallocate (ChartProcNum)

  end subroutine AtlasLocalize

    

  subroutine AtlasOrder(Atlas)
    ! Sort an atlas so that it is in canonical order.  That means
    ! ordered by cells, and then for each cell by face number.
    use PGSLib_module, ONLY: PGSLib_GRADE_UP, PGSLib_Global_EOSHIFT, PGSLib_PARITY_PREFIX
    type(BC_Atlas), intent(INOUT), target :: Atlas

    ! Local variables
    integer, allocatable, dimension(:) :: Rank
    integer, allocatable, dimension(:) :: CellSegNum
    logical, allocatable, dimension(:) :: CellSeg
    logical, allocatable, dimension(:) :: CellMask
    integer, POINTER,     dimension(:) :: Cells, Faces
    integer :: status

    allocate (rank(SIZE(Atlas)), STAT=status)
    if (status /= 0) call TLS_panic ('AtlasOrder: allocate failed: Rank')
    allocate (CellSegNum(SIZE(Atlas)), STAT=status)
    if (status /= 0) call TLS_panic ('AtlasOrder: allocate failed: CellSegNum')
    allocate (CellSeg(SIZE(Atlas)), STAT=status)
    if (status /= 0) call TLS_panic ('AtlasOrder: allocate failed: CellSeg')
    allocate (CellMask(SIZE(Atlas)), STAT=status)
    if (status /= 0) call TLS_panic ('AtlasOrder: allocate failed: CellMask')

    ! First we rank each data item by cell number
    ! Then we permute the atlas to put it into order.
    Cells => BC_Get_Cell(Atlas)

    Rank  = PGSLib_GRADE_UP(Cells)

    call PERMUTE(Atlas, Rank, REINDEX = .FALSE.)

    ! Then we segment rank each data item in a cell by face number.
    ! Then we permute the atlas to put it into order.
    Faces => BC_Get_Face(Atlas)

    CellMask = Cells /= PGSLib_Global_EOSHIFT(Cells, SHIFT = -1, BOUNDARY = -1)

    CellSeg  = PGSLib_PARITY_PREFIX(CellMask)

    Rank     = PGSLib_GRADE_UP(Faces, SEGMENT = CellSeg)
    ! Since this is the last permutation, recompute the Offset and Length fields
    call PERMUTE(Atlas, Rank, REINDEX = .TRUE.)

    deallocate (CellMask)
    deallocate (CellSeg)
    deallocate (CellSegNum)
    deallocate (rank)

  end subroutine AtlasOrder

  subroutine AtlasPermute(Atlas, Rank, REINDEX)
    ! Permute an atlas by permute vector rank.  By default permuting means
    ! to re-order the data, and also recompute the Length and Offset fields.
    ! If REINDEX is present and false, then the Length and Offset fields
    ! are not defined and not correct after this call.
    type(BC_Atlas), intent(INOUT), target :: Atlas
    integer, dimension(:), intent(IN) :: Rank
    logical, OPTIONAL,     intent(IN) :: REINDEX

    ! Local variables
    type(BC_Atlas_Data), POINTER :: AtlasData
    logical :: ReCompute

    AtlasData => BC_Get_Data(Atlas)
    call PERMUTE(AtlasData, Rank)

    ! Do we recompute indices?
    if (PRESENT(REINDEX)) then
       ReCompute = REINDEX
    else
       ReCompute = .TRUE.
    end if
    if (ReCompute) then
       call ComputeIndex(Atlas)
    end if

  end subroutine AtlasPermute

  subroutine AtlasRedistribute(Dest, Source)
    ! Redistribute the source atlas into the dest atlas.  The redistribution
    ! does not do anything with the AtlasSpec component.  This routine
    ! must be followed by a ComputeIndex call to finish constructing the atlas.
    type(BC_Atlas), intent(INOUT), target :: Dest
    type(BC_Atlas), intent(INOUT), target :: Source

    ! Local variables
    type(BC_Atlas_Data), POINTER :: AtlasDataDest
    type(BC_Atlas_Data), POINTER :: AtlasDataSource
    integer, POINTER, dimension(:  )  :: Cells_Dest

    AtlasDataDest   => BC_Get_Data(Dest)
    AtlasDataSource => BC_Get_Data(Source)
    Cells_Dest      => BC_Get_Cell(Dest)
    ! Make sure this is initialized, so we can find out how many items get put into each atlas
    Cells_Dest      = BC_INVALID_CELL

    call REDISTRIBUTE(AtlasDataDest, AtlasDataSource)

    ! Need to set AtlasDataSize in 
    call BC_Set_Data_Size(Dest, SIZE=COUNT(Cells_Dest /= BC_INVALID_CELL))

    ! And need to set the scope of the cell pointers.
    call SET_SCOPE(Dest, GET_SCOPE(Source))
    
  end subroutine AtlasRedistribute

  subroutine ComputeIndex(Atlas)
    ! Compute the length and offset fields based on the atlas data.
    type(BC_Atlas), intent(INOUT) :: Atlas

    ! Local variables
    integer :: numcharts, datasize, Current_cell, Current_Face, Item, ChartStart
    logical :: NewCell, NewFace, NewChart, ChartComplete
    integer, POINTER, dimension(:) :: Cells
    integer, POINTER, dimension(:) :: Faces
    integer, POINTER, dimension(:) :: ValueIndex
    logical, POINTER, dimension(:) :: UseFunction
    real(r8), POINTER, dimension(:,:) :: Values
    real(r8), POINTER, dimension(:,:) :: Positions
    type(BC_Atlas) :: Temp_Atlas

    ! We plan to totally clobber Atlas, so we need to store the data someplace
    call CLONE(Temp_Atlas, Atlas)

    ! Now we can work on Temp_Atlas, so free up Atlas until we are ready to re-allocate it.
    call FREE(Atlas)

    ! How much data do we have?
    datasize  = SIZE(Temp_Atlas)
    
    ! This routine may change the number of charts in the atlas, since
    ! it may merge charts which were inserted as distinct entries.

    Cells      => BC_Get_Cell(Temp_Atlas)
    Faces      => BC_Get_Face(Temp_Atlas)
    ValueIndex => BC_Get_ValueIndex(Temp_Atlas)
    UseFunction=> BC_Get_UseFunction(Temp_Atlas)
    Values     => BC_Get_Values(Temp_Atlas)
    Positions  => BC_Get_Positions(Temp_Atlas)

    ! Count the number of charts
    numcharts = 0

    current_cell = BC_INVALID_CELL
    current_face = BC_INVALID_FACE

    do Item = 1, datasize
       ! A transition from one chart to the next occurs if either:
       ! A) Cell number changes
       ! B) Face number changes
       NewCell = (Cells(Item) /= Current_Cell)
       NewFace = (Faces(Item) /= Current_Face)
       NewChart = NewCell .OR. NewFace
       if (NewChart) then 
          numcharts = numcharts + 1
          Current_Face = Faces(Item)
          Current_Cell = Cells(Item)
       end if
          
    end do

    ! Allocate Atlas to have enough space
    call ALLOC(Atlas, MaxNumberOfCharts   = NumCharts, &
                      MaxAtlasDataSize    = DataSize,  &
                      DIMENSIONALITY = DIMENSIONALITY(Temp_Atlas),&
                      DOF = BC_Get_DOF(Temp_Atlas))
    ! The scope for the cells doesn't change in this routine
    call SET_SCOPE(Atlas, SCOPE = GET_SCOPE(Temp_Atlas))
    call BC_Set_Name(Atlas, BC_Get_Name(Temp_Atlas))

    ! Put data into Atlas one chart at a time
    ChartStart =  1
    do Item = 1, DataSize
       Current_Cell = Cells(Item)
       Current_Face = Faces(Item)
       ! A transition from one chart to the next occurs if either:
       ! A) Cell number changes
       ! B) Face number changes
       if (Item < DataSize) then
          NewCell = (Cells(Item+1) /= Current_Cell)
          NewFace = (Faces(Item+1) /= Current_Face)
          ChartComplete = NewCell .OR. NewFace
       else
          ! End of data means chart must be complete
          ChartComplete = .TRUE.
       end if

       if (ChartComplete) then
          ! Then put data from complete chart into Atlas.
          call APPEND(Atlas, VALUES    = Values(:,ChartStart:Item),    &
                             POSITIONS = Positions(:,ChartStart:Item), &
                             VALUEINDEX= ValueIndex(ChartStart:Item),  &
                             USEFUNCTION=UseFunction(ChartStart:Item), &
                             CELL      = Current_Cell,                 & 
                             FACE      = Current_Face)
          ChartStart = Item + 1
       end if

    end do

    ! Finally, done with Temp_Atlas, so free it
    call FREE(Temp_Atlas)

  end subroutine ComputeIndex
    
    
  subroutine AtlasRenumberCells(Atlas, SCOPE)
    ! Renumber the cells based on the the input scope
    use legacy_mesh_api, only: ncells
    use pgslib_module, ONLY: PGSLib_SUM_PREFIX
    type(BC_Atlas), intent(INOUT), target :: Atlas
    integer, intent(IN) :: SCOPE

    ! Local variables
    integer :: c
    integer, dimension(ncells) :: Global_Cell_Number
    integer :: Global_Offset
    integer, POINTER, dimension(:) :: Cells

    ! If output scope should be local, then convert global to local,
    ! unless scope is already local
    if (SCOPE == BC_Cells_Scope_Local) then
       if (GET_SCOPE(Atlas) == BC_Cells_Scope_Global) then
          ! Current atlas scope is global, so there is something to do

          ! Need the global cell number to add on or subtract off
          ! global cell number = local cell number + global offset
          Global_Cell_Number = PGSLib_SUM_PREFIX( (/ (1, c = 1, ncells) /) )
          Global_Offset      = Global_Cell_Number(1) - 1

          Cells      => BC_Get_Cell(Atlas)
          ! Subtract that offset from all the cell numbers to get local cell number
          Cells = Cells - Global_Offset
          ! Set the scope to be local
          call Set_Scope(Atlas, SCOPE = BC_Cells_Scope_Local)
       end if
    end if
    
    ! If output scope should be global, then convert local to global,
    ! unless scope is already global
    if (SCOPE == BC_Cells_Scope_Global) then
       if (GET_SCOPE(Atlas) == BC_Cells_Scope_Local) then
          ! Current atlas scope is local, so there is something to do

          ! Need the global cell number to add on or subtract off
          ! global cell number = local cell number + global offset
          Global_Cell_Number = PGSLib_SUM_PREFIX( (/ (1, c = 1, ncells) /) )
          Global_Offset      = Global_Cell_Number(1) - 1

          Cells      => BC_Get_Cell(Atlas)
          ! Add that offset to all the cell numbers to get global cell number
          Cells = Cells + Global_Offset
          ! Set the scope to be global
          call Set_Scope(Atlas, SCOPE = BC_Cells_Scope_Global)
       end if
    end if
    
  end subroutine AtlasRenumberCells

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!   Output routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AtlasCollate(collated_atlas, local_atlas)
    ! Collated  Local_Atlas into Collated_Atlas.
    ! This routine assumes that some fields of Collated_Atlas have
    ! already been setup, so Collated_Atlas is an INOUT argument.
    use parallel_util_module, only: Is_IO_PE
    use pgslib_module,        ONLY: PGSLib_Global_SUM, PGSlib_Collate
    type(BC_Atlas), intent(INOUT) :: collated_atlas
    type(BC_Atlas), intent(IN) :: local_atlas

    ! local variables
    integer :: d, p, Collated_Size
    integer,  pointer, dimension(:) :: Collated_Cells, Local_Cells
    integer,  pointer, dimension(:) :: Collated_Faces, Local_Faces
    integer,  pointer, dimension(:) :: Collated_ValueIndex, Local_ValueIndex
    logical,  pointer, dimension(:) :: Collated_UseFunction, Local_UseFunction
    real(r8), pointer, dimension(:,:) :: Collated_Values, Local_Values
    real(r8), pointer, dimension(:,:) :: Collated_Positions, Local_Positions
    integer :: DIMS
    real(r8), ALLOCATABLE, dimension(:,:) :: ChartPositions

    ! Would like to have this on the stack, but Fujitsu doesn't like
    ! using DIMENSIONALITY in the declaration list
    DIMS = DIMENSIONALITY(Local_Atlas)
    ALLOCATE(ChartPositions(DIMS, 1))

    !!! THIS CODE IS BROKEN.  I ASSUME ONE CELL PER CHART WHICH IS BOGUS!!!!

    Collated_Size = PGSLib_Global_SUM(SIZE(Local_Atlas))
    if (.NOT. Is_IO_PE()) then
       Collated_Size = 0
    end if

    ALLOCATE(Collated_Cells(Collated_Size))
    ALLOCATE(Collated_Faces(Collated_Size))
    ALLOCATE(Collated_ValueIndex(Collated_Size))
    ALLOCATE(Collated_UseFunction(Collated_Size))
    ALLOCATE(Collated_Values(BC_Get_DOF(Local_Atlas),Collated_Size))
    ALLOCATE(Collated_Positions(DIMENSIONALITY(Collated_Atlas),Collated_Size))

    ! Collate the data onto the IO processor
    Local_Cells    => BC_Get_Cell(Local_Atlas)
    call pgslib_collate(Collated_Cells, Local_Cells)

    Local_Faces    => BC_Get_Face(Local_Atlas)
    call pgslib_collate(Collated_Faces, Local_Faces)

    Local_ValueIndex    => BC_Get_ValueIndex(Local_Atlas)
    call pgslib_collate(Collated_ValueIndex, Local_ValueIndex)

    Local_UseFunction    => BC_Get_UseFunction(Local_Atlas)
    call pgslib_collate(Collated_UseFunction, Local_UseFunction)

    Local_Values    => BC_Get_Values(Local_Atlas)
    do d = 1, BC_Get_DOF(Local_Atlas)
       call pgslib_collate(Collated_Values(d,:), Local_Values(d,:))
    end do

    Local_Positions    => BC_Get_Positions(Local_Atlas)
    do d = 1, SIZE(Local_Positions, 1)
       call pgslib_collate(Collated_Positions(d,:), Local_Positions(d,:))
    end do

    ! Now append it into the collated atlas
    if (Is_IO_PE()) then
       do p = 1, Collated_Size
          call APPEND(ATLAS      = Collated_Atlas,             &
                      VALUES     = Collated_Values(:,p:p), &
                      POSITIONS  = Collated_Positions(:,p:p),&
                      ValueIndex = Collated_ValueIndex(p:p),&
                      UseFunction= Collated_UseFunction(p:p),&
                      Cell       = Collated_Cells(p),&
                      Face       = Collated_Faces(p))
       end do
    end if

    call BC_Set_Name(Collated_Atlas, BC_Get_Name(Local_Atlas))
    DEALLOCATE(Collated_Positions)
    DEALLOCATE(Collated_ValueIndex)
    DEALLOCATE(Collated_UseFunction)
    DEALLOCATE(Collated_Values)
    DEALLOCATE(Collated_Faces)
    DEALLOCATE(Collated_Cells)

    DEALLOCATE(ChartPositions)

  end subroutine AtlasCollate

end Module BC_ATLASES_Util
