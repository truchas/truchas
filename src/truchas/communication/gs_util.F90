MODULE GS_UTIL
  !=======================================================================
  ! Purpose(s):
  !
  !   Define utilities used to support gather & scatter.
  !
  !=======================================================================
  use gs_info_module, only: EE_TRACE,               &
                            EE_All_Ngbr_Trace,      &
                            EE_Mask_Initialized,    &
                            EL_Nbr_Mask,            &
                            NN_All_Ngbr_Trace
  implicit none

  ! Private Module
  private

  ! Public procedures
  Public :: GS_INIT_EE_MASK, &
            EE_GS_INIT,      &
            EN_GS_INIT,      &
            NN_GS_INIT

  ! Arrays and variables used only inside this module

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE GS_INIT_EE_MASK()
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use truchas_logging_services, only: TLS_fatal_if_any
    use kind_module,      only: int_kind
    use mesh_module,      only: Mesh, DEGENERATE_FACE
    use parameter_module, only: ncells, nfc

    implicit none 

    ! Local variables
    integer(KIND = int_kind) :: memerror, f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (El_Nbr_MASK(nfc, ncells), STAT=memerror)
    call TLS_fatal_if_any (memerror /= 0, 'EE_GATHER_INT: could not allocate El_Nbr_Mask')
    
    do f = 1,nfc
       El_Nbr_Mask(f,:) = Mesh%Ngbr_cell(f) /= 0 .and. &
                          Mesh%Ngbr_cell(f) /= DEGENERATE_FACE
    end do

    return

  END SUBROUTINE GS_INIT_EE_MASK

  SUBROUTINE EE_GS_INIT()
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use ArrayAllocate_Module
    use kind_module,      only: int_kind
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nfc
    use two_level_partition,    only: Cell_Two_Level_Partitioning, &
                                      Cell_Cell_Two_Level_Partitioned, &
                                      Set_Of_Cell_Cell_Edges, &
                                      Cell_Cell_Two_Level_Edges_Part
    use var_vector_module
    use pgslib_module,    only: PGSLIB_SETUP_TRACE,           &
                                      GRAPH_HEAD,             &
                                      GRAPH_TAIL,             &
                                      INITIALIZE,             &
                                      SET, &
                                      PGSLib_Sum_Prefix, &
                                      ALLOC, &
                                      LOOKUP_TAIL, &
                                      Get_Num_Edges_Available,&
                                      PGSLib_Local, &
                                      PGSLib_Global_SUM
    
    implicit none

    ! Local variables
    integer(KIND = int_kind), dimension(nfc,ncells) :: Mesh_Ngbr_Cell
    integer(KIND = int_kind), dimension(:,:),       &
                              POINTER               :: Mesh_Ngbr_Cell_PE
    integer(KIND = int_kind), dimension(:),         &
                              POINTER               :: Mesh_Ngbr_Cells_All
    integer(KIND = int_kind), dimension(:),         &
                              POINTER               :: Tail_Partition
    ! This stuff is for setting up the two_level partitioned graph
    integer(int_kind),       dimension(2, ncells*nfc)  :: Cell_Cell_Edges
    integer(int_kind),       dimension(ncells)         :: Global_Cell_Number
    integer(int_kind)                                  :: edge
    integer(int_kind)                                  :: c
    integer(int_kind)                                  :: f
    integer(int_kind)                                  :: Number_Edges_Avail
    integer(int_kind)                                  :: Number_Edges

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.NOT. EE_MASK_Initialized) then
       EE_MASK_Initialized = .true.
       call GS_INIT_EE_MASK()
    end if
    do f = 1,nfc
       Mesh_Ngbr_Cell(f,:) = Mesh%Ngbr_cell(f)
    end do

    ! We want pgslib to return values in Mesh_Ngbr_Cell_PE.  To signal that,
    ! we send it an unallocated pointer.  Since pointers are not NULL by defualt
    ! (in F90, fixed in F95), we must NULLIFY.
    NULLIFY(Mesh_Ngbr_Cell_PE)

    EE_Trace => PGSLib_Setup_Trace(INDEX        = Mesh_Ngbr_Cell,   &
                                   SIZE_OF_DEST = ncells,           &
                                   PE_ARRAY     = Mesh_Ngbr_Cell_PE,&
                                   MASK         = El_Nbr_Mask)
    do f = 1, nfc
       ! Save original Ngbr_Cell
       Mesh%Ngbr_Cell_Orig(f) = Mesh%Ngbr_Cell(f) 
       ! Extract the new index numbers
       Mesh%Ngbr_Cell(f)      = Mesh_Ngbr_Cell(f,:) 
       ! Save the PE number, for global/local stuff
       Mesh%Ngbr_Cell_PE(f)   = MERGE(Mesh_Ngbr_Cell_PE(f,:), -1, EL_Nbr_Mask(f,:))
    end do

    ! This is also an excellent place to establish the cell-cell (across faces) partitioned graph.
    CALL INITIALIZE(Cell_Cell_Two_Level_Partitioned)
    Global_Cell_Number = PGSLib_Sum_Prefix( (/ (1, c=1,ncells) /) )
    edge = 0
    do c = 1, ncells
       do f = 1, nfc
          edge = edge + 1
          ! Graph points from tail to head.
          Cell_Cell_Edges(GRAPH_TAIL, edge) = Global_Cell_Number(c)
          Cell_Cell_Edges(GRAPH_HEAD, edge) = Mesh(c)%Ngbr_Cell_Orig(f)
       end do
    end do
    
    call SET(Cell_Cell_Two_Level_Partitioned, &
            NUMBER_AVAILABLE_EDGES=SIZE(Cell_Cell_Edges,2), &
            EDGES = Cell_Cell_Edges, &
            HEAD_PARTITIONED_SET = Cell_Two_Level_Partitioning, &
            TAIL_PARTITIONED_SET = Cell_Two_Level_Partitioning)


    ! We also need to know how the edges are partitioned.  The partitioning
    ! of edges is according to the partition of the tail of each vertex.
    ! We need a set of edges.
    call INITIALIZE(Set_Of_Cell_Cell_Edges)
    Number_Edges_Avail = Get_Num_Edges_Available(Cell_Cell_Two_Level_Partitioned)
    Number_Edges = PGSlib_Global_Sum(Number_Edges_Avail)
    call ALLOC(Set_Of_Cell_Cell_Edges, Number_Edges)

    ! Use the partitions of the tails to partition the edges
    ALLOCATE(Tail_Partition(Number_Edges_Avail))
    do edge = 1, Number_Edges_Avail
       Tail_Partition(edge) = LOOKUP_TAIL(Cell_Cell_Two_Level_Partitioned, edge)
    end do
       

    call INITIALIZE(Cell_Cell_Two_Level_Edges_Part)
    call SET(Cell_Cell_Two_Level_Edges_Part, Set_Of_Cell_Cell_Edges, &
             Tail_Partition, SCOPE=PGSLib_LOCAL)
    DEALLOCATE(Tail_Partition)

    ! This code is for gathering/scattering from/to all neighbor cells
    ! Notice that here we don't need a mask, since the length of the
    ! neighbor list varies.  
    
    ! Preserve the original global cell numbers
    Call CREATE(Mesh%Ngbr_Cells_All_Orig, SIZES = SIZES(Mesh%Ngbr_Cells_All))
    Mesh%Ngbr_Cells_All_Orig = FLATTEN(Mesh%Ngbr_Cells_All)

    ! We need to move the list of neighbors into a single long array
    ! in order to pass it to pgslib for initialization.
    ! That is provided by FLATTENing Mesh%Ngbr_Cells_All
    Mesh_Ngbr_Cells_All => FLATTEN(Mesh%Ngbr_Cells_All)
    
    EE_All_Ngbr_Trace => PGSLib_Setup_Trace(INDEX = Mesh_Ngbr_Cells_All, &
                                            SIZE_OF_DEST = ncells)
    
!!$    Don't need this, because we are already pointing at permanent storage slot.
!!$    ! Move renumbered index back to permanent storage
!!$    Mesh%Ngbr_Cells_All = Mesh_Ngbr_Cells_ALL


    return

  END SUBROUTINE EE_GS_INIT

  SUBROUTINE EN_GS_INIT()
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use gs_info_module,   only: EN_TRACE
    use kind_module,      only: int_kind
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nnodes, nvc
    use pgslib_module,    only: PGSLIB_SETUP_TRACE


    implicit none

    ! Local variables
    integer(KIND = int_kind) :: v
    integer(KIND = int_kind), dimension(nvc, ncells) :: Mesh_Ngbr_Vrtx
    integer(KIND = int_kind), dimension(:,:) ,       &
                              POINTER                :: Mesh_Ngbr_Vrtx_PE

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do v = 1, nvc
       Mesh_Ngbr_Vrtx(v,:) = Mesh%Ngbr_Vrtx(v)
    end do

    ! We want pgslib to return values in Mesh_Ngbr_Vrtx_PE.  To signal that,
    ! we send it an unallocated pointer.  Since pointers are not NULL by defualt
    ! (in F90, fixed in F95), we must NULLIFY.
    NULLIFY(Mesh_Ngbr_Vrtx_PE)

    EN_Trace => PGSLib_Setup_Trace(INDEX        = Mesh_Ngbr_Vrtx,    &
                                   SIZE_OF_DEST = NNodes,            &
                                   PE_ARRAY     = Mesh_Ngbr_Vrtx_PE)
                                
    do v = 1,nvc
       ! Save original vertex number in case we need it
       Mesh%Ngbr_Vrtx_Orig(v) = Mesh%Ngbr_Vrtx(v)
       ! Extract the new vertex number into our data structure
       Mesh%Ngbr_Vrtx(v)      = Mesh_Ngbr_Vrtx(v,:)
       ! Save the PE number, for global/local stuff
       Mesh%Ngbr_Vrtx_PE(v)   = Mesh_Ngbr_Vrtx_PE(v,:)
    end do

    return
  
  END SUBROUTINE EN_GS_INIT

  SUBROUTINE NN_GS_INIT()
    !=======================================================================
    ! Purpose(s): Initialize the trace for node-node gather & scatter
    !
    !=======================================================================
    use ArrayAllocate_Module
    use kind_module,      only: int_kind
    use mesh_module,      only: Vertex_Ngbr_All, Vertex_Ngbr_All_Orig
    use parameter_module, only: nnodes
    use var_vector_module
    use pgslib_module,    only: PGSLIB_SETUP_TRACE
    
    implicit none

    ! Local variables
    integer(KIND = int_kind), dimension(:),         &
                              POINTER               :: Vertex_Ngbr_List

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! This code is for gathering/scattering from/to all neighbor vertices
    ! Notice that here we don't need a mask, since the length of the
    ! neighbor list varies.  
    
    ! Preserve the original global cell numbers
    ALLOCATE(Vertex_Ngbr_All_Orig(SIZE(Vertex_Ngbr_All)))
    CALL CREATE(Vertex_Ngbr_All_Orig, SIZES = SIZES(Vertex_Ngbr_All))
    Vertex_Ngbr_All_Orig = FLATTEN(Vertex_Ngbr_All)

    ! We need to move the list of neighbors into a single long array
    ! in order to pass it to pgslib for initialization.
    Vertex_Ngbr_List => FLATTEN(Vertex_Ngbr_All)
    
    NN_All_Ngbr_Trace => PGSLib_Setup_Trace(INDEX = Vertex_Ngbr_List, &
                                            SIZE_OF_DEST = nnodes)
    
    return

  END SUBROUTINE NN_GS_INIT

END MODULE GS_UTIL
