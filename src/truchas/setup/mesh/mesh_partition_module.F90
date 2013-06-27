Module MESH_PARTITION_MODULE
  !======================================================================
  ! Purpose(s):
  !
  !   Partition a mesh for parallel computation.  The routines in this
  !   module determine which partition each node and element belongs to,
  !   determine how many elements and how many nodes on each PE, and
  !   provide a permutation vector to map original order to partitioned
  !   order for both mesh and vertex. 
  !
  ! Contains: MESH_PARTITIONS
  !
  ! Author(s) : Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  use truchas_logging_services
  implicit none
  private

  public :: MESH_PARTITIONS

  ! Interface blocks
  Interface
     Subroutine CHACO_F90_WRAPPER (nnodes_tot, nPE, start, edges, node_colors_tot, status)
       Implicit None
       Integer               :: nnodes_tot, nPE, status
       Integer, Dimension(*) :: start, edges, node_colors_tot
     End Subroutine CHACO_F90_WRAPPER
  End Interface

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE MESH_PARTITIONS (MeshPermute, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Return new values for nnodes and ncells, which are the values
    !   when mesh is partitioned.  MeshPermute and VertexPermute, which
    !   tell how to take input data and permute it into partitioned data.  
    !
    !=======================================================================
    use mesh_gen_data,        only: Generated_Mesh
    use partitioner_data,     only: get_partitioner, &
                                    PART_None,       &
                                    PART_Cartesian,  &
                                    PART_Chaco,      &
                                    PART_Metis

    ! Arguments
    integer, dimension(:), pointer :: MeshPermute
    integer, dimension(:), pointer :: VertexPermute

    ! Local Variables
    integer, dimension(:),   pointer :: Node_Colors
    integer, dimension(:),   pointer :: Cell_Colors
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call TLS_info ('')
    call TLS_info (' Determining mesh partitions and permutation vectors.')

    ! Do the desired partitioning
    ! After this SELECT clause Cell_Colors and Node_Colors are defined
    SELECT case(get_partitioner())
       CASE(PART_None)
          ! For no partitioning, assign all cells same color.  Same for all nodes
          Cell_Colors => Identity_Cell_Colors()
          Node_Colors => Identity_Node_Colors()

       CASE(PART_Cartesian)
          ! The cartesian partitioner works only for meshes we generated.
          if (.NOT. Generated_Mesh()) then
             call TLS_fatal ('Mesh_Partition: asked for Cartesian partitioning on mesh not generated internally')
          end if
          
          ! For cartesian mesh, intuitively want to color cells
          Cell_Colors => CARTESIAN_PARTITION ()    
          ! Color nodes based on cell color.
          Node_Colors => NODE_COLORS_FROM_CELL_COLORS (Cell_Colors)

       CASE(PART_Chaco)
#ifndef USE_CHACO
          ! Fatal error to ask for Chaco if not built with Chaco
          call TLS_fatal ('Mesh_Partition: chaco requested for paritioner, but code built with use_Chaco = NO')
#endif
          ! General partitioner colors nodes
          node_colors => GENERAL_PARTITION ()

          ! Based on node color, color cells
          Cell_colors => CELL_COLORS_FROM_NODE_COLORS (node_colors)

       CASE(PART_Metis)
#ifndef USE_METIS
          ! Fatal error to ask for Metis if not built with Metis
          call TLS_fatal ('Mesh_Partition: Metis requested for paritioner, but code built with use_Metis = NO')
#endif
          
       CASE DEFAULT
          call TLS_fatal ('Mesh_Partition: Unrecognized partitioner requested.')

    END SELECT

    ! We need to find the number of nodes number of elements which should be assigned
    ! to each processor.
    ! We also need to determine the permutation vectors which take the old
    ! layout to the new layout.
    call COMPUTE_PERMUTATIONS (Cell_Colors, Node_Colors, MeshPermute, VertexPermute)

    ! Deallocate temporaries.
    DEALLOCATE (Cell_Colors, Node_Colors)

    return

  END SUBROUTINE MESH_PARTITIONS

  FUNCTION GENERAL_PARTITION () RESULT( Node_Colors )
    !=======================================================================
    ! Purpose(s):
    !
    !   Given a general mesh (through use association of Mesh and Vertex)
    !   color the nodes.  The result is pointer to an array of node colors.
    !   The array is distributed, with the same distribution as (conformal to)
    !   Vertex.
    !
    !=======================================================================
    use mesh_gen_data,        only: Partitions_Total
    use mesh_utilities,       only: NODE_CONNECTIVITY
    use parallel_info_module, only: p_info
    use parameter_module,     only: nnodes, nnodes_tot
    use pgslib_module,        only: PGSLIB_GLOBAL_SUM, PGSLIB_COLLATE, &
                                    PGSLIB_DIST, PGSLib_SUM_PREFIX
    use var_vector_module
    use string_utilities, only: i_to_c

    ! Arguments
    integer, dimension(:), pointer :: Node_Colors

    ! Local Variables
    integer :: edge, n_edges, n_edges_tot, head, thisnode, node
    integer, dimension(:), pointer  :: Edges_Tot1
    integer, dimension(:), pointer  :: Edges_Tot2
    integer, dimension(:), pointer  :: Start
    integer, dimension(:), pointer  :: Node_Colors_Tot
    integer, dimension(:,:), pointer :: Edges
    type (int_var_vector), dimension(nnodes)   :: All_Edges
    integer :: status
    integer :: StartEdge, EndEdge
    integer, dimension(1) :: GlobalStartNode
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    ! Need some arrays for coloring
    if (p_info%IOP) then
       ALLOCATE(Node_Colors_Tot(nnodes_tot))
    else
       ALLOCATE(Node_Colors_Tot(0))
    end if

    ! We may get here even if we only want a single partition.  In that case, skip
    ! the whole Chaco thing and return a single color.
    call TLS_info (' Using ' // i_to_c(Partitions_Total) // ' partitions.')

    CALL_PARTITIONER: if (Partitions_Total == 1) then
      
       Node_Colors_Tot = 0

    ELSE
       ! The partitioner wants node<->node connectivity.  We don't have that, so
       ! we have to construct it from the cell->node connectivity, which we do have.

       All_Edges = NODE_CONNECTIVITY()
       n_edges   = SUM(SIZES(All_Edges))

       ALLOCATE(Edges(2, n_edges))
       Edges(2,:) = FLATTEN(All_Edges)
       
       ! Need to put in global node numbers
       GlobalStartNode = PGSLib_SUM_PREFIX((/ nnodes /))
       GlobalStartNode = GlobalStartNode - nnodes + 1

       StartEdge = 1
       do node = 1, nnodes
          EndEdge = StartEdge + SIZES(All_Edges(node)) - 1
          Edges(1,StartEdge:EndEdge) = node + GlobalStartNode(1)
          StartEdge = EndEdge + 1
       end do

       ! Done with All_Edges
       call DESTROY(All_Edges)

       ! For Chaco, we need to pass in a single array on the IO processor,
       ! since it is a serial code.
       n_edges_tot = PGSLib_GLOBAL_SUM (SIZE(Edges,2))
       if (p_info%IOP) then
          ALLOCATE (Edges_Tot1(n_edges_tot), Edges_Tot2(n_edges_tot))
       else
          ALLOCATE (Edges_Tot1(0), Edges_Tot2(0))
       endif
       call PGSLib_COLLATE (Edges_Tot1(:), Edges(1,:))
       call PGSLib_COLLATE (Edges_Tot2(:), Edges(2,:))

       DEALLOCATE (Edges)

       ! Prepare to pass data to Chaco
       SERIAL_CODE: if (p_info%IOP) then

          ALLOCATE(Start(nnodes_tot+1))
          ! Chaco assumes that the edges are listed in order.
          ! Chaco wants a single array containing all the tails of the edges,
          ! and another array, with one entry per node, indexing into
          ! the first edge for each node.

          Head = -1
          ThisNode = 0
          do edge = 1, SIZE(Edges_Tot1,1)
             ! If this edge has the same head (originating node) as the 
             ! previous edge, then we don't do anything.
             ! If this edge has a new head, then we need to notice that,
             ! and process accordingly.
             if (Edges_Tot1(Edge) /= Head) then
                if (Edges_Tot1(Edge) <= Head) then
                   write(message,*) 'MESH_PARTITION: Edge, Edges_Tot1(Edge), Head = ', &
                        Edge, Edges_Tot1(Edge), Head
                   call TLS_panic (message)
                endif
                ThisNode = ThisNode + 1
                ! Since we are passing to C, start has to be 0 based to index into Edges_Tot2.
                Start(ThisNode) = Edge - 1 
                Head = Edges_Tot1(Edge)
             end if
          end do
       
          ! We should have found a segment for each node (we are assuming a connected graph).
          if (ThisNode < (nnodes_tot-1)) then
             call TLS_panic ('Mesh_Partition: not enough nodes in Mesh_Partition_Module')
          else
             if (ThisNode == (nnodes_tot-1)) then
                ThisNode = ThisNode+1
                Start(ThisNode) = SIZE(Edges_Tot2,1)
             end if
          endif
          Start(ThisNode+1) = SIZE(Edges_Tot2,1)

          ! If built with USE_CHACO, this calls Chaco, others
          ! returns 0 in Node_Colors_Tot (which means a single color).
          call TLS_info (' Partitioning with Chaco.')
          Call CHACO_F90_WRAPPER (nnodes_tot,      &
                                  Partitions_Total,&
                                  Start,           &
                                  Edges_Tot2,      &
                                  Node_Colors_Tot, &
                                  status)
          DEALLOCATE(Start)
       else
          status = 0
       end if SERIAL_CODE

       Call TLS_fatal_if_any ((status /= 0), 'MESH_PARTITION: CHACO_F90_WRAPPER failure: error code 1')
       DEALLOCATE (Edges_Tot1, Edges_Tot2)

    end if CALL_PARTITIONER

    ALLOCATE (Node_Colors(nnodes))

    call PGSLib_DIST (Node_Colors, Node_Colors_Tot)

    DEALLOCATE(Node_Colors_Tot)

  END FUNCTION GENERAL_PARTITION

  FUNCTION CARTESIAN_PARTITION () RESULT( Cell_Colors )
    !=======================================================================
    ! Purpose(s):
    !
    !   Given a cartesian mesh (generated on-the-fly only?), 
    !   color the Cells.  The result is pointer to an array of cell colors.
    !   The array is distributed, with the same distribution as (conformal to)
    !   Mesh.
    !
    !=======================================================================
    use partitioner_data,     only: get_Processor_Array
    use parameter_module,     only: ncells, Nx_Tot, ndim
    use pgslib_module,        only: PGSLIB_SUM_PREFIX

    ! Arguments
    integer, dimension(:), pointer :: Cell_Colors

    ! Local Variables
    logical                          :: fatal
    integer                          :: c, d, d1
    integer, dimension(ndim, ncells) :: BlockNumbers, Indices
    integer, dimension(ncells)       :: GlobalIndex
    integer, dimension(ndim)         :: Strides, PE_Strides, BlockSizes
    integer,        dimension(ndim)         :: Processor_Array
    character(128) :: message
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Find the cartesian indices for each cell
    GlobalIndex = PGSLib_SUM_PREFIX ( (/ (1, c=1,ncells) /) ) - 1
    Strides = 1
    do d = 2, ndim
       Strides(d) = Strides(d-1)* Nx_Tot(d-1)
    end do
    
    Indices(ndim, :) = GlobalIndex / Strides(ndim)
    do d = ndim-1,1,-1
       Indices(d,:) = GlobalIndex
       do d1 = d+1,ndim
          Indices(d,:) = Indices(d,:) - Indices(d1,:)*Strides(d1)
       end do
       Indices(d,:) = Indices(d,:)/Strides(d)
    end do
    
    ! Convert this to cartesian block number
    ! First compute block sizes
    ! Need the processor_array for that
    Processor_Array = get_Processor_Array()
    fatal = .false.
    do d = 1, ndim
       BlockSizes(d) = Nx_Tot(d)/ Processor_Array(d)
       ! Check that mesh and Processor_Array are commensurate
       if (Processor_Array(d) * BlockSizes(d) /= Nx_Tot(d) ) then
          write(message,*) 'Processor array not commensurate with mesh: Dimension, Nx_Tot, Processor_Array, BlockSizes', &
                                  d, Nx_Tot(d), Processor_Array(d), BlockSizes(d)
          call TLS_error (message)
          fatal = .true.
       end if
    end do
    call TLS_fatal_if_any (fatal, 'Cartesian_Partition: Processor array not commensurate with mesh')

    do d = 1, ndim
       BlockNumbers(d,:) = Indices(d,:)/BlockSizes(d)
    end do
    
    ! Finally, cell color is just total processor number (assuming 
    ! *FORTRAN* (is that row or column major?) ordering of processors).
    ! Need to use Processor_Array strides
    PE_Strides = 1
    do d = 2, ndim
       PE_Strides(d) = PE_Strides(d-1) * Processor_Array(d-1)
    end do
    
    ALLOCATE (Cell_Colors(ncells))
    Cell_Colors = 0
    do d = 1, ndim
       Cell_Colors = Cell_Colors + BlockNumbers(d,:) * PE_Strides(d)
    end do

  END FUNCTION CARTESIAN_PARTITION

  FUNCTION CELL_COLORS_From_NODE_COLORS (Node_Colors) RESULT (Cell_Colors)
    !=======================================================================
    ! Purpose(s):
    !
    !   Given a general mesh (through use association of Mesh and Vertex)
    !   and a coloring for each node, assign a coloring to the cells (elements).
    !
    !=======================================================================
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nvc
    use pgslib_module,    only: PGSLib_GATHER

    ! Arguments
    integer, dimension(:), pointer :: CELL_COLORS
    integer, dimension(:)          :: Node_Colors

    ! Local Variables
    integer :: node
    integer, dimension(nvc,ncells) :: Cell_Colors_All, Mesh_Cell_Nodes

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! We need space for the result.
    ALLOCATE (Cell_Colors(ncells))

    ! Now we have a partition (color) assignment for each node (vertex)
    ! We will assign element colors based on the color of surrounding nodes
    do node = 1, nvc
       Mesh_Cell_Nodes(node,:) = Mesh%Ngbr_Vrtx(node)
    end do

    Cell_Colors_ALL = -1
    call PGSLib_GATHER (DEST   = Cell_Colors_All, &
                        SOURCE = Node_Colors,     &
                        INDEX  = Mesh_Cell_Nodes)

    Cell_Colors = MAXVAL (Cell_Colors_All, DIM = 1)

  END FUNCTION CELL_COLORS_From_NODE_COLORS
  
  FUNCTION NODE_COLORS_From_Cell_Colors (Cell_Colors) RESULT( Node_Colors )
    !=======================================================================
    ! Purpose(s):
    !
    !   Given a general mesh (through use association of Mesh and Vertex)
    !   and a coloring for each cell, assign a coloring to the nodes (vertices).
    !
    !=======================================================================
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nnodes, nvc
    use pgslib_module,    only: PGSLib_SCATTER_MAX

    ! Arguments
    integer, dimension(:), pointer :: Node_Colors
    integer, dimension(:)          :: Cell_Colors

    ! Local Variables
    integer :: node
    integer, dimension(nvc,ncells) :: Mesh_Cell_Nodes, Temp_Cell_Colors

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! We need space for the result
    ALLOCATE (Node_Colors(nnodes))

    ! Now we have a partition (color) assignment for each cell (elements)
    ! We will assign vertex colors based on the color of the connected cells.
    do node = 1, nvc
       Mesh_Cell_Nodes(node,:)  = Mesh%Ngbr_Vrtx(node)
       Temp_Cell_Colors(node,:) = Cell_Colors
    end do
    Node_Colors = 0

    call PGSLib_SCATTER_MAX (DEST   = Node_Colors,      &
                             SOURCE = Temp_Cell_Colors, &
                             INDEX  = Mesh_Cell_Nodes)

  END FUNCTION NODE_COLORS_From_Cell_Colors
  
  SUBROUTINE COMPUTE_PERMUTATIONS (Cell_Colors, Node_Colors, MeshPermute, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Input: Cell_Colors
    !          Node_Colors
    !
    !   Return:
    !          MeshPermute and VertexPermute, which tell how 
    !          to take input data and permute it into partitioned data.  
    !          Their sizes are the new ncells and nnodes.
    !   
    !
    !=======================================================================
    use parallel_util_module, only: p_info, Is_IO_PE
    use two_level_partition,  only: Cell_Two_Level_Partitioning, &
                                    Set_Of_Cells, &
                                    Precond_2level_Active
    use parameter_module,     only: ncells, nnodes, ncells_tot
    use pgslib_module,        only: PGSLib_REDISTRIBUTE, PGSLib_PERMUTE, &
                                    PGSLib_Collate, &
                                    ALLOC, INITIALIZE, SET, PGSLib_Local,&
                                    Get_Num_Partitions_Available,        &
                                    Get_Num_Partitions
    use string_utilities, only: i_to_c

    ! Arguments
    integer, dimension(:), intent(INOUT) :: Node_Colors
    integer, dimension(:), intent(INOUT) :: Cell_Colors
    integer, dimension(:), pointer       :: MeshPermute
    integer, dimension(:), pointer       :: VertexPermute

    ! Local Variables
    integer :: ncells_new, nnodes_new
    integer :: pe
    integer :: Partitions_This_Processor
    integer, dimension(:), POINTER :: Cell_Colors_New
    integer, dimension(ncells)              :: MeshPermute_Orig_Layout
    integer, dimension(nnodes)              :: VertexPermute_Orig_Layout
    integer, dimension(p_info%nPE)          :: num_avail_part_tot
#ifdef USE_OLD_PERMUTE_WAY
    integer :: partition
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! First take care of the cells
    ! Find the permutation and distribution
    call PERMUTE_Colors_To_Processors(PERMUTER = MeshPermute_Orig_Layout, &
                                      ItemsThisProcessor = ncells_new,    &
                                      Colors = Cell_Colors)

    ! Now create the permutation vector with the new layout
    ALLOCATE(MeshPermute(ncells_new))
    call PGSLib_REDISTRIBUTE(DEST   = MeshPermute, &
                             SOURCE = MeshPermute_Orig_Layout)

    ! Finally, before we finish with Cell_Colors, keep track of the partitioning.
    ! Typically this means multiple partitions on each processor.
    ! This got the colors into the right place.

    ALLOCATE(Cell_Colors_New(ncells_new))
    call PGSLib_PERMUTE(DEST   = Cell_Colors_New, &
                        SOURCE = Cell_Colors,     &
                        INDEX  = MeshPermute_Orig_Layout)

    Partitions_This_Processor = MAXVAL(Cell_Colors_New) - MINVAL(Cell_Colors_New) + 1

    ! Establish the set of cells 
    ! This really needs to go someplace else, but put it here since we first need it hear.
    call INITIALIZE(Set_Of_Cells)
    call ALLOC(Set_Of_Cells, ncells_tot)

    ! Now set up the 2-level partitioning
    ! What we've done is split the mesh into many partitions.  Already we've determined
    ! how many partitions go on each processor.  Now we are using the same partitioning
    ! for our 2-level work.  That is, now we are noting that we may have more than
    ! one partition per process.
    call INITIALIZE(Cell_Two_Level_Partitioning)
    call SET(Cell_Two_Level_Partitioning, Set_Of_Cells, Cell_Colors_New, SCOPE=PGSLib_LOCAL)

    if (precond_2level_active) then
       call TLS_info (' Two-Level Partition Parameters ')
       call TLS_info ('    Total of ' // i_to_c(Get_Num_Partitions(Cell_Two_Level_partitioning)) // ' partitions.')
       ! Collate number of partitions/process for output
       call pgslib_collate(num_avail_part_tot, Get_Num_Partitions_Available(Cell_Two_Level_Partitioning))
       if (Is_IO_PE()) then
          do pe = 1, SIZE(num_avail_part_tot)
             call TLS_info ('     On processor ' // i_to_c(pe) // ' partitions.')
          end do
       end if
    end if

    ! Now take care of the nodes
    ! Find permuation and distribution
    call PERMUTE_Colors_To_Processors(PERMUTER = VertexPermute_Orig_Layout, &
                                      ItemsThisProcessor = nnodes_new,      &
                                      COLORS   = Node_Colors)

    ! Now creat the permutaiton vector with the new layout
    ALLOCATE(VertexPermute(nnodes_new))
    call PGSLib_REDISTRIBUTE(DEST   = VertexPermute, &
                             SOURCE = VertexPermute_Orig_Layout)

#ifdef USE_OLD_PERMUTE_WAY


    ! We need to find the number of nodes number of elements which should be assigned
    ! to each processor.
    ! First have to get partitions into compact segments
    ! Mesh still has "old" layout

    ncells_old = SIZE (Cell_Colors,1)
    nnodes_old = SIZE (Node_Colors,1)
    ALLOCATE (Mesh_Color_Rank(ncells_old), Mesh_Temp(ncells_old))
    ! Vertex still has "old" layout
    ALLOCATE (Vertex_Color_Rank(nnodes_old), Vertex_Temp(nnodes_old))

    Mesh_Color_Rank   = PGSLib_GRADE_UP (Cell_Colors)
    ALLOCATE(cell_colors_temp(SIZE(Cell_colors)))
    call PGSLib_PERMUTE (DEST   = Cell_Colors_Temp, &
                         SOURCE = Cell_Colors,      &
                         INDEX  = Mesh_Color_Rank)
    Cell_Colors = Cell_Colors_Temp
    DEALLOCATE(Cell_Colors_Temp)


    Vertex_Color_Rank = PGSLib_GRADE_UP (Node_Colors)

    ALLOCATE(Node_Colors_Temp(SIZE(Node_Colors)))
    call PGSLib_PERMUTE (DEST   = Node_Colors_Temp, &
                         SOURCE = Node_Colors,      &
                         INDEX  = Vertex_Color_Rank)
    Node_Colors = Node_Colors_Temp
    DEALLOCATE(Node_Colors_Temp)


    npartitions = PGSLib_GLOBAL_MAXVAL (Node_Colors) 
    npartitions = npartitions - PGSLib_Global_MINVAL(Node_Colors) + 1

    if (p_info%IOP) then
       ALLOCATE (Node_Colors_Tot(nnodes_tot), Nnodes_Per_Color(npartitions), Nnodes_List(p_info%nPE))
    else
       ALLOCATE (Node_Colors_Tot(1), Nnodes_Per_Color(1), Nnodes_List(1))
    end if

    call PGSLib_COLLATE (Node_Colors_Tot, Node_Colors)

    if (p_info%IOP) then
       Nnodes_Per_Color = 1
       color = 1
       do node = 2, nnodes_tot
          if (Node_Colors_Tot(node) == Node_Colors_Tot(node-1)) then
             Nnodes_Per_Color(color) = Nnodes_Per_Color(color) + 1
          else
             color = color + 1
          end if
       end do
       ! We now know how many nodes of each color, so we can figure
       ! out how to map the colors to the processors.
       partition = 1
       do pe = 1, p_info%nPE
          colors_this_pe = PGSLib_BLOCK_SIZE(npartitions, pe)
          Nnodes_List(pe) = 0
          do color = 1, colors_this_pe
             Nnodes_List(pe) = Nnodes_List(pe) + Nnodes_Per_Color(partition)
             partition = partition + 1
          end do
       end do
    end if
    
    call PGSLib_DIST (nnodes_new, Nnodes_List)

    DEALLOCATE (Nnodes_Per_Color, Node_Colors_Tot, Nnodes_List)

    npartitions = PGSLib_GLOBAL_MAXVAL(Cell_Colors) - PGSLib_Global_MINVAL(Cell_Colors) + 1

    if (p_info%IOP) then
       ALLOCATE(Cell_Colors_Tot(ncells_tot), Ncells_Per_Color(npartitions), Ncells_List(p_info%nPE))
    else
       ALLOCATE (Cell_Colors_Tot(1), Ncells_Per_Color(1), Ncells_List(1))
    end if

    ! Figure out distribution of cell colors to processors
    ! May not have the same number of cell partitions as node partitions
    npartitions = PGSLib_GLOBAL_MAXVAL (Cell_Colors) - PGSLib_Global_MINVAL(Cell_Colors) + 1

    call PGSLib_COLLATE (Cell_Colors_Tot, Cell_Colors)
    Ncells_Per_Color = 1
    color = 1
    if (p_info%IOP) then
       do cell = 2, ncells_tot
          if (Cell_Colors_Tot(cell) == Cell_Colors_Tot(cell-1)) then
             Ncells_Per_Color(color) = Ncells_Per_Color(color) + 1
          else
             color = color + 1
          end if
       end do
       ! We now know how many nodes of each color, so we can figure
       ! out how to map the colors to the processors.
       partition = 1
       do pe = 1, p_info%nPE
          colors_this_pe = PGSLib_BLOCK_SIZE (npartitions, pe)
          Ncells_List(pe) = 0
          do color = 1, colors_this_pe
             Ncells_List(pe) = Ncells_List(pe) + Ncells_Per_Color(partition)
             partition = partition + 1
          end do
       end do
    end if
    
    call PGSLib_DIST (ncells_new, Ncells_List)

    DEALLOCATE (Ncells_Per_Color, Cell_Colors_Tot, Ncells_List)

    ! Finally ncells and nnodes are the appropriate sizes
    ! We now need the permutation vector for the nodes and cells
    ! We generate vectors with the "old" (original) layout
    ! Then we permute those into vectors with the new layout
    ALLOCATE (MeshPermute(ncells_new), VertexPermute(nnodes_new))
    ! Although  these are identity permutation the dest has different layout from the source
    Mesh_Temp = 1

    Mesh_Temp = PGSLib_SUM_PREFIX (Mesh_Temp)
    call PGSLib_PERMUTE (MeshPermute, Mesh_Color_Rank, Mesh_Temp)
    Vertex_Temp = 1
    Vertex_Temp = PGSLib_SUM_PREFIX (Vertex_Temp)
    call PGSLib_PERMUTE (VertexPermute, Vertex_Color_Rank, Vertex_Temp)

    DEALLOCATE (Mesh_Color_Rank, Mesh_Temp, Vertex_Color_Rank, Vertex_Temp)
    
#endif

  END SUBROUTINE COMPUTE_PERMUTATIONS

  SUBROUTINE PERMUTE_Colors_To_Processors (Permuter,                       &
                                           ItemsThisProcessor,             &
                                           Colors)
    !=======================================================================
    ! Purpose(s):
    !
    !   Input: Colors
    !
    !   Output: Permuter
    !           ItemsThisProcessor
    !
    !   Return:
    !          Given the distributed array Colors, figure out
    !          ther permutation array Permuter which puts the Colors
    !          array into order.  Also figure out how many items are
    !          on each processors, and return that in ItermsThisProcessor.
    !   
    !
    !=======================================================================
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLIB_PERMUTE,        &
                                    PGSLIB_GRADE_UP,       &
                                    PGSLIB_GLOBAL_ALL,     &
                                    PGSLIB_GLOBAL_MAXVAL,  &
                                    PGSLIB_GLOBAL_MINVAL,  &
                                    PGSLib_BLOCK_SIZE,     &
                                    PGSLib_SCATTER_SUM,    &
                                    PGSLib_SUM_PREFIX

    ! Arguments
    integer, dimension(:), intent(IN   ) :: Colors
    integer, dimension(:), intent(  OUT) :: Permuter
    integer,               intent(  OUT) :: ItemsThisProcessor

    ! Local Variables
    integer :: NColors, NDomains, color, c, d, ColorsThisDomain
    integer :: Max_Color, Min_Color

    integer, dimension(SIZE(Colors)) :: OrderedColors
    integer, dimension(SIZE(Colors)) :: Domain

    integer, dimension(1)            :: ItemsPerDomain

    integer, dimension(:),   pointer :: DomainForColor

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Sometimes we call this with all items having the same color.
    ! This means only a single partition, and the permutation is the
    ! identity.  Test for and deal with this case first.
    
    ! We have only a single partition if all colors have same value
    Max_Color = PGSLib_Global_MAXVAL(Colors)
    Min_Color = PGSLib_Global_MINVAL(Colors)

    if (PGSLib_Global_ALL(Max_Color == Min_Color)) then
       ItemsThisProcessor = SIZE(Permuter)
       Permuter = PGSLib_SUM_PREFIX ( (/ (1, c=1,SIZE(Permuter)) /) )
       RETURN
    end if
    
    ! Only get here if NOT identity permutation

    ! First put the colors in order
    Permuter = PGSLib_Grade_Up(Colors)
    call PGSLib_Permute(DEST   = OrderedColors, &
                        SOURCE = Colors,        &
                        Index  = Permuter)

    ! Now find out how many colors on each domain (processor)
    ! Construct a table which maps from color to domain.  Note that this
    ! table is the same (replicated) on every processor
    NColors = PGSLib_Global_MAXVAL(OrderedColors) - PGSLib_Global_MINVAL(OrderedColors) + 1
    ALLOCATE(DomainForColor(NColors))

    NDomains = p_info%nPE

    color = 1
    DO d = 1, NDomains
       ColorsThisDomain = PGSLib_Block_Size(NColors, d)
       DO c = 1, ColorsThisDomain
          DomainForColor(color) = d
          Color = Color + 1
       END DO
    END DO

    ! Assign a domain to each color
    DO C = 1, SIZE(OrderedColors)
       Domain(c) = DomainForColor(OrderedColors(c) + 1)
    END DO

    ! By construction the domains are ordered.  Since the sort is stable we
    ! could have assigned domains and then ordered.  I think the result is the same.

    ! Now we count how many items on each domain

    ! By construction we have one domain per processor, so we can 
    ! make an array, ItemsPerDomain of size NDomains, with one item per processor.  Then,
    ! we scatter add 1 from Domain to this ItemsPerDomain.  That is the total
    ! total number of items for the domain.


    ItemsPerDomain = 0
    call PGSLib_SCATTER_SUM(DEST   = ItemsPerDomain,                  &
                            SOURCE = (/ (1, d=1,SIZE(Domain) ) /),  &
                            INDEX  = Domain)
    

    ! Finally, the number of items on this processor has been found
    ItemsThisProcessor = ItemsPerDomain(1)

    ! The desired distribution is obtained by permuting according to
    ! the array PERMUTER, onto a distribution with ItemsThisProcessor
    ! items on each processor.

    ! Clean up and go home
    DEALLOCATE(DomainForColor)
    
  END SUBROUTINE PERMUTE_Colors_To_Processors

  FUNCTION Identity_Cell_Colors () RESULT( Cell_Colors )
    !=======================================================================
    ! Purpose(s):
    !
    !   Return same color for all cells.  Effect is no partitioning
    !
    !   Input: <NONE>
    !
    !   Return:
    !          Cell_Colors is ALLOCATED and assigned the value 0.
    !      
    !=======================================================================
    use parameter_module, only: ncells

    ! Arguments
    integer, dimension(:), pointer :: Cell_Colors

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Allocate the result
    ALLOCATE (Cell_Colors(ncells))
    
    ! Assign it, all the same value for a single partition
    Cell_Colors = 0

  END FUNCTION Identity_Cell_Colors
  
  FUNCTION Identity_Node_Colors () RESULT( Node_Colors )
    !=======================================================================
    ! Purpose(s):
    !
    !   Return same color for all nodes.  Effect is no partitioning
    !
    !   Input: <NONE>
    !
    !   Return:
    !          Node_Colors is ALLOCATED and assigned the value 0.
    !      
    !=======================================================================
    use parameter_module, only: nnodes

    ! Arguments
    integer, dimension(:), pointer :: Node_Colors

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Allocate the result
    ALLOCATE (Node_Colors(nnodes))
    
    ! Assign it, all the same value for a single partition
    Node_Colors = 0

  END FUNCTION Identity_Node_Colors


  SUBROUTINE Identity_Partitons(MeshPermute, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Input: <NONE>
    !
    !   Return:
    !          MeshPermute and VertexPermute, which tell how 
    !          to take input data and permute it into partitioned data.  
    !          In this case they are the identity, and their sizes
    !          are the current ncells & nnodes.
    !   
    !
    !=======================================================================
    use parameter_module, only: ncells, nnodes
    use pgslib_module, only: PGSLib_SUM_PREFIX

    ! Arguments
    integer, dimension(:), pointer :: MeshPermute
    integer, dimension(:), pointer :: VertexPermute

    ! Local Variables
    integer :: c

    ! Allocate the arrays
    ALLOCATE(MeshPermute(ncells))
    ALLOCATE(VertexPermute(nnodes))

    ! Fill them with global counts
    MeshPermute   = PGSLib_SUM_PREFIX ( (/ (1, c=1,ncells) /) ) 
    VertexPermute = PGSLib_SUM_PREFIX ( (/ (1, c=1,nnodes) /) ) 

  END SUBROUTINE Identity_Partitons

END MODULE MESH_PARTITION_MODULE


