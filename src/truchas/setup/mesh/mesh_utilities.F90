!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MESH_UTILITIES
  !=======================================================================
  ! Purpose(s):
  !
  !   Define utilities for the mesh
  !
  !   Public Interface:
  !
  !     * call MESH_DIAGNOSTICS
  !            Run some diagnostics to check mesh integrity
  !        
  !
  ! Contains: 
  !           MESH_DIAGNOSTICS
  !           NODE_CONNECTIVITY
  !
  ! Author(s): Robert Ferrell (ferrell.cpca.com)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Subroutines
  public :: MESH_DIAGNOSTICS
  PUBLIC :: NODE_CONNECTIVITY
  PUBLIC :: SET_DEGENERATE_FACES

CONTAINS

  SUBROUTINE MESH_DIAGNOSTICS ()
    !=======================================================================
    ! Purpose(s):
    !
    !    Compute and report mesh diagnostics. Check for connectivity errors.
    !
    !=======================================================================
    use gs_module,              only: EE_GATHER, EN_SUM_SCATTER
    use mesh_module,            only: Face_Vrtx, Mesh, degenerate_points, &
                                      degenerate_lines, triangle_faces,   &
                                      quad_faces, DEGENERATE_FACE
    use mesh_tests,             only: Test_All_Neighbors
    use mesh_parameter_module,  only: ncells, nnodes, nfc, nvf, ndim, nvc
    use pgslib_module,          only: PGSLib_GLOBAL_COUNT, PGSLib_GLOBAL_MAXVAL,&
                                      PGSLib_SUM_PREFIX
    use mesh_quality_module,    only: TWISTED_CELL_TEST

    ! Local Variables
    integer :: f, n, v1, v2, v3, v4
    logical :: fatal, mismatch
    logical,  dimension(ncells)     :: Match
    integer,  dimension(nfc,ncells) :: Conn, Neighbor
    integer,  dimension(ncells)     :: Global_Cell_Number
    real(r8), dimension(ncells)     :: X1, X2
    real(r8), dimension(nvc,ncells) :: Xn
    real(r8), dimension(nnodes)     :: Coord
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Count and report the total number of hexes, prisms,
    ! pyramids, and tets comprising the mesh. These element
    ! types are defined as follows:
    !  Element Type  Number/Type of Face Degeneracy  Flag
    !  ------------  ------------------------------  ----
    !    hex                     0                     0
    !    prism                 1/line                  1
    !    pyramid               1/point                 3
    !    tet               1/point, 1/line             4

    X2 = 0.0_r8
    do f = 1,nfc
       X1 = 0.0_r8
       do v1 = 1,nvf
          do v2 = 1,nvf
             if (v2 == v1) cycle
             Match = Mesh%Ngbr_Vrtx_Orig(Face_Vrtx(f,v2)) == Mesh%Ngbr_Vrtx_Orig(Face_Vrtx(f,v1))
             ! Vertices match: we have a degeneracy.
             where (Match) X1 = X1 + 1.0_r8
          end do
       end do
       X1 = X1/REAL(nvf)
       ! Add in X1 if this is face is degenerate (line or point).
       if (ndim == 3) then
          X2 = MERGE(X2 + X1, X2, X1 == 1 .or. X1 == 3)
       end if
    end do
    if (ndim == 3) then
       v1 = PGSLib_GLOBAL_COUNT(X2 == 0.0_r8)
       v2 = PGSLib_GLOBAL_COUNT(X2 == 1.0_r8)
       v3 = PGSLib_GLOBAL_COUNT(X2 == 3.0_r8)
       v4 = PGSLib_GLOBAL_COUNT(X2 == 4.0_r8)
       call TLS_info ('')
       call TLS_info ('                               Mesh Diagnostics')
       call TLS_info ('                               ----------------')
       call TLS_info ('')
       call TLS_info ('                   Element Type   Number     Face Type   Number')
       call TLS_info ('                   ------------   ------     ---------   ------')
       write(message,'(22x,"Hex",6x,i8,7x,"Point",4x,i8)') v1, degenerate_points
       call TLS_info (message)
       write(message,'(22x,"Prism",4x,i8,7x,"Line",5x,i8)') v2, degenerate_lines
       call TLS_info (message)
       write(message,'(22x,"Pyramid",2x,i8,7x,"Triangle",1x,i8)') v3, triangle_faces
       call TLS_info (message)
       write(message,'(22x,"Tet",6x,i8,7x,"Quad",5x,i8)') v4, quad_faces
       call TLS_info (message)
    end if

    ! Sum-scatter a value of 1.0 to nodes, making sure
    ! that redundant scatters don't occur at degenerate faces.
    Xn = 1.0_r8
    do v1 = 1,nvc
       Match = Xn(v1,:) == 0.0_r8
       do v2 = 1,nvc
          if (v2 == v1) cycle
          where (Mesh%Ngbr_Vrtx_Orig(v2) == Mesh%Ngbr_Vrtx_Orig(v1) .and. &
                 .not.Match) Xn(v2,:) = 0.0_r8
       end do
    end do
    call EN_SUM_SCATTER (Coord, Xn)
 
    ! Write out unstructured mesh statistics:
    ! the number of vertices having a given
    ! number of cells as neighbors.
    call TLS_info ('')
    call TLS_info ('                            Nodes               Cells')
    call TLS_info ('                            -----               -----')
    v1 = PGSLib_GLOBAL_MAXVAL(Coord)
    MESH_STAT: do f = 1,v1
       n = PGSLib_GLOBAL_COUNT(Coord == f)
       if (n == 0) cycle MESH_STAT
       write (message, 4) n, f
4      format (25x,i8,' are shared by ',1x,i2)
       call TLS_info (message)
    end do MESH_STAT

    ! Make sure connectivity is self-consistent.
    ! Need to use original cell numbers here
    do f = 1,nfc
       Conn(f,:) = MERGE (0, Mesh%Ngbr_Cell_Orig(f), Mesh%Ngbr_Cell_Orig(f) == DEGENERATE_FACE)
    end do
    ! Gather cell connectivity.
    call EE_GATHER (Neighbor, Conn)
    mismatch = .false.
    Global_Cell_Number = PGSLib_SUM_PREFIX((/ (1,n=1,ncells) /))

    do n = 1,ncells
       do f = 1,nfc
          if (Mesh(n)%Ngbr_Cell_Orig(f) == DEGENERATE_FACE .or. &
              Mesh(n)%Ngbr_Cell_Orig(f) == 0) cycle
          mismatch = mismatch .or. Neighbor(f,n) /= Global_Cell_Number(n)
       end do
    end do
    fatal = mismatch
    do f = 1,nfc
       Conn(f,:) = MERGE (0, Mesh%Ngbr_Face(f), Mesh%Ngbr_Face(f) == DEGENERATE_FACE)
    end do
    ! Gather face connectivity.
    call EE_GATHER (Neighbor, Conn)
    mismatch = .false.
    do n = 1,ncells
       do f = 1,nfc
          if (Mesh(n)%Ngbr_Face(f) == DEGENERATE_FACE .or. &
              Mesh(n)%Ngbr_Face(f) == 0) cycle
          mismatch = mismatch .or. Neighbor(f,n) /= f
       end do   
    end do
    fatal = fatal .or. mismatch

    ! Connectivity is inconsistent; report and punt
    call TLS_fatal_if_any (fatal, 'MESH_DIAGNOSTICS: Inconsistent connectivity')

    ! Test the all-neighbor connectivity
    fatal = .NOT. Test_All_Neighbors()
    call TLS_fatal_if_any (fatal, 'MESH_DIAGNOSTICS: Mesh All Neighbors data structure not correct')

    ! Check mesh quality and ensure no tangled or folded cells
    ! via the Jacobian determinant criterion
    call TWISTED_CELL_TEST ()

  END SUBROUTINE MESH_DIAGNOSTICS

  FUNCTION NODE_CONNECTIVITY () RESULT(NODE_EDGES)
    !=======================================================================
    ! Purpose(s):
    !
    !   Given a general mesh (through use association of Mesh and Vertex)
    !   with cell->vertex (element->node) connectivity, construct
    !   the node<->node (vertex<->vertex) connectivity.
    !   The result is a pointer to a 2d array of node connections.
    !   Each element of the array is an edge in the node graph.  An edge
    !   is a (start_node, end_node) pair.
    !   The array is distributed.
    !
    !=======================================================================
    use mesh_module,          only: Mesh, Cell_Edge
    use mesh_parameter_module, only: nec, nnodes
    use pgslib_module,        only: PGSLIB_PERMUTE, PGSLIB_GRADE_UP,             &
                                    PGSLIB_GLOBAL_EOSHIFT, PGSLIB_PARITY_PREFIX, &
                                    PGSLIB_SUM_PREFIX, PGSLib_Scatter_SUM
    use var_vector_module

    ! Arguments
    type(int_var_vector), dimension(nnodes) :: Node_Edges

    ! Local Variables
    integer, dimension(:,:), pointer :: Edges
    integer, dimension(:),   pointer :: ITemp, Edge_Rank, Edge_Pack_Index
    logical, dimension(:),   pointer :: Edge_Mask, Edge_Segment
    integer, dimension(:,:), pointer :: Edges_Packed
    integer, dimension(nnodes) :: Edges_Per_Node
    integer :: edge, offset, n_edges, max_n_edges
    integer, dimension(:), allocatable :: Temp_Edge_Src


    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! The partitioner wants node<->node connectivity.  We don't have that, so
    ! we have to construct it from the cell->node connectivity, which we do have.
    ! Each cell is surrounded by nvc nodes.  Some of the nodes surrounding a cell 
    ! are (implicitly) connected.  Pairs of nodes which are adjacent (connected)
    ! are given by the array Cell_Edge([1,2],edge), where edge in [1..nec] (nec is
    ! the number of edges surrounding a cell, 4 in 2 dimensions, 12 
    ! in 3 dimensions), and Cell_Edge is initialized in the routine Assign_Cell_Edge
    ! (in mesh_module).  For a given edge the vertex Cell_Edge(1,edge) is 
    ! connected to the vertex Cell_Edge(2,edge).  The first step of this routine
    ! constructs all the possible edges (pairs of connected vertices).  There
    ! are nec*ncells possible edges.

    ! This over counts the edges, however, since a typical edge is shared by 4
    ! cells.  So, once we've constructed all possible edges we construct a mask
    ! which deselects all redundant edges.  Read the comments below for that algorithm
    ! (uses a sort algorithm).

    ! Partioners such as Chaco want a symmetrized graph of unique edges.  
    ! We pack down to unique
    ! edges, and symmetrize the graph, in the same step.  The packing is done
    ! by enumerating all the unique edges (using Edge_Mask), then permuting
    ! Edges into Edges_Packed, using Edge_Mask and the enumerated index.


    ! ALLOCATE the space we will need.
    max_n_edges = nec*SIZE(Mesh,1)
    ALLOCATE (Edges(2, max_n_edges), Edge_Rank(max_n_edges), ITemp(max_n_edges), &
              Edge_Mask(max_n_edges), Edge_Segment(max_n_edges), Edge_Pack_Index(max_n_edges))
    ! Initialize Edges with all possible edges (nec per cell)
    offset = 1
    do edge= 1,nec
       Edges(1,offset:offset+SIZE(Mesh,1)-1) = Mesh%Ngbr_Vrtx(Cell_Edge(1,edge))
       Edges(2,offset:offset+SIZE(Mesh,1)-1) = Mesh%Ngbr_Vrtx(Cell_Edge(2,edge))
       offset = edge*SIZE(Mesh,1) + 1
    end do

    ! We need to identify identical edges.  For that, we specify that
    ! any edge (v1,v2) must be sepcified with v1<= v2.
    ITemp      = MERGE (Edges(1,:), Edges(2,:), Edges(1,:) <= Edges(2,:))
    Edges(2,:) = MERGE (Edges(1,:), Edges(2,:), Edges(1,:) > Edges(2,:))
    Edges(1,:) = ITemp
    
    ! Sort to get all edges with the same v1 adjacent in the list.
    Edge_Rank = PGSLib_GRADE_UP (Edges(1,:))
    call PGSLib_PERMUTE (Edges(1,:), Edges(1,:), Edge_Rank)
    call PGSLib_PERMUTE (Edges(2,:), Edges(2,:), Edge_Rank)

    ! All edges in the same segment have the same v1.
    Edge_Segment = Edges(1,:) /= PGSLib_GLOBAL_EOSHIFT (Edges(1,:), SHIFT=-1, BOUNDARY=-1)
    Edge_Segment = PGSLib_Parity_PREFIX (Edge_Segment)

    ! Sort the v2 so that identical edges will be adjacent.
    Edge_Rank = PGSLib_GRADE_UP (Edges(2,:), SEGMENT = Edge_Segment)
    call PGSLib_PERMUTE (Edges(2,:), Edges(2,:), Edge_Rank)

    ! Mask unique edges.  All redundant edges are masked out.
    Edge_Mask = ((Edges(1,:) /= PGSLib_GLOBAL_EOSHIFT(Edges(1,:), SHIFT=-1, BOUNDARY=-1)) .or.  &
                 (Edges(2,:) /= PGSLib_GLOBAL_EOSHIFT(Edges(2,:), SHIFT=-1, BOUNDARY=-1))) .and. &
                 (Edges(1,:) /= Edges(2,:))

!!$ Replaced with chunk of code below July 9/02, RCF
!!$    ! n_edges_tot is the total number of unique edges
!!$    ! n_edges is how many edges on this processor (will vary from PE to PE).
!!$    n_edges_tot = 2*PGSLib_GLOBAL_COUNT (Edge_Mask)
!!$    n_edges     =   PGSLib_BLOCK_SIZE (n_edges_tot)
!!$ Remove when code below is working.

    ! The Edges_Packed(1, edge) is the tail node of "edge".  We want Edges_Packed
    ! distributed so that each edge is owned by the processor which owns the tail node.
    ! All we need to do is count the number of edges that should be on each processor,

    ! Since we need a symmetrized graph, edge (a,b) will also generate edge (b,a).
    ! So, in the counting we need to note that for each edge (a,b), that will
    ! generate an edge at both nodes a & b.

    Edges_Per_Node = 0
    ALLOCATE(Temp_Edge_Src(SIZE(Edges,2)))
    Temp_Edge_Src = (/ (1, edge = 1, SIZE(Edges,2))/)

    call PGSLib_Scatter_Sum(DEST   = Edges_Per_Node,                 &
                            SOURCE = Temp_Edge_Src,                  &
                            INDEX  = Edges(1,:),                     &
                            MASK   = Edge_Mask)

    call PGSLib_Scatter_Sum(DEST   = Edges_Per_Node,                 &
                            SOURCE = Temp_Edge_Src,                  &
                            INDEX  = Edges(2,:),                     &
                            MASK   = Edge_Mask)

    DEALLOCATE(Temp_Edge_Src)

    ! Total number of edges on each processor is the sum of all the edges on all the nodes
    ! on that processor.  That is, a local sum over edges_per_node.
    n_edges = SUM(Edges_Per_Node)
    
    ! Chaco needs a symmetrized graph.  We provide that by
    ! doubling the edges we just found.
    ALLOCATE (Edges_Packed(2,n_edges))
    Edges_Packed = 0

    ! Construct the index for packing by enumerate the unique edges.
    ! Use Edge_Rank as temporary.
    Edge_Rank = 2
    Edge_Pack_Index = PGSLib_SUM_PREFIX (Edge_Rank, Mask=Edge_Mask)
    
    ! This is the even sites
    ! Index assigns to 2:n_edges:2 in destination

    call PGSLib_PERMUTE (Edges_Packed(1,:), Edges(1,:), Edge_Pack_Index, Edge_Mask)
    call PGSLib_PERMUTE (Edges_Packed(2,:), Edges(2,:), Edge_Pack_Index, Edge_Mask)

    ! This is the odd sites
    ! Index assigns to 1:n_edges:2 in destination
    Edge_Pack_Index = Edge_Pack_Index - 1

    call PGSLib_PERMUTE (Edges_Packed(2,:), Edges(1,:), Edge_Pack_Index, Edge_Mask)
    call PGSLib_PERMUTE (Edges_Packed(1,:), Edges(2,:), Edge_Pack_Index, Edge_Mask)

    ! Now Edges_Packed has a symmetrized graph of unique edges, and we
    ! can get rid of a bunch of unneeded temporaries.
    DEALLOCATE (Edges, Edge_Rank, ITemp, Edge_Mask, Edge_Segment, Edge_Pack_Index)

    ! We have to get the edges in order again.
    ! Use same algorithm as above, except
    ! this time don't want to enforce v1 <= v2, since we want a symmetric graph.
    ALLOCATE (Edge_Rank(n_edges), Edge_Segment(n_edges))

    ! Sort the first vertex of each edge.
    Edge_Rank = 0
    Edge_Rank = PGSLib_GRADE_UP (Edges_Packed(1,:))
    call PGSLib_PERMUTE (Edges_Packed(1,:), Edges_Packed(1,:), Edge_Rank)
    call PGSLib_PERMUTE (Edges_Packed(2,:), Edges_Packed(2,:), Edge_Rank)

    ! Segments have the same first vertex.
    Edge_Segment = Edges_Packed(1,:) /= PGSLib_GLOBAL_EOSHIFT (Edges_Packed(1,:), SHIFT=-1, BOUNDARY=-1)
    Edge_Segment = PGSLib_PARITY_PREFIX (Edge_Segment)

    Edge_Rank = PGSLib_GRADE_UP (Edges_Packed(2,:), Edge_Segment)

    ! Within each segment, sort the second vertex of each edge.
    call PGSLib_PERMUTE (Edges_Packed(1,:), Edges_Packed(1,:), Edge_Rank)
    call PGSLib_PERMUTE (Edges_Packed(2,:), Edges_Packed(2,:), Edge_Rank)
    DEALLOCATE (Edge_Rank, Edge_Segment)

    ! Want to return a distributed VarVector, since that is the natural format.
    ! Each node gets a list of connected nodes, so need a vector for each node.
    
    call CREATE(Node_Edges, SIZES = Edges_Per_Node)

    ! Stuff data into Connectivity
    ! We only need Edges_Packed(2,:), since Edges_Packed(1,:) is implicit.
    ! That is, each edge is owned by the node that owns the tail of that edge.
    ITemp => Edges_Packed(2,:)
    Node_Edges = ITemp

    ! Check that the tails are owned by the appropriate nodes.
!    call CREATE(NodeTail, SIZES = Edges_Per_Node)

!    ITemp => Edges_Packed(1,:)
!    NodeTail = ITemp
!    do node = 1, SIZE(NodeTail)
!       Tails => FLATTEN(NodeTail(node))
!       if (ANY(Tails /= node)) then
!          error_code = tty_out_err_stream .OUT. &
!                       'ERROR: In NODE_CONNECTIVITY, some edge tails are not on the right processor' .OUT. &
!                       END_LINE
!       end if
!       Nullify(Tails)
!    end do

!    call DESTROY(NodeTail)

    DEALLOCATE(Edges_Packed)

  END FUNCTION NODE_CONNECTIVITY

  SUBROUTINE SET_DEGENERATE_FACES(Mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Find the degenerate faces for any degenerate hex elements and set
    !   the appropriate flages in mesh%ngbr_cell and mesh%ngbr_face.
    !
    !=======================================================================
    use mesh_module, only: MESH_CONNECTIVITY, DEGENERATE_FACE

    type(MESH_CONNECTIVITY), dimension(:), intent(INOUT) :: Mesh

    integer :: i

    do i = 1,SIZE(Mesh,1)
       ! Tets, prisms and pyramids all have a degenerate face 6
       if (Mesh(i)%Ngbr_Vrtx(5) == Mesh(i)%Ngbr_Vrtx(8)) then
          Mesh(i)%Ngbr_Face(6) = DEGENERATE_FACE
          Mesh(i)%Ngbr_Cell(6) = DEGENERATE_FACE
       end if

       ! Tets also have a degenerate face 2
       if(Mesh(i)%Ngbr_Vrtx(1) == Mesh(i)%Ngbr_Vrtx(2)) then
          Mesh(i)%Ngbr_Face(2) = DEGENERATE_FACE
          Mesh(i)%Ngbr_Cell(2) = DEGENERATE_FACE
       end if
    end do

  END SUBROUTINE SET_DEGENERATE_FACES

END MODULE MESH_UTILITIES
