!!
!! GRAPH_TYPE
!!
!! This module defines a data structure for representing an N-node graph.  The
!! implemented functionality is limited to dynamically building a graph and then
!! extracting a compact static description of the graph in a standard format.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2014.  Adaptation of GraphModule (1998) for Fortran 2008.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines the GRAPH derived type. Objects of this type are properly
!! finalized when deallocated or otherwise cease to exist.  Objects of this type
!! should not be used in assignment statements, as only intrinsic assignment is
!! defined. The type has the following type bound procedures.
!!
!!  INIT (N [,DIRECTED] [,SELF_EDGE]) initializes the object to describe an
!!    N-node graph.  This must be the first call.  The graph starts with no
!!    edges defined.  The remaining optional logical arguments specify the
!!    type of graph.  If SELF_EDGE is present with value true, edges linking
!!    a node to itself are allowed; otherwise such self edges are not stored
!!    (but may be "added" -- see ADD_EDGE).  If DIRECTED is present with value
!!    true, the graph is directed; otherwise it is undirected, and whenever an
!!    edge from one node to another is added, the oppositely directed edge is
!!    automatically added.
!!
!!  INIT (N, M) initializes the object to describe a bipartite graph where
!!    N is the number of nodes in the domain set and M the number of nodes
!!    in the range set, with directed edges from domain nodes to range nodes.
!!    Only the methods ADD_EDGE and GET_ADJACENCY are applicable to a bipartite
!!    graph.
!!
!!  ADD_EDGE (FROM, TO) adds the edge from node FROM to node TO. TO may be a
!!    scalar or rank-1 array, and in the latter case a sequence of edges is
!!    added. When the graph is not a bipartite graph, edges that are not
!!    consistent with the type of graph (SELF_EDGE) are silently ignored, and
!!    for undirected graph, the edge or edges from TO to FROM are also added.
!!
!!  ADD_CLIQUE (CLIQUE) adds a clique of edges.  CLIQUE is a rank-1 array, and
!!    an edge from CLIQUE(j) to CLIQUE(k) is added for all pairs (j,k).  Edges
!!    that are not consistent with the type of graph (SELF_EDGE) are silently
!!    ignored.
!!
!!  GET_ADJACENCY (XADJ, ADJNCY) returns the adjacency structure of the graph
!!    in a commonly-used commpressed format.  XADJ and ADJNCY are intent(out)
!!    rank-1 allocatable integer arrays, which are allocated and initialized
!!    by the procedure.  ADJNCY(XADJ(j):XADJ(j+1)-1) is the list of neighbor
!!    nodes of node j.  The nodes in each of the neighbor lists are sorted in
!!    ascending order.  The returned XADJ will have size N+1, and ADJNCY size
!!    XADJ(N+1)-1.
!!
!!  VERTEX_COLORING (VCOLOR, NCOLOR) returns a vertex coloring of the graph.
!!    A vertex coloring is a mapping from vertices to a finite set (the colors)
!!    such that adjacent vertices are assigned different colors.  NCOLOR is an
!!    intent(out) integer specifying the total number of colors used.  VCOLOR is
!!    an intent(out) rank-1 array of integers in the range [1, NCOLOR], where
!!    VCOLOR(j) is the color of vertex j.  This algorithm guarantees that NCOLOR
!!    <= D+1, where D is the maximum degree of any vertex.  Vertex coloring is
!!    undefined for SELF_EDGE graphs, and such graphs are silently ignored.
!!

#include "f90_assert.fpp"

module graph_type

  use integer_set_type
  implicit none
  private

  type, public :: graph
    private
    integer :: n = 0, m = 0
    logical :: directed = .false.
    logical :: self_edge = .false.
    logical :: bigraph = .false.
    type(integer_set), allocatable :: nbrs(:)
  contains
    procedure, private :: graph_init
    procedure, private :: bigraph_init
    generic :: init => graph_init, bigraph_init
    procedure, private :: graph_add_edge_one
    procedure, private :: graph_add_edge_many
    generic :: add_edge => graph_add_edge_one, graph_add_edge_many
    procedure, private :: graph_add_clique_one, graph_add_clique_many
    generic :: add_clique => graph_add_clique_one, graph_add_clique_many
    procedure :: get_adjacency => graph_get_adjacency
    procedure, private :: graph_get_components
    procedure, private :: graph_get_components_mask
    generic :: get_components => graph_get_components, graph_get_components_mask
    procedure :: vertex_coloring => graph_vertex_coloring
  end type graph

contains

  !! Initialize the graph object.
  subroutine graph_init (this, n, directed, self_edge)
    class(graph), intent(out) :: this
    integer, intent(in) :: n
    logical, intent(in), optional :: directed, self_edge
    ASSERT(n >= 0)
    this%n = n
    this%m = n
    allocate(this%nbrs(n))
    if (present(directed))  this%directed  = directed
    if (present(self_edge)) this%self_edge = self_edge
  end subroutine graph_init

  !! Initialize the bipartite graph object
  subroutine bigraph_init (this, n, m)
    class(graph), intent(out) :: this
    integer, intent(in) :: n, m
    ASSERT(n >= 0)
    ASSERT(m >= 0)
    this%n = n
    this%m = m
    allocate(this%nbrs(n))
    this%bigraph   = .true.
    this%directed  = .true.
    this%self_edge = .true.
  end subroutine bigraph_init

  !! Add the specified edge to the graph.
  subroutine graph_add_edge_one (this, from, to)
    class(graph), intent(inout) :: this
    integer, intent(in) :: from, to
    ASSERT(from >= 1 .and. from <= this%n)
    ASSERT(to >= 1 .and. to <= this%m)
    if (.not.this%self_edge .and. from == to) return
    call this%nbrs(from)%add (to)
    if (.not.this%directed) call this%nbrs(to)%add (from)
  end subroutine graph_add_edge_one

  !! Add the specified set of edges to the graph.
  subroutine graph_add_edge_many (this, from, to)
    class(graph), intent(inout) :: this
    integer, intent(in) :: from, to(:)
    integer :: j
    ASSERT(from >= 1 .and. from <= this%n)
    do j = 1, size(to)
      ASSERT(to(j) >= 1 .and. to(j) <= this%m)
      if (.not.this%self_edge .and. from == to(j)) cycle
      call this%nbrs(from)%add (to(j))
      if (.not.this%directed) call this%nbrs(to(j))%add (from)
    end do
  end subroutine graph_add_edge_many

  !! Add the specified clique of edges to the graph.
  subroutine graph_add_clique_one(this, clique)
    class(graph), intent(inout) :: this
    integer, intent(in) :: clique(:)
    integer :: j, k, from, to
    INSIST(.not.this%bigraph)
    do j = 1, size(clique)
      from = clique(j)
      ASSERT(from >= 1 .and. from <= this%n)
      do k = 1, size(clique)
        to = clique(k)
        ASSERT(to >= 1 .and. to <= this%n)
        if (to == from .and. .not.this%self_edge) cycle
        call this%nbrs(from)%add (to)
      end do
    end do
  end subroutine

  !! Add the specified set of cliques to the graph.
  subroutine graph_add_clique_many(this, clique)
    class(graph), intent(inout) :: this
    integer, intent(in) :: clique(:,:)
    integer :: j
    do j = 1, size(clique,dim=2)
      call graph_add_clique_one(this, clique(:,j))
    end do
  end subroutine

  !! Return the adjacency structure of the graph.
  subroutine graph_get_adjacency (this, xadj, adjncy)
    class(graph), intent(in) :: this
    integer, allocatable, intent(out) :: xadj(:), adjncy(:)
    integer :: j
    allocate(xadj(this%n+1))
    xadj(1) = 1
    do j = 1, this%n
      xadj(j+1) = xadj(j) + this%nbrs(j)%size()
    end do
    allocate(adjncy(xadj(this%n+1)-1))
    do j = 1, this%n
      call this%nbrs(j)%copy_to_array(adjncy(xadj(j):))
    end do
  end subroutine graph_get_adjacency

  !! Return the connected components of the graph.
  subroutine graph_get_components (this, ncomp, xcomp, comp)
    class(graph), intent(in) :: this
    integer, intent(out) :: ncomp
    integer, allocatable, intent(out) :: xcomp(:), comp(:)
    integer :: j
    integer :: tag(this%n)
    INSIST(.not.this%bigraph)
    if (this%directed) return
    ncomp = 0
    tag = 0
    do j = 1, this%n
      if (tag(j) /= 0) cycle
      ncomp = ncomp + 1
      call tag_component (j)
    end do
    allocate(xcomp(ncomp+1),comp(this%n))
    !! Prepare XCOMP; nodes in component N will be COMP(XCOMP(N):XCOMP(N+1)-1)
    xcomp = 0
    do j = 1, size(tag)
      xcomp(1+tag(j)) = xcomp(1+tag(j)) + 1
    end do
    xcomp(1) = 1
    do j = 2, size(xcomp)
      xcomp(j) = xcomp(j) + xcomp(j-1)
    end do
    !! Fill COMP; XCOMP(N) stores the next free location for component N.
    do j = 1, size(tag)
      comp(xcomp(tag(j))) = j
      xcomp(tag(j)) = xcomp(tag(j)) + 1
    end do
    !! Restore XCOMP; XCOMP(N) is now the start of component N+1 instead of N
    do j = size(xcomp), 2, -1
      xcomp(j) = xcomp(j-1)
    end do
    xcomp(1) = 1
  contains
    !! Tag all the nodes connected to ROOT with the current component number
    recursive subroutine tag_component (root)
      integer, intent(in) :: root
      type(integer_set_iterator) :: iter
      tag(root) = ncomp
      iter = this%nbrs(root)%begin()
      do while (.not.iter%at_end())
        if (tag(iter%value()) == 0) call tag_component (iter%value())
        call iter%next
      end do
    end subroutine
  end subroutine graph_get_components

  !! Return the connected components of the subgraph defined by mask.
  subroutine graph_get_components_mask (this, mask, ncomp, xcomp, comp)
    class(graph), intent(in) :: this
    logical, intent(in)  :: mask(:)
    integer, intent(out) :: ncomp
    integer, allocatable, intent(out) :: xcomp(:), comp(:)
    integer :: j
    integer :: tag(this%n)
    INSIST(.not.this%bigraph)
    ASSERT(size(mask) == this%n)
    if (this%directed) return
    ncomp = 0
    tag = 0
    do j = 1, this%n
      if (mask(j) .and. tag(j) == 0) then
        ncomp = ncomp + 1
        call tag_component (j)
      end if
    end do
    ASSERT(all((tag > 0) .eqv. mask))
    allocate(xcomp(ncomp+1),comp(count(mask)))
    !! Prepare XCOMP; nodes in component N will be COMP(XCOMP(N):XCOMP(N+1)-1)
    xcomp = 0
    do j = 1, size(tag)
      if (mask(j)) xcomp(1+tag(j)) = xcomp(1+tag(j)) + 1
    end do
    xcomp(1) = 1
    do j = 2, size(xcomp)
      xcomp(j) = xcomp(j) + xcomp(j-1)
    end do
    !! Fill COMP; XCOMP(N) stores the next free location for component N.
    do j = 1, size(tag)
      if (.not.mask(j)) cycle
      comp(xcomp(tag(j))) = j
      xcomp(tag(j)) = xcomp(tag(j)) + 1
    end do
    !! Restore XCOMP; XCOMP(N) is now the start of component N+1 instead of N
    do j = size(xcomp), 2, -1
      xcomp(j) = xcomp(j-1)
    end do
    xcomp(1) = 1
  contains
    !! Tag all the nodes connected to ROOT with the current component number
    recursive subroutine tag_component (root)
      integer, intent(in) :: root
      type(integer_set_iterator) :: iter
      tag(root) = ncomp
      iter = this%nbrs(root)%begin()
      do while (.not.iter%at_end())
        if (mask(iter%value()) .and. tag(iter%value()) == 0) call tag_component (iter%value())
        call iter%next
      end do
    end subroutine
  end subroutine graph_get_components_mask

  !! Returns a vertex coloring of the graph.
  subroutine graph_vertex_coloring (this, vcolor, ncolor)
    use sort_utilities
    class(graph), intent(in) :: this
    integer, allocatable, intent(out) :: vcolor(:)
    integer, intent(out) :: ncolor
    integer :: degree(this%n), perm(this%n)
    logical, allocatable :: avail(:)  ! whether a color is available to a vertex
    type(integer_set_iterator) :: iter
    integer :: i, v, n, c, max_ncolor
    INSIST(.not.this%bigraph)
    if (this%self_edge) return
    allocate(vcolor(this%n))
    vcolor = -1
    !! Sort vertices by degree
    degree = this%nbrs%size()
    call heap_sort (degree, perm)
    !! At most D + 1 colors needed
    max_ncolor = degree(perm(this%n)) + 1
    allocate(avail(max_ncolor))
    !! Color vertices in descending order of degree
    ncolor = 0
    do i = this%n, 1, -1
      v = perm(i)  ! current vertex
      !! Get availabe colors
      avail = .true.
      iter = this%nbrs(v)%begin()
      do while (.not. iter%at_end())
        n = iter%value()  ! neighbor of v
        if (vcolor(n) /= -1) avail(vcolor(n)) = .false.
        call iter%next()
      end do
      !! Assign color
      ! TODO: balance colors rather than use first available?
      c = findloc(avail, .true., dim=1)
      ASSERT(c /= 0)
      vcolor(v) = c
      ncolor = max(c, ncolor)
    end do
  end subroutine graph_vertex_coloring

end module graph_type
