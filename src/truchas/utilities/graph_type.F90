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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
!!  ADD_EDGE (FROM, TO) adds the edge from node FROM to node TO.  TO may be
!!    a scalar or rank-1 array, and in the latter case a sequence of edges is
!!    added.  Edges that are not consistent with the type of graph (SELF_EDGE)
!!    are silently ignored.  In the case of an undirected graph, the edge or
!!    edges from TO to FROM are also added.
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

#include "f90_assert.fpp"

module graph_type

  use integer_set_type
  implicit none
  private

  type, public :: graph
    private
    integer :: n = 0
    logical :: directed = .false.
    logical :: self_edge = .false.
    type(integer_set), allocatable :: nbrs(:)
  contains
    procedure :: init => graph_init
    procedure, private :: graph_add_edge_one
    procedure, private :: graph_add_edge_many
    generic :: add_edge => graph_add_edge_one, graph_add_edge_many
    procedure :: add_clique => graph_add_clique
    procedure :: get_adjacency => graph_get_adjacency
  end type graph

contains

  !! Initialize the graph object.
  subroutine graph_init (this, n, directed, self_edge)
    class(graph), intent(out) :: this
    integer, intent(in) :: n
    logical, intent(in), optional :: directed, self_edge
    ASSERT(n >= 0)
    this%n = n
    allocate(this%nbrs(n))
    if (present(directed))  this%directed  = directed
    if (present(self_edge)) this%self_edge = self_edge
  end subroutine graph_init

  !! Add the specified edge to the graph.
  subroutine graph_add_edge_one (this, from, to)
    class(graph), intent(inout) :: this
    integer, intent(in) :: from, to
    ASSERT(from >= 1 .and. from <= this%n)
    ASSERT(to >= 1 .and. to <= this%n)
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
      ASSERT(to(j) >= 1 .and. to(j) <= this%n)
      if (.not.this%self_edge .and. from == to(j)) cycle
      call this%nbrs(from)%add (to(j))
      if (.not.this%directed) call this%nbrs(to(j))%add (from)
    end do
  end subroutine graph_add_edge_many

  !! Add the specified clique of edges to the graph.
  subroutine graph_add_clique (this, clique)
    class(graph), intent(inout) :: this
    integer, intent(in) :: clique(:)
    integer :: j, k, from, to
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
  end subroutine graph_add_clique

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

end module graph_type
