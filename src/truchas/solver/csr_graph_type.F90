!!
!! CSR_GRAPH_TYPE
!!
!! This module defines the derived types CSR_GRAPH and PCSR_GRAPH which
!! describe the non-zero structure of sparse matrices stored in a compressed
!! sparse row (CSR) format. Either type can be used for distributed sparse
!! matrices, but the PCSR_GRAPH type stores additional data that describes
!! the parallel distribution of the row and column index sets.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module csr_graph_type

  use graph_type
  use index_map_type
  implicit none
  private

  type, private :: base_csr_graph
    integer, allocatable :: xadj(:), adjncy(:)
    type(graph), allocatable, private :: g ! temporary during construction
  contains
    procedure, private :: add_edge_one, add_edge_many
    generic   :: add_edge => add_edge_one, add_edge_many
    procedure, private :: add_clique_one, add_clique_many
    generic :: add_clique => add_clique_one, add_clique_many
    procedure :: add_complete
    procedure :: defined
    procedure :: index => index_linear 
    !procedure :: index => index_binary
  end type

  type, extends(base_csr_graph), public :: csr_graph
    integer :: nrow, ncol
  contains
    generic :: init => csr_graph_init, csr_bigraph_init
    procedure, private :: csr_graph_init, csr_bigraph_init
  end type

  type, extends(base_csr_graph), public :: pcsr_graph
    type(index_map), pointer :: row_imap => null() ! unowned reference
    type(index_map), pointer :: col_imap => null() ! unowned reference
  contains
    generic :: init => pcsr_graph_init, pcsr_bigraph_init
    procedure, private :: pcsr_graph_init, pcsr_bigraph_init
  end type

contains

  subroutine csr_graph_init(this, nrow)
    class(csr_graph), intent(out) :: this
    integer, intent(in) :: nrow
    this%nrow = nrow
    this%ncol = nrow
    allocate(this%g)
    call this%g%init(this%nrow, self_edge=.true.)
  end subroutine

  subroutine csr_bigraph_init(this, nrow, ncol)
    class(csr_graph), intent(out) :: this
    integer, intent(in) :: nrow, ncol
    this%nrow = nrow
    this%ncol = ncol
    allocate(this%g)
    call this%g%init(this%nrow, this%ncol)
  end subroutine

  subroutine pcsr_graph_init(this, row_imap)
    class(pcsr_graph), intent(out) :: this
    type(index_map), intent(in), target :: row_imap
    this%row_imap => row_imap
    this%col_imap => row_imap
    allocate(this%g)
    call this%g%init(row_imap%local_size, self_edge=.true.)
  end subroutine

  subroutine pcsr_bigraph_init(this, row_imap, col_imap)
    class(pcsr_graph), intent(out) :: this
    type(index_map), intent(in), target :: row_imap, col_imap
    this%row_imap => row_imap
    this%col_imap => col_imap
    allocate(this%g)
    call this%g%init(row_imap%local_size, col_imap%local_size)
  end subroutine

  subroutine add_edge_one(this, row, col)
    class(base_csr_graph), intent(inout) :: this
    integer, intent(in) :: row, col
    ASSERT(allocated(this%g))
    call this%g%add_edge(row, col)
  end subroutine

  subroutine add_edge_many(this, row, col)
    class(base_csr_graph), intent(inout) :: this
    integer, intent(in) :: row, col(:)
    ASSERT(allocated(this%g))
    call this%g%add_edge(row, col)
  end subroutine

  subroutine add_clique_one(this, indices)
    class(base_csr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine add_clique_many(this, indices)
    class(base_csr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:,:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine add_complete(this)
    class(base_csr_graph), intent(inout) :: this
    ASSERT(allocated(this%g))
    call this%g%get_adjacency(this%xadj, this%adjncy)
    deallocate(this%g)
  end subroutine

  logical function defined(this)
    class(base_csr_graph), intent(in) :: this
    defined = allocated(this%xadj)
  end function

  !! Return the ADJNCY array index corresponding to matrix element (row, col),
  !! or 0 if that element is not stored by the sparse matrix. This version does
  !! a linear search of the appropriate sorted segment of the array.

  pure integer function index_linear(this, row, col) result(k)
    class(base_csr_graph), intent(in) :: this
    integer, intent(in) :: row, col
    do k = this%xadj(row), this%xadj(row+1)-1
      if (this%adjncy(k) == col) return
    end do
    k = 0
  end function

  !! Return the ADJNCY array index corresponding to matrix element (row, col),
  !! or 0 if that element is not stored by the sparse matrix. This version does
  !! a binary search of the appropriate sorted segment of the array. This does
  !! not appear to be significantly more or less efficient that the linear
  !! search for the lengths of adjacency lists encountered in practice (FEM).

  pure integer function index_binary(this, row, col) result(k)
    class(base_csr_graph), intent(in) :: this
    integer, intent(in) :: row, col
    integer :: k1, k2
    k1 = this%xadj(row)
    k2 = this%xadj(row+1) - 1
    do while (k1 <= k2)
      k = (k1 + k2) / 2
      if (this%adjncy(k) < col) then
        k1 = k + 1
      else if (this%adjncy(k) > col) then
        k2 = k - 1
      else ! this%adjncy(k) == col
        return
      end if
    end do
    k = 0
  end function

end module csr_graph_type
