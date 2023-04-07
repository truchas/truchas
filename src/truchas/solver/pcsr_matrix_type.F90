!!
!! PCSR_MATRIX_TYPE
!!
!! This module defines distributed data structures for sparse matrices stored
!! in compressed sparse row (CSR) format.  The implementation is limited to
!! square matrices with a symmetric non-zero structure, but general
!! non-symmetric values.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, February 2015.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines the derived types PCSR_MATRIX and PCSR_GRAPH that
!! together describe a distributed sparse matrix using the CSR format.  A
!! PCSR_GRAPH object describes the symmetric non-zero structure of a sparse
!! matrix.  A PCSR_MATRIX object composes a PCSR_GRAPH object with the
!! corresponding non-zero matrix values to complete the description of a
!! sparse matrix.  Multiple PCSR_MATRIX objects can share the same PCSR_GRAPH
!! object.
!!
!! The derived type PCSR_GRAPH has the following type bound procedures.  The
!! process for defining the value of the object begins with a call to INIT to
!! specify the size of the matrix, followed by one or more calls to ADD_EDGE
!! and/or ADD_CLIQUE to specify the non-zero elements of the matrix, concluding
!! with a call to ADD_COMPLETE.  The structure of an N-by-N sparse matrix A is
!! described by an N-node graph: the graph contains the edge (i,j) if and only
!! if A_ij is a non-zero element (potentially).  This implementation assumes a
!! symmetric non-zero structure, and thus the constructed graph is undirected,
!! containing the edge (j,i) whenever it contains the edge (i,j).
!!
!!  INIT (ROW_IP) initializes the graph object; ROW_IP is a type IP_DESC object
!!    that describes how the row indices are distributed across processes.  The
!!    column indices are assumed to have the same distribution.  This must be
!!    the first call.
!!
!!  ADD_EDGE (ROW, COL) adds the edge from node ROW to node COL to the sparse
!!    matrix graph.  This marks element (ROW, COL) as a non-zero element of
!!    the matrix.  ROW is a scalar, and COL a scalar or rank-1 array.  The
!!    edge(s) from COL to ROW is also added automatically.  Adding an edge
!!    already in the graph has no effect and is allowed.
!!
!!  ADD_CLIQUE (INDICES) adds the edges from i to j, where i and j are values
!!    in the rank-1 array INDICES, to the sparse matrix graph.  Adding an edge
!!    already in the graph has no effect and is allowed.
!!
!!  ADD_COMPLETE () performs the final initialization of the object.  After
!!    this call the graph object is ready for use.  Note that no additonal
!!    edges can be added at this point.
!!
!!  DEFINED() returns true if ADD_COMPLETE() has been called and the object
!!    is thus fully initialized.
!!
!! The derived type PCSR_MATRIX is initialized using one of the following type
!! bound subroutines.  This allocates storage for the non-zero matrix values,
!! but the actual values are not yet defined.
!!
!!  INIT (GRAPH [,TAKE_GRAPH]) initializes the sparse matrix object.  GRAPH is
!!    a pointer to a PCSR_GRAPH object; GRAPH%DEFINED() must return true.  If
!!    the optional logical argument TAKE_GRAPH is present and true, then the
!!    matrix object assumes ownership of the graph, and the graph will be
!!    finalized when the matrix is finalized.  Otherwise the matrix object
!!    considers the graph object, which it maintains a reference to, to be
!!    a shared resource, and will not finalize it with the matrix is finalized.
!!    Storage for the non-zero matrix values is allocated but not initialized.
!!
!!  INIT (MOLD) initializes the sparse matrix object using the same graph used
!!    by the PCSR_MATRIX object MOLD.  The matrix object regards the graph as a
!!    shared resource.  Storage for the non-zero matrix values is allocated but
!!    not initialized.
!!
!! The sparse matrix values can be accessed (and modified) directly via the
!! data components of the object: if ARRAY is a PCSR_MATRIX object then
!!    ARRAY%VALUES(ARRAY%GRAPH%XADJ(J):ARRAY%GRAPH%XADJ(J+1)-1)
!! are the values of the the non-zero elements in row J with corresponding
!! column indices
!!    ARRAY%GRAPH%ADJNCY(ARRAY%GRAPH%XADJ(J):ARRAY%GRAPH%XADJ(J+1)-1)
!! Such direct access can usually be avoided by using one of the following
!! type bound procedures.
!!
!!  SET_ALL(VALUE) sets all matrix values to VALUE.
!!
!!  SET(ROW, COL, VALUE) sets the value of the specified matrix element to
!!    VALUE.  It is an error if the matrix element does not belong to the
!!    non-zero structure.
!!
!!  INCREMENT(ROW, COL, VALUE) adds VALUE to the existing value of the specified
!!    matrix element.  It is an error if the matrix element does not belong to
!!    the non-zero structure.
!!
!!  PROJECT_OUT(INDEX) sets the value of all matrix elements in row and column
!!    INDEX to zero.
!!
!!  GET_ROW_VIEW(ROW, VALUES, INDICES) returns a view of the values in the
!!    specified row of the sparse matrix.  The target of the returned real
!!    rank-1 array pointer VALUES contains the values, and the target of the
!!    returned integer rank-1 array pointer INDICES contains the corresponding
!!    column indices.  The object must have the target attribute.  Use with
!!    caution.  The pointers cannot be deallocated. The matrix values may be
!!    modified, but the column indices should not be.
!!
!!  GET_DIAG_COPY(DIAG) puts a copy of the diagonal elements of the sparse
!!    matrix into the rank-1 real array DIAG.
!!
!!  GRAPH() returns a TYPE(PCRS_GRAPH) pointer to the sparse matrix graph.
!!

#include "f90_assert.fpp"

module pcsr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use index_map_type
  use graph_type
  implicit none
  private

  type, public :: pcsr_graph
    type(index_map), pointer :: row_imap => null()  ! reference only -- do not own
    integer, allocatable :: xadj(:), adjncy(:)  ! static graph
    type(graph), allocatable, private :: g      ! dynamic graph -- temporary
  contains
    procedure :: init => pcsr_graph_init
    procedure, private :: pcsr_graph_add_edge_one
    procedure, private :: pcsr_graph_add_edge_many
    generic   :: add_edge => pcsr_graph_add_edge_one, pcsr_graph_add_edge_many
    procedure :: add_clique => pcsr_graph_add_clique
    procedure :: add_complete => pcsr_graph_add_complete
    procedure :: defined => pcsr_graph_defined
  end type pcsr_graph

  type, public :: pcsr_matrix
    real(r8), allocatable :: values(:), x_(:)
    type(pcsr_graph), pointer :: graph => null()
    logical, private :: graph_dealloc = .false.
    integer :: nrow = 0, nrow_onP = 0
    real(r8), allocatable :: diag(:)
  contains
    procedure, private :: pcsr_matrix_init
    procedure, private :: pcsr_matrix_init_mold
    generic   :: init => pcsr_matrix_init, pcsr_matrix_init_mold
    procedure :: graph_ptr => pcsr_matrix_graph
    procedure :: set_all => pcsr_matrix_set_all
    procedure :: set => pcsr_matrix_set
    procedure :: add_to => pcsr_matrix_add_to
    procedure :: project_out => pcsr_matrix_project_out
    procedure :: get_row_view => pcsr_matrix_get_row_view
    procedure :: get_diag_copy => pcsr_matrix_get_diag_copy
    procedure :: copy_to_ijmatrix
    procedure :: matvec
    procedure :: update_diag
    final :: pcsr_matrix_delete
  end type pcsr_matrix

contains

  subroutine pcsr_graph_init (this, row_imap)
    class(pcsr_graph), intent(out) :: this
    type(index_map), intent(in), target :: row_imap
    this%row_imap => row_imap
    allocate(this%g)
    call this%g%init (row_imap%local_size, self_edge=.true.)
  end subroutine pcsr_graph_init

  subroutine pcsr_graph_add_clique (this, indices)
    class(pcsr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:)
    ASSERT(allocated(this%g))
    call this%g%add_clique (indices)
  end subroutine pcsr_graph_add_clique

  subroutine pcsr_graph_add_edge_one (this, row, col)
    class(pcsr_graph), intent(inout) :: this
    integer, intent(in) :: row, col
    ASSERT(allocated(this%g))
    call this%g%add_edge (row, col)
  end subroutine pcsr_graph_add_edge_one

  subroutine pcsr_graph_add_edge_many (this, row, col)
    class(pcsr_graph), intent(inout) :: this
    integer, intent(in) :: row, col(:)
    ASSERT(allocated(this%g))
    call this%g%add_edge (row, col)
  end subroutine pcsr_graph_add_edge_many

  subroutine pcsr_graph_add_complete (this)
    class(pcsr_graph), intent(inout) :: this
    ASSERT(allocated(this%g))
    call this%g%get_adjacency (this%xadj, this%adjncy)
    deallocate(this%g)
  end subroutine pcsr_graph_add_complete

  logical function pcsr_graph_defined (this)
    class(pcsr_graph), intent(in) :: this
    pcsr_graph_defined = allocated(this%xadj)
  end function pcsr_graph_defined

  !! Final subroutine for PCSR_MATRIX type objects.
  subroutine pcsr_matrix_delete (this)
    type(pcsr_matrix), intent(inout) :: this
    if (this%graph_dealloc) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine pcsr_matrix_delete

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine pcsr_matrix_init (this, graph, take_graph)
    class(pcsr_matrix), intent(out) :: this
    type(pcsr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%graph => graph
    if (present(take_graph)) this%graph_dealloc = take_graph
    allocate(this%values(size(graph%adjncy)))
    this%nrow = graph%row_imap%local_size
    this%nrow_onP = graph%row_imap%onP_size
  end subroutine pcsr_matrix_init

  !! Initialize the sparse matrix using the non-zero structure from MOLD.
  subroutine pcsr_matrix_init_mold (this, mold)
    class(pcsr_matrix), intent(out) :: this
    class(pcsr_matrix), intent(in)  :: mold
    call pcsr_matrix_init (this, mold%graph, take_graph=.false.)
  end subroutine pcsr_matrix_init_mold

  !! Return a pointer to the sparse matrix graph.
  function pcsr_matrix_graph (this) result (graph)
    class(pcsr_matrix), intent(in) :: this
    type(pcsr_graph), pointer :: graph
    graph => this%graph
  end function pcsr_matrix_graph

  !! Set all matrix elements to the given value.
  subroutine pcsr_matrix_set_all (this, value)
    class(pcsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    ASSERT(allocated(this%values))
    this%values = value
  end subroutine pcsr_matrix_set_all

  !! Set the specified matrix element to the given value.
  subroutine pcsr_matrix_set (this, row, col, value)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    n = array_index(this%graph, row, col)
    ASSERT(n /= 0)
    this%values(n) = value
  end subroutine pcsr_matrix_set

  !! Increment the specified matrix element by the given value.
  subroutine pcsr_matrix_add_to (this, row, col, value)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    n = array_index(this%graph, row, col)
    ASSERT(n /= 0)
    this%values(n) = this%values(n) + value
  end subroutine pcsr_matrix_add_to

  !! Zero-out the values of the specified row and column elements.
  subroutine pcsr_matrix_project_out (this, index)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    m = index
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%values(lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      do lnm = this%graph%xadj(n), this%graph%xadj(n+1)-1
        if (this%graph%adjncy(lnm) == m) exit
      end do
      ASSERT(lnm < this%graph%xadj(n+1))
      this%values(lnm) = 0.0_r8
    end do
  end subroutine pcsr_matrix_project_out

  !! Return pointers to the values and corresponding column indices
  !! of the specified sparse matrix row.
  subroutine pcsr_matrix_get_row_view (this, row, values, indices)
    class(pcsr_matrix), intent(in), target :: this
    integer, intent(in) :: row
    real(r8), pointer :: values(:)
    integer, pointer :: indices(:)
    ASSERT(row >= 1 .and. row <= this%nrow)
    values => this%values(this%graph%xadj(row):this%graph%xadj(row+1)-1)
    indices => this%graph%adjncy(this%graph%xadj(row):this%graph%xadj(row+1)-1)
  end subroutine pcsr_matrix_get_row_view

  !! Put a copy of the matrix diagonal into the given rank-1 array D.
  !! NB: Elements not written retain their original value.  This means that D
  !!     retains its input value for rows not containing a diagonal element.
  !!     It would more accurate to set those elements to 0, but I'm not sure
  !!     of the current use case.
  subroutine pcsr_matrix_get_diag_copy (this, diag)
    class(pcsr_matrix), intent(in) :: this
    real(r8), intent(inout) :: diag(:)
    integer :: i, k
    ASSERT(size(diag) >= this%nrow_onP)
    do i = 1, this%nrow_onP
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        if (this%graph%adjncy(k) == i) then
          diag(i) = this%values(k)
          exit
        end if
      end do
    end do
  end subroutine pcsr_matrix_get_diag_copy

  !! Auxillary function that computes the index within the compressed
  !! value array corresponding to the specified matrix element.
  pure integer function array_index (graph, row, col) result (k)
    type(pcsr_graph), intent(in) :: graph
    integer, intent(in) :: row, col
    do k = graph%xadj(row), graph%xadj(row+1)-1
      if (graph%adjncy(k) == col) return
    end do
    k = 0
  end function array_index

  !! This auxillary routine copies a PCSR_MATRIX object SRC to an equivalent
  !! HYPRE_IJMatrix object.  The HYPRE matrix is created if it does not exist.
  !! Otherwise the elements of the existing HYPRE matrix are overwritten with
  !! the values from SRC.  In the latter case the sparsity pattern of the two
  !! matrices must be identical.
  !!
  !! The parallel CSR matrix is defined over local indices, both on-process
  !! and off-process.  By nature of our particular construction, the on-process
  !! rows of this matrix are complete and describe a partitioning of the global
  !! matrix by rows.  The off-process rows, however, are partial and extraneous
  !! and should be ignored.
  subroutine copy_to_ijmatrix(src, matrix)

    use fhypre

    class(pcsr_matrix), intent(in) :: src
    type(hypre_obj), intent(inout) :: matrix

    integer :: j, ierr, ilower, iupper, nrows, nnz
    integer, allocatable :: ncols_onP(:), ncols_offP(:), ncols(:), rows(:), cols(:)

    nrows  = src%graph%row_imap%onp_size
    ilower = src%graph%row_imap%first_gid
    iupper = src%graph%row_imap%last_gid

    call fHYPRE_ClearAllErrors

    if (.not.hypre_associated(matrix)) then
      call fHYPRE_IJMatrixCreate(ilower, iupper, ilower, iupper, matrix, ierr)
      !! For each row we know how many column entries are on-process and how many
      !! are off-process.  HYPRE is allegedly much faster at forming its CSR matrix
      !! if it knows this info up front.
      allocate(ncols_onP(nrows), ncols_offP(nrows))
      do j = 1, nrows
        ncols_offP(j) = count(src%graph%adjncy(src%graph%xadj(j):src%graph%xadj(j+1)-1) > nrows)
        ncols_onP(j)  = src%graph%xadj(j+1) - src%graph%xadj(j) - ncols_offP(j)
      end do
      call fHYPRE_IJMatrixSetDiagOffdSizes(matrix, ncols_onP, ncols_offP, ierr)
      deallocate(ncols_onP, ncols_offP)
      !! Let HYPRE know that we won't be setting any off-process matrix values.
      call fHYPRE_IJMatrixSetMaxOffProcElmts(matrix, 0, ierr)
      INSIST(ierr == 0)
    end if

    !! After initialization the HYPRE matrix elements can be set.
    call fHYPRE_IJMatrixInitialize(matrix, ierr)
    INSIST(ierr == 0)

    !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    !! nonzero structure of the matrix and the values of those elements. HYPRE
    !! expects global row and column indices.
    nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    allocate(ncols(nrows), rows(nrows), cols(nnz))
    rows = [ (j, j = ilower, iupper) ]
    ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    cols = src%graph%row_imap%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    call fHYPRE_IJMatrixSetValues(matrix, nrows, ncols, rows, cols, src%values, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble(matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

  subroutine matvec(this, x, b)

    class(pcsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: x(:)
    real(r8), intent(out) :: b(:)

    integer :: i, j, xj
    real(r8), pointer :: values(:) => null()
    integer, pointer :: indices(:) => null()

    ASSERT(this%nrow_onP <= size(x))

    if (.not.allocated(this%x_)) allocate(this%x_(this%nrow))
    this%x_(:this%nrow_onP) = x(:this%nrow_onP)
    call this%graph%row_imap%gather_offp(this%x_)

    do i = 1, this%nrow_onP
      call this%get_row_view(i, values, indices)
      b(i) = 0
      do xj = 1, size(indices)
        j = indices(xj)
        b(i) = b(i) + values(xj) * this%x_(j)
      end do
    end do

  end subroutine matvec


  subroutine update_diag(this)

    class(pcsr_matrix), intent(inout) :: this

    integer :: i, k
    integer, pointer :: indices(:) => null()
    real(r8), pointer :: values(:) => null()

    if (.not.allocated(this%diag)) allocate(this%diag(this%nrow))

    do i = 1, this%nrow
      this%diag(i) = 0
      call this%get_row_view(i, values, indices)
      do k = 1, size(indices)
        if (indices(k) == i) this%diag(i) = values(k)
      end do
    end do

  end subroutine update_diag

end module pcsr_matrix_type
