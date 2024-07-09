!!
!! PCSR_MATRIX_TYPE
!!
!! This module defines a distributed data structure for sparse matrices stored
!! in compressed sparse row (CSR) format.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, February 2015.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module pcsr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use csr_graph_type, only: pcsr_graph
  implicit none
  private

  public :: pcsr_graph

  type, public :: pcsr_matrix
    integer :: nrow = 0, ncol = 0, nrow_onP = 0, ncol_onP = 0
    real(r8), allocatable :: values(:)!, x_(:)
    type(pcsr_graph), pointer :: graph => null()
    logical, private :: graph_is_owned = .false.
    integer, allocatable :: kdiag(:)
  contains
    generic :: init => init_graph,init_mold
    procedure, private :: init_graph, init_mold
    procedure :: set_all
    procedure :: set
    generic :: add_to => add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure, private :: add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure :: project_out
    procedure :: matvec
    procedure :: kdiag_init
    procedure :: get_row_view !TODO? eliminate
    procedure :: copy_to_ijmatrix
    final :: pcsr_matrix_delete
  end type pcsr_matrix

contains

  !! Final subroutine for PCSR_MATRIX type objects.
  subroutine pcsr_matrix_delete(this)
    type(pcsr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine init_graph(this, graph, take_graph)
    class(pcsr_matrix), intent(out) :: this
    type(pcsr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%nrow = graph%row_imap%local_size
    this%ncol = graph%col_imap%local_size
    this%nrow_onP = graph%row_imap%onP_size
    this%ncol_onP = graph%col_imap%onP_size
    allocate(this%values(size(graph%adjncy)), source=0.0_r8)
    this%graph => graph
    if (present(take_graph)) then
      this%graph_is_owned = take_graph
      if (take_graph) nullify(graph)
    end if
  end subroutine

  !! Initialize the sparse matrix using the non-zero structure from MOLD.
  subroutine init_mold(this, mold)
    class(pcsr_matrix), intent(out) :: this
    class(pcsr_matrix), intent(in)  :: mold
    call init_graph(this, mold%graph, take_graph=.false.)
  end subroutine

  !! Set all matrix elements to the given value.
  subroutine set_all(this, value)
    class(pcsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    ASSERT(allocated(this%values))
    this%values = value
  end subroutine

  !! Set the specified matrix element to the given value.
  subroutine set(this, row, col, value)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    n = this%graph%index(row, col)
    ASSERT(n /= 0)
    this%values(n) = value
  end subroutine

  !! Increment the specified matrix element by the given value.
  subroutine add_to_element (this, row, col, value)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    n = this%graph%index(row, col)
    ASSERT(n /= 0)
    this%values(n) = this%values(n) + value
  end subroutine

  !! Increment the principal submatrix specified by the list of row indices
  !! with the corresponding values in the given matrix.

  subroutine add_to_submatrix(this, rows, matrix)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:)
    integer :: i, j, n
    INSIST(associated(this%graph%row_imap, this%graph%col_imap))
    ASSERT(size(matrix,1) == size(matrix,2))
    ASSERT(size(rows) == size(matrix,1))
    ASSERT(minval(rows) >= 1 .and. maxval(rows) <= this%nrow)
    do i = 1, size(rows)
      do j = 1, size(rows)
        n = this%graph%index(rows(i),rows(j))
        ASSERT(n /= 0)
        this%values(n) = this%values(n) + matrix(i,j)
      end do
    end do
  end subroutine

  !! Increment the principal submatrix specified by the list of row indices
  !! with the corresponding values in the given symmetric matrix in upper
  !! packed storage format.

  subroutine add_to_submatrix_upm(this, rows, matrix)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:)
    integer :: i, j, l, n
    INSIST(associated(this%graph%row_imap, this%graph%col_imap))
    ASSERT(size(matrix) == (size(rows)*(size(rows)+1))/2)
    ASSERT(minval(rows) >= 1 .and. maxval(rows) <= this%nrow)
    l = 0
    do j = 1, size(rows)
      do i = 1, j-1
        l = l + 1
        n = this%graph%index(rows(i),rows(j))
        this%values(n) = this%values(n) + matrix(l)
        n = this%graph%index(rows(j),rows(i))
        this%values(n) = this%values(n) + matrix(l)
      end do
      l = l + 1
      n = this%graph%index(rows(j),rows(j))
      this%values(n) = this%values(n) + matrix(l)
    end do
  end subroutine

  !! Zero-out the values of the specified row and column elements.
  !! NB: Assumes symmetric structure

  subroutine project_out(this, index)
    class(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    INSIST(associated(this%graph%row_imap, this%graph%col_imap))
    m = index
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%values(lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      lnm = this%graph%index(n, m)
      ASSERT(lnm > 0)
      this%values(lnm) = 0.0_r8
    end do
  end subroutine

  subroutine matvec(this, x, b, incr)

    class(pcsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: x(:)
    real(r8), intent(inout) :: b(:)
    logical, intent(in), optional :: incr

    integer :: i, k
    real(r8) :: s
    logical :: incr_

    ASSERT(size(x) == this%ncol)
    ASSERT(size(b) >= this%nrow_onP)

    incr_ = .false.
    if (present(incr)) incr_ = incr

    do i = 1, this%nrow_onP
      s = merge(b(i), 0.0_r8, incr_)
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + this%values(k) * x(this%graph%adjncy(k))
      end do
      b(i) = s
    end do

  end subroutine matvec

  subroutine kdiag_init(this)
    class(pcsr_matrix), intent(inout) :: this
    integer :: i
    if (allocated(this%kdiag)) return
    INSIST(this%nrow == this%ncol)
    allocate(this%kdiag(this%nrow))
    do i = 1, this%nrow
      this%kdiag(i) = this%graph%index(i,i)
    end do
    ASSERT(all(this%kdiag > 0))
  end subroutine

  !! Return pointers to the values and corresponding column indices
  !! of the specified sparse matrix row.
  subroutine get_row_view(this, row, values, indices)
    class(pcsr_matrix), intent(in), target :: this
    integer, intent(in) :: row
    real(r8), pointer :: values(:)
    integer, pointer :: indices(:)
    ASSERT(row >= 1 .and. row <= this%nrow)
    values => this%values(this%graph%xadj(row):this%graph%xadj(row+1)-1)
    indices => this%graph%adjncy(this%graph%xadj(row):this%graph%xadj(row+1)-1)
  end subroutine

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

    integer :: j, ierr, ilower, iupper, nrows, jlower, jupper, nnz, ncols
    integer, allocatable :: ncols_onP(:), ncols_offP(:), counts(:), rows(:), cols(:)

    nrows  = src%graph%row_imap%onp_size
    ilower = src%graph%row_imap%first_gid
    iupper = src%graph%row_imap%last_gid
    ncols  = src%graph%col_imap%onp_size
    jlower = src%graph%col_imap%first_gid
    jupper = src%graph%col_imap%last_gid

    call fHYPRE_ClearAllErrors

    if (.not.hypre_associated(matrix)) then
      call fHYPRE_IJMatrixCreate(ilower, iupper, jlower, jupper, matrix, ierr)
      !! For each row we know how many column entries are on-process and how many
      !! are off-process.  HYPRE is allegedly much faster at forming its CSR matrix
      !! if it knows this info up front.
      allocate(ncols_onP(nrows), ncols_offP(nrows))
      do j = 1, nrows
        ncols_offP(j) = count(src%graph%adjncy(src%graph%xadj(j):src%graph%xadj(j+1)-1) > ncols)
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
    allocate(counts(nrows), rows(nrows), cols(nnz))
    rows = [ (j, j = ilower, iupper) ]
    counts = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    cols = src%graph%col_imap%global_index(src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    call fHYPRE_IJMatrixSetValues(matrix, nrows, counts, rows, cols, src%values, ierr)
    deallocate(counts, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble(matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module pcsr_matrix_type
