!!
!! PBSR_MATRIX_TYPE
!!
!! This module defines a distributed data structure for sparse matrices stored
!! in block compressed sparse row (BSR) format.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

module pbsr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use csr_graph_type, only: pcsr_graph
  implicit none
  private

  public :: pcsr_graph

  type, public :: pbsr_matrix
    integer :: nrow, ncol, bsize, nrow_onP, ncol_onP
    real(r8), allocatable :: values(:,:,:)
    type(pcsr_graph), pointer :: graph => null()
    logical, private :: graph_is_owned = .false.
  contains
    generic :: init => init_graph, init_mold
    procedure, private :: init_graph, init_mold
    procedure :: set_all
    procedure :: set
    generic :: add_to => add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure, private :: add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure :: project_out
    procedure :: matvec
    final :: pbsr_matrix_delete
  end type

contains

  !! Final subroutine for pbsr_MATRIX type objects.
  elemental subroutine pbsr_matrix_delete(this)
    type(pbsr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine init_graph(this, bsize, graph, take_graph)
    class(pbsr_matrix), intent(out) :: this
    integer, intent(in) :: bsize
    type(pcsr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%graph => graph
    if (present(take_graph)) this%graph_is_owned = take_graph
    this%nrow = graph%row_imap%local_size
    this%ncol = graph%col_imap%local_size
    this%bsize = bsize
    this%nrow_onP = graph%row_imap%onP_size
    this%ncol_onP = graph%col_imap%onP_size
    allocate(this%values(bsize,bsize,size(this%graph%adjncy)), source=0.0_r8)
  end subroutine

  !! Initialize the sparse matrix using the non-zero graph structure from MOLD.
  subroutine init_mold(this, mold)
    class(pbsr_matrix), intent(out) :: this
    class(pbsr_matrix), intent(in)  :: mold
    call init_graph(this, mold%bsize, mold%graph, take_graph=.false.)
  end subroutine

  !! Set all matrix elements to the given value.
  subroutine set_all(this, value)
    class(pbsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    this%values = value
  end subroutine

  !! Set the value of the specified matrix element to the given value.
  subroutine set(this, row, col, value)
    class(pbsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value(:,:)
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    ASSERT(size(value,1) == this%bsize)
    ASSERT(size(value,2) == this%bsize)
    n = this%graph%index(row,col)
    ASSERT(n /= 0)
    this%values(:,:,n) = value
  end subroutine

  !! Increment the specified matrix element by the given value.
  subroutine add_to_element(this, row, col, value)
    class(pbsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value(:,:)
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    ASSERT(size(value,1) == this%bsize)
    ASSERT(size(value,2) == this%bsize)
    n = this%graph%index(row,col)
    ASSERT(n /= 0)
    this%values(:,:,n) = this%values(:,:,n) + value
  end subroutine

  !! Increment the principle submatrix specified by the list of row indices
  !! with the corresponding values in the given matrix.
  subroutine add_to_submatrix(this, rows, matrix)
    class(pbsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:,:,:)
    integer :: i, j, n
    INSIST(associated(this%graph%row_imap, this%graph%col_imap))
    ASSERT(size(matrix,1) == this%bsize)
    ASSERT(size(matrix,2) == this%bsize)
    ASSERT(size(matrix,3) == size(matrix,4))
    ASSERT(size(rows) == size(matrix,3))
    ASSERT(minval(rows) >= 1 .and. maxval(rows) <= this%nrow)
    do i = 1, size(rows)
      do j = 1, size(rows)
        n = this%graph%index(rows(i),rows(j))
        ASSERT(n /= 0)
        this%values(:,:,n) = this%values(:,:,n) + matrix(:,:,i,j)
      end do
    end do
  end subroutine

  !! Increment the principle submatrix specified by the list of row indices
  !! with the corresponding values in the given symmetric matrix in upper
  !! packed storage format.

  subroutine add_to_submatrix_upm(this, rows, matrix)
    class(pbsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:,:)
    integer :: i, j, l, n
    INSIST(associated(this%graph%row_imap, this%graph%col_imap))
    ASSERT(size(matrix,1) == this%bsize)
    ASSERT(size(matrix,2) == this%bsize)
    ASSERT(size(matrix,3) == (size(rows)*(size(rows)+1))/2)
    ASSERT(minval(rows) >= 1 .and. maxval(rows) <= this%nrow)
    l = 0
    do j = 1, size(rows)
      do i = 1, j-1
        l = l + 1
        n = this%graph%index(rows(i),rows(j))
        this%values(:,:,n) = this%values(:,:,n) + matrix(:,:,l)
        n = this%graph%index(rows(j),rows(i))
        this%values(:,:,n) = this%values(:,:,n) + matrix(:,:,l)
      end do
      l = l + 1
      n = this%graph%index(rows(j),rows(j))
      this%values(:,:,n) = this%values(:,:,n) + matrix(:,:,l)
    end do
  end subroutine

  !! Zero-out the values of the specified row and column elements.
  !! NB: Assumes symmetric structure
  subroutine project_out(this, index)
    class(pbsr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    INSIST(associated(this%graph%row_imap, this%graph%col_imap))
    m = index
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%values(:,:,lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      lnm = this%graph%index(n, m)
      ASSERT(lnm > 0)
      this%values(:,:,lnm) = 0.0_r8
    end do
  end subroutine

  !! Returns the matrix-vector product with vector X.
  subroutine matvec(this, x, y, incr)

    class(pbsr_matrix), intent(in) :: this
    real(r8), intent(in) :: x(:,:)
    real(r8), intent(inout) :: y(:,:)
    logical, intent(in), optional :: incr

    integer :: i, k
    real(r8) :: s(this%bsize)
    logical :: incr_

    ASSERT(size(x,1) == this%bsize)
    ASSERT(size(x,2) == this%ncol)
    ASSERT(size(y,1) == this%bsize)
    ASSERT(size(y,2) >= this%nrow_onP)

    incr_ = .false.
    if (present(incr)) incr_ = incr

    do i = 1, this%nrow_onP
      if (incr_) then
        s = y(:,i)
      else
        s = 0.0_r8
      end if
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + matmul(this%values(:,:,k), x(:,this%graph%adjncy(k)))
      end do
      y(:,i) = s
    end do
  end subroutine

end module pbsr_matrix_type
