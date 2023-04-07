!!
!! COMPLEX_PCSR_MATRIX_TYPE
!!
!! This module defines a distributed data structure for sparse complex-valued
!! matrices stored in compressed sparse row (CSR) format.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module complex_pcsr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use csr_graph_type, only: pcsr_graph
  implicit none
  private

  public :: pcsr_graph

  type, public :: complex_pcsr_matrix
    integer :: nrow = 0, ncol = 0, nrow_onP = 0, ncol_onP = 0
    complex(r8), allocatable :: values(:)
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
    final :: complex_pcsr_matrix_delete
  end type

  !public :: gs_relaxation
  interface gs_relaxation
    module procedure gs_relaxation
  end interface

contains

  !! Final subroutine for COMPLEX_PCSR_MATRIX type objects.
  subroutine complex_pcsr_matrix_delete(this)
    type(complex_pcsr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine init_graph(this, graph, take_graph)
    class(complex_pcsr_matrix), intent(out) :: this
    type(pcsr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%nrow = graph%row_imap%local_size
    this%ncol = graph%col_imap%local_size
    this%nrow_onP = graph%row_imap%onP_size
    this%ncol_onP = graph%col_imap%onP_size
    allocate(this%values(size(graph%adjncy)), source=cmplx(0,0,kind=r8))
    this%graph => graph
    if (present(take_graph)) then
      this%graph_is_owned = take_graph
      if (take_graph) nullify(graph)
    end if
  end subroutine

  !! Initialize the sparse matrix using the non-zero structure from MOLD.
  subroutine init_mold(this, mold)
    class(complex_pcsr_matrix), intent(out) :: this
    class(complex_pcsr_matrix), intent(in)  :: mold
    call init_graph(this, mold%graph, take_graph=.false.)
  end subroutine

  !! Set all matrix elements to the given value.
  subroutine set_all(this, value)
    class(complex_pcsr_matrix), intent(inout) :: this
    complex(r8), intent(in) :: value
    ASSERT(allocated(this%values))
    this%values = value
  end subroutine

  !! Set the specified matrix element to the given value.
  subroutine set(this, row, col, value)
    class(complex_pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    complex(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    n = this%graph%index(row, col)
    ASSERT(n /= 0)
    this%values(n) = value
  end subroutine

  !! Increment the specified matrix element by the given value.
  subroutine add_to_element(this, row, col, value)
    class(complex_pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    complex(r8), intent(in) :: value
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
    class(complex_pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    complex(r8), intent(in) :: matrix(:,:)
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
    class(complex_pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    complex(r8), intent(in) :: matrix(:)
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
    class(complex_pcsr_matrix), intent(inout) :: this
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

    class(complex_pcsr_matrix), intent(in) :: this
    complex(r8), intent(in) :: x(:)
    complex(r8), intent(inout) :: b(:)
    logical, intent(in), optional :: incr

    integer :: i, k
    complex(r8) :: s
    logical :: incr_

    ASSERT(size(x) == this%ncol)
    ASSERT(size(b) >= this%nrow_onP)

    incr_ = .false.
    if (present(incr)) incr_ = incr

    do i = 1, this%nrow_onP
      s = merge(b(i), cmplx(0,0,kind=r8), incr_)
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + this%values(k) * x(this%graph%adjncy(k))
      end do
      b(i) = s
    end do

  end subroutine matvec

  subroutine kdiag_init(this)
    class(complex_pcsr_matrix), intent(inout) :: this
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
    class(complex_pcsr_matrix), intent(in), target :: this
    integer, intent(in) :: row
    complex(r8), pointer :: values(:)
    integer, pointer :: indices(:)
    ASSERT(row >= 1 .and. row <= this%nrow)
    values => this%values(this%graph%xadj(row):this%graph%xadj(row+1)-1)
    indices => this%graph%adjncy(this%graph%xadj(row):this%graph%xadj(row+1)-1)
  end subroutine

  !!
  !! Gauss-Seidel relaxation
  !!

  subroutine gs_relaxation(a, f, u, pattern)

    type(complex_pcsr_matrix), intent(inout) :: a
    complex(r8), intent(in) :: f(:)
    complex(r8), intent(inout) :: u(:)
    character(*), intent(in) :: pattern

    integer :: i, i1, i2, di, j, k, n
    complex(r8) :: s

    ASSERT(a%nrow == a%ncol)
    ASSERT(size(u) >= a%ncol)

    n = min(a%nrow, size(f))

    if (.not.allocated(a%kdiag)) call a%kdiag_init

    do j = 1, len(pattern)
      call loop_range(pattern(j:j), n, i1, i2, di)
      do i = i1, i2, di
        s = f(i)
        do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
          s = s - a%values(k) * u(a%graph%adjncy(k))
        end do
        u(i) = u(i) + s / a%values(a%kdiag(i))
      end do
    end do

  end subroutine gs_relaxation

  subroutine loop_range(direction, len, i1, i2, di)
    character(1), intent(in) :: direction
    integer, intent(in) :: len
    integer, intent(out) :: i1, i2, di
    select case (direction)
    case ('f', 'F') ! forward sweep
      i1 = 1
      i2 = len
      di = 1
    case ('b', 'B') ! backward sweep
      i1 = len
      i2 = 1
      di = -1
    case default
      INSIST(.false.)
    end select

  end subroutine

end module complex_pcsr_matrix_type
