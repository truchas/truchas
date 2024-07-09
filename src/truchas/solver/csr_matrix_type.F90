!!
!! CSR_MATRIX_TYPE
!!
!! This module defines a data structure for sparse matrices stored in compressed
!! sparse row (CSR) format.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module csr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use csr_graph_type, only:csr_graph
  implicit none
  private

  public :: csr_graph

  type, public :: csr_matrix
    integer :: nrow, ncol
    real(r8), allocatable :: values(:)
    type(csr_graph), pointer :: graph => null()
    logical, private :: graph_is_owned = .false.
    integer, allocatable :: kdiag(:)
  contains
    generic :: init => init_graph, init_mold
    procedure, private :: init_graph, init_mold
    procedure :: set_all
    procedure :: set
    generic :: add_to => add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure, private :: add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure :: project_out
    procedure :: matvec
    procedure :: kdiag_init
    procedure :: is_symmetric
    final :: csr_matrix_delete
  end type

  public :: gs_relaxation
  interface gs_relaxation
    module procedure gs_relaxation
  end interface

contains

  !! Final subroutine for CSR_MATRIX type objects.
  subroutine csr_matrix_delete(this)
    type(csr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine init_graph(this, graph, take_graph)
    class(csr_matrix), intent(out) :: this
    type(csr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%nrow = graph%nrow
    this%ncol = graph%ncol
    allocate(this%values(size(graph%adjncy)), source=0.0_r8)
    this%graph => graph
    if (present(take_graph)) then
      this%graph_is_owned = take_graph
      if (take_graph) nullify(graph)
    end if
  end subroutine

  !! Initialize the sparse matrix using the non-zero graph structure from MOLD.
  subroutine init_mold(this, mold)
    class(csr_matrix), intent(out) :: this
    class(csr_matrix), intent(in)  :: mold
    call init_graph(this, mold%graph, take_graph=.false.)
  end subroutine

  !! Set all matrix elements to the given value.
  subroutine set_all(this, value)
    class(csr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    this%values = value
  end subroutine

  !! Set the value of the specified matrix element to the given value.
  subroutine set(this, row, col, value)
    class(csr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    n = this%graph%index(row,col)
    ASSERT(n /= 0)
    this%values(n) = value
  end subroutine

  !! Increment the specified matrix element by the given value.
  subroutine add_to_element(this, row, col, value)
    class(csr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    ASSERT(col >= 1 .and. col <= this%ncol)
    n = this%graph%index(row,col)
    ASSERT(n /= 0)
    this%values(n) = this%values(n) + value
  end subroutine

  !! Increment the principal submatrix specified by the list of row indices
  !! with the corresponding values in the given matrix.

  subroutine add_to_submatrix(this, rows, matrix)
    class(csr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:)
    integer :: i, j, n
    INSIST(this%nrow == this%ncol)
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
    class(csr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:)
    integer :: i, j, l, n
    INSIST(this%nrow == this%ncol)
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
    class(csr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    INSIST(this%nrow == this%ncol)
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

  !! Returns the matrix-vector product with vector X.
  function matvec(this, x) result(y)
    class(csr_matrix), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: y(this%nrow)
    integer :: i, k
    real(r8) :: s
    ASSERT(size(x) == this%ncol)
    do i = 1, this%nrow
      s = 0.0_r8
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + this%values(k) * x(this%graph%adjncy(k))
      end do
      y(i) = s
    end do
  end function

  !! Returns true if the matrix is exactly symmetric (bit-for-bit); otherwise false.
  logical function is_symmetric(this)
    class(csr_matrix), intent(in) :: this
    logical :: checked(size(this%values))
    integer :: m, n, lmn, lnm
    is_symmetric = .false.
    checked = .false.
    do m = 1, this%nrow
      do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
        if (checked(lmn)) cycle
        n = this%graph%adjncy(lmn)
        if (n == m) cycle ! diagonal element
        lnm = this%graph%index(n, m)
        if (lnm == 0) return ! non-symmetric structure
        if (this%values(lnm) /= this%values(lmn)) return
        checked(lmn) = .true.
        checked(lnm) = .true.
      end do
    end do
    is_symmetric = .true.
  end function

  subroutine kdiag_init(this)
    class(csr_matrix), intent(inout) :: this
    integer :: i
    if (allocated(this%kdiag)) return
    INSIST(this%nrow == this%ncol)
    allocate(this%kdiag(this%nrow))
    do i = 1, this%nrow
      this%kdiag(i) = this%graph%index(i,i)
    end do
    ASSERT(all(this%kdiag > 0))
  end subroutine

  !!
  !! Gauss-Seidel relaxation
  !!

  subroutine gs_relaxation(a, f, u, pattern)

    type(csr_matrix), intent(inout) :: a
    real(r8), intent(in) :: f(:)
    real(r8), intent(inout) :: u(:)
    character(*), intent(in) :: pattern

    integer :: i, i1, i2, di, j, k
    real(r8) :: s

    ASSERT(a%nrow == a%ncol)
    ASSERT(size(f) <= a%nrow)
    ASSERT(size(u) >= a%ncol)

    if (.not.allocated(a%kdiag)) call a%kdiag_init

    do j = 1, len(pattern)
      call loop_range(pattern(j:j), size(f), i1, i2, di)
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

end module csr_matrix_type
