!!
!! BSR_MATRIX_TYPE
!!
!! This module defines a data structure for sparse matrices stored in block
!! compressed sparse row (BSR) format.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module bsr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use csr_graph_type, only: csr_graph
  implicit none
  private

  public :: csr_graph

  type, public :: bsr_matrix
    integer :: nrow, ncol, bsize
    real(r8), allocatable :: values(:,:,:)
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
    final :: bsr_matrix_delete
  end type

  public :: gs_relaxation
  interface gs_relaxation
    module procedure gs_relaxation_bsr
    !module procedure test_gs_relaxation
  end interface

contains

  !! Final subroutine for BSR_MATRIX type objects.
  elemental subroutine bsr_matrix_delete(this)
    type(bsr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine init_graph(this, bsize, graph, take_graph)
    class(bsr_matrix), intent(out) :: this
    integer, intent(in) :: bsize
    type(csr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%nrow = graph%nrow
    this%ncol = graph%ncol
    this%bsize = bsize
    allocate(this%values(bsize,bsize,size(graph%adjncy)), source=0.0_r8)
    this%graph => graph
    if (present(take_graph)) then
      this%graph_is_owned = take_graph
      if (take_graph) nullify(graph)
    end if
  end subroutine

  !! Initialize the sparse matrix using the non-zero graph structure from MOLD.
  subroutine init_mold(this, mold)
    class(bsr_matrix), intent(out) :: this
    class(bsr_matrix), intent(in)  :: mold
    call init_graph(this, mold%bsize, mold%graph, take_graph=.false.)
  end subroutine

  !! Set all matrix elements to the given value.
  subroutine set_all(this, value)
    class(bsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    this%values = value
  end subroutine

  !! Set the value of the specified matrix element to the given value.
  subroutine set(this, row, col, value)
    class(bsr_matrix), intent(inout) :: this
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
    class(bsr_matrix), intent(inout) :: this
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
    class(bsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:,:,:)
    integer :: i, j, n
    INSIST(this%nrow == this%ncol)
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
    class(bsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:,:)
    integer :: i, j, l, n
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
    class(bsr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    INSIST(this%nrow == this%ncol)
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
  function matvec(this, x) result(y)
    class(bsr_matrix), intent(in) :: this
    real(r8), intent(in) :: x(:,:)
    real(r8) :: y(this%bsize,this%nrow)
    integer :: i, k
    real(r8) :: s(this%bsize)
    ASSERT(size(x,1) == this%bsize)
    ASSERT(size(x,2) == this%ncol)
    do i = 1, this%nrow
      s = 0.0_r8
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + matmul(this%values(:,:,k), x(:,this%graph%adjncy(k)))
      end do
      y(:,i) = s
    end do
  end function

  subroutine kdiag_init(this)
    class(bsr_matrix), intent(inout) :: this
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

  !! NB: hard coded for 2x2 blocks. This does the standard GS
  !! algorithm with the 2x2 blocks regarded as the "numbers".
  !! The original does standard GS but with the block structure
  !! flattened.

  subroutine test_gs_relaxation(a, f, u, pattern)

    type(bsr_matrix), intent(inout) :: a
    real(r8), intent(in) :: f(:,:)
    real(r8), intent(inout) :: u(:,:)
    character(len=*), intent(in) :: pattern

    integer :: i, i1, i2, di, j, k, ld
    real(r8) :: s(a%bsize)

    ASSERT(a%nrow == a%ncol)
    ASSERT(size(f,1) == a%bsize)
    ASSERT(size(f,2) <= a%nrow)
    ASSERT(size(u,1) == a%bsize)
    ASSERT(size(u,2) >= a%ncol)

    if (.not.allocated(a%kdiag)) call a%kdiag_init

    do j = 1, len(pattern)
      call loop_range(pattern(j:j), size(f,2), i1, i2, di)
      do i = i1, i2, di
        s = f(:,i)
        do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
          s = s - matmul(a%values(:,:,k), u(:,a%graph%adjncy(k)))
        end do
        u(:,i) = u(:,i) + matmul(inverse(a%values(:,:,a%kdiag(i))), s)
      end do
    end do

  contains

    function inverse(a) result(inv_a)
      real(r8), intent(in) :: a(:,:)
      real(r8) :: inv_a(2,2)
      inv_a(1,1) = a(2,2)
      inv_a(2,2) = a(1,1)
      inv_a(1,2) = -a(1,2)
      inv_a(2,1) = -a(2,1)
      inv_a = inv_a / (a(1,1)*a(2,2)-a(2,1)*a(1,2))
    end function

  end subroutine test_gs_relaxation

  subroutine gs_relaxation_bsr(a, f, u, pattern)

    type(bsr_matrix), intent(inout) :: a
    real(r8), intent(in) :: f(:,:)
    real(r8), intent(inout) :: u(:,:)
    character(len=*), intent(in) :: pattern

    integer :: i, i1, i2, di, j, k, m1, m2, dm
    real(r8) :: s(a%bsize)

    ASSERT(a%nrow == a%ncol)
    ASSERT(size(f,1) == a%bsize)
    ASSERT(size(f,2) <= a%nrow)
    ASSERT(size(u,1) == a%bsize)
    ASSERT(size(u,2) >= a%ncol)

    if (.not.allocated(a%kdiag)) call a%kdiag_init

    do j = 1, len(pattern)
      call loop_range(pattern(j:j), size(f,2), i1, i2, di)
      call loop_range(pattern(j:j), a%bsize, m1, m2, dm)
      do i = i1, i2, di
        s = f(:,i)
        do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
          if (a%graph%adjncy(k) == i) cycle
          s = s - matmul(a%values(:,:,k), u(:,a%graph%adjncy(k)))
        end do
        call dense_gs(m1, m2, dm, a%values(:,:,a%kdiag(i)), s, u(:,i))
      end do
    end do

  contains

    pure subroutine dense_gs(i1, i2, di, a, f, u)
      integer, intent(in) :: i1, i2, di
      real(r8), intent(in) :: a(:,:), f(:)
      real(r8), intent(inout) :: u(:)
      integer :: i, j
      real(r8) :: s
      do i = i1, i2, di
        s = f(i)
        do j = 1, size(a,2)
          s = s - a(i,j)*u(j)
        end do
        u(i) = u(i) + s / a(i,i)
      end do
    end subroutine

  end subroutine gs_relaxation_bsr

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

end module bsr_matrix_type
