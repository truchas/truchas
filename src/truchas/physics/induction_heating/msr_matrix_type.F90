!!
!! Neil Carlson <nnc@newmexico.com>
!! Last revised 4 Apr 2004
!!
!! This module encapsulates just the essential sparse matrix bits needed
!! to implement the Hiptmair relaxation procedure used as a CG preconditioner.
!! This in NOT a general sparse matrix module.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sparse_matrix

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use graph_type
  implicit none
  private

  type, public :: msr_graph
    integer :: nrow
    integer, allocatable :: xadj(:), adjncy(:)
    type(graph), allocatable, private :: g  ! temporary during construction
  contains
    procedure :: init => msr_graph_init
    generic :: add_clique => add_clique_one, add_clique_many
    procedure, private :: add_clique_one, add_clique_many
    procedure :: add_complete => add_complete
    procedure :: defined => msr_graph_defined
    procedure :: index => index_linear ! alternatively, index_binary
  end type

  type, public :: msr_matrix
    integer :: nrow, ncol
    real(r8), allocatable :: diag(:), nonz(:)
    type(msr_graph), pointer :: graph => null()
    logical, private :: graph_is_owned = .false.
  contains
    generic :: init => msr_matrix_init, msr_matrix_init_mold
    procedure, private :: msr_matrix_init, msr_matrix_init_mold
    procedure :: add_to
    procedure :: project_out
    procedure :: matvec
    procedure :: is_symmetric
    final :: msr_matrix_delete
  end type

  public :: gs_relaxation
  interface gs_relaxation
    module procedure gs_relaxation_msr
  end interface

contains

  subroutine msr_graph_init(this, n)
    class(msr_graph), intent(out) :: this
    integer, intent(in) :: n
    this%nrow = n
    allocate(this%g)
    call this%g%init(this%nrow)
  end subroutine

  subroutine add_clique_one(this, indices)
    class(msr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine add_clique_many(this, indices)
    class(msr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:,:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine add_complete(this)
    class(msr_graph), intent(inout) :: this
    ASSERT(allocated(this%g))
    call this%g%get_adjacency(this%xadj, this%adjncy)
    deallocate(this%g)
  end subroutine

  pure logical function msr_graph_defined(this)
    class(msr_graph), intent(in) :: this
    msr_graph_defined = allocated(this%xadj)
  end function

  !! Return the ADJNCY array index corresponding to matrix element (row, col),
  !! or 0 if that element is not stored by the sparse matrix. This version does
  !! a linear search of the appropriate sorted segment of the array.

  pure integer function index_linear(this, row, col) result(k)
    class(msr_graph), intent(in) :: this
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
    class(msr_graph), intent(in) :: this
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

  !! Final subroutine for MSR_MATRIX type objects.
  subroutine msr_matrix_delete(this)
    type(msr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine msr_matrix_init(this, graph, take_graph)
    class(msr_matrix), intent(out) :: this
    type(msr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    integer :: n
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%graph => graph
    if (present(take_graph)) this%graph_is_owned = take_graph
    n = graph%nrow
    this%nrow = n
    this%ncol = n
    allocate(this%diag(n), this%nonz(size(this%graph%adjncy)))
    this%diag = 0.0_r8  !TODO: make caller set initial values?
    this%nonz = 0.0_r8
  end subroutine

  !! Initialize the sparse matrix using the non-zero graph structure from MOLD.
  subroutine msr_matrix_init_mold(this, mold)
    class(msr_matrix), intent(out) :: this
    class(msr_matrix), intent(in)  :: mold
    call msr_matrix_init(this, mold%graph, take_graph=.false.)
  end subroutine

  !! Increment the specified matrix element by the given value.
  subroutine add_to(this, row, col, value)
    class(msr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    !if (row > this%nrow) return  !TODO: can we delete this? is it used?
    !if (col > this%ncol) return
    ASSERT(row >= 1 .or. row <= this%nrow)
    if (row == col) then  ! diagonal element
      this%diag(row) = this%diag(row) + value
    else  ! off-diagonal element
      n = this%graph%index(row,col)
      ASSERT(n /= 0)
      this%nonz(n) = this%nonz(n) + value
    end if
  end subroutine

  !! Zero-out the values of the specified row and column elements.
  !! NB: Assumes symmetric structure
  subroutine project_out(this, index)
    class(msr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    m = index
    this%diag(m) = 0.0_r8
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%nonz(lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      lnm = this%graph%index(n, m)
      ASSERT(lnm > 0)
      this%nonz(lnm) = 0.0_r8
    end do
  end subroutine

  !! Returns the matrix-vector product with vector X.
  function matvec(this, x) result(y)
    class(msr_matrix), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: y(size(x))
    integer :: i, k
    real(r8) :: s
    ASSERT(size(x) == this%nrow)
    do i = 1, this%nrow
      s = this%diag(i) * x(i)
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + this%nonz(k) * x(this%graph%adjncy(k))
      end do
      y(i) = s
    end do
  end function

  !! Returns true if the matrix is exactly symmetric (bit-for-bit); otherwise false.
  logical function is_symmetric(this)
    class(msr_matrix), intent(in) :: this
    logical :: checked(size(this%nonz))
    integer :: m, n, lmn, lnm
    is_symmetric = .false.
    checked = .false.
    do m = 1, this%nrow
      do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
        if (checked(lmn)) cycle
        n = this%graph%adjncy(lmn)
        lnm = this%graph%index(n, m)
        if (lnm == 0) return ! non-symmetric structure
        if (this%nonz(lnm) /= this%nonz(lmn)) return
        checked(lmn) = .true.
        checked(lnm) = .true.
      end do
    end do
    is_symmetric = .true.
  end function

  !!
  !! Gauss-Seidel relaxation
  !!

  subroutine gs_relaxation_msr(a, f, u, pattern)

    type(msr_matrix), intent(in) :: a
    real(r8), intent(in) :: f(:)
    real(r8), intent(inout) :: u(:)
    character(len=*), intent(in) :: pattern

    integer :: i, j, k
    real(r8) :: s

    ASSERT(a%nrow == a%ncol)
    ASSERT(size(f) <= a%nrow)
    ASSERT(size(u) >= a%ncol)

    do j = 1, len(pattern)
      select case (pattern(j:j))
      case ('f', 'F') ! FORWARD SWEEP
        do i = 1, size(f)
          s = f(i)
          do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
            s = s - a%nonz(k) * u(a%graph%adjncy(k))
          end do
          u(i) = s / a%diag(i)
        end do

      case ('b', 'B') ! BACKWARD SWEEP
        do i = size(f), 1, -1
          s = f(i)
          do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
            s = s - a%nonz(k) * u(a%graph%adjncy(k))
          end do
          u(i) = s / a%diag(i)
        end do
      end select
    end do

  end subroutine gs_relaxation_msr

end module sparse_matrix
