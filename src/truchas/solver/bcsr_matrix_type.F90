!!
!! BCSR_MATRIX_TYPE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

module bcsr_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use graph_type
  implicit none
  private

  type, public :: csr_graph
    integer :: nrow
    integer, allocatable :: xadj(:), adjncy(:)
    type(graph), allocatable, private :: g  ! temporary during construction
  contains
    procedure :: init => csr_graph_init
    generic :: add_clique => add_clique_one, add_clique_many
    procedure, private :: add_clique_one, add_clique_many
    procedure :: add_complete => add_complete
    procedure :: defined => csr_graph_defined
    procedure :: index => index_linear ! alternatively, index_binary
  end type

  type, public :: bcsr_matrix
    integer :: nrow, ncol, bsize
    real(r8), allocatable :: diag(:,:,:), nonz(:,:,:)
    type(csr_graph), pointer :: graph => null()
    logical, private :: graph_is_owned = .false.
  contains
    generic :: init => bcsr_matrix_init, bcsr_matrix_init_mold
    procedure, private :: bcsr_matrix_init, bcsr_matrix_init_mold
    procedure :: set_all
    procedure :: set
    generic :: add_to => add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure, private :: add_to_element, add_to_submatrix, add_to_submatrix_upm
    procedure :: project_out
    procedure :: matvec
    !procedure :: is_symmetric
    final :: bcsr_matrix_delete
  end type

  public :: gs_relaxation
  interface gs_relaxation
    module procedure gs_relaxation_msr
  end interface

contains

  subroutine csr_graph_init(this, n)
    class(csr_graph), intent(out) :: this
    integer, intent(in) :: n
    this%nrow = n
    allocate(this%g)
    call this%g%init(this%nrow)
  end subroutine

  subroutine add_clique_one(this, indices)
    class(csr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine add_clique_many(this, indices)
    class(csr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:,:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine add_complete(this)
    class(csr_graph), intent(inout) :: this
    ASSERT(allocated(this%g))
    call this%g%get_adjacency(this%xadj, this%adjncy)
    deallocate(this%g)
  end subroutine

  pure logical function csr_graph_defined(this)
    class(csr_graph), intent(in) :: this
    csr_graph_defined = allocated(this%xadj)
  end function

  !! Return the ADJNCY array index corresponding to matrix element (row, col),
  !! or 0 if that element is not stored by the sparse matrix. This version does
  !! a linear search of the appropriate sorted segment of the array.

  pure integer function index_linear(this, row, col) result(k)
    class(csr_graph), intent(in) :: this
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
    class(csr_graph), intent(in) :: this
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
  subroutine bcsr_matrix_delete(this)
    type(bcsr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  !! Initialize the sparse matrix with non-zero structure GRAPH.
  subroutine bcsr_matrix_init(this, bsize, graph, take_graph)
    class(bcsr_matrix), intent(out) :: this
    integer, intent(in) :: bsize
    type(csr_graph), pointer :: graph
    logical, intent(in), optional :: take_graph
    integer :: n
    ASSERT(associated(graph))
    ASSERT(graph%defined())
    this%graph => graph
    if (present(take_graph)) this%graph_is_owned = take_graph
    n = graph%nrow
    this%nrow = n
    this%ncol = n
    this%bsize = bsize
    allocate(this%diag(bsize,bsize,n), this%nonz(bsize,bsize,size(this%graph%adjncy)))
    this%diag = 0.0_r8  !TODO: make caller set initial values?
    this%nonz = 0.0_r8
  end subroutine

  !! Initialize the sparse matrix using the non-zero graph structure from MOLD.
  subroutine bcsr_matrix_init_mold(this, mold)
    class(bcsr_matrix), intent(out) :: this
    class(bcsr_matrix), intent(in)  :: mold
    call bcsr_matrix_init(this, mold%bsize, mold%graph, take_graph=.false.)
  end subroutine

  !! Set all matrix elements to the given value.
  subroutine set_all(this, value)
    class(bcsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    this%diag = value
    this%nonz = value
  end subroutine

  !! Set the value of the specified matrix element to the given value.
  subroutine set(this, row, col, value)
    class(bcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value(:,:)
    integer :: n
    ASSERT(min(row,col) >= 1 .and. max(row,col) <= this%nrow)
    ASSERT(size(value,1) == this%bsize)
    ASSERT(size(value,2) == this%bsize)
    if (row == col) then
      this%diag(:,:,row) = value
    else
      n = this%graph%index(row,col)
      ASSERT(n /= 0)
      this%nonz(:,:,n) = value
    end if
  end subroutine

  !! Increment the specified matrix element by the given value.
  subroutine add_to_element(this, row, col, value)
    class(bcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value(:,:)
    integer :: n
    ASSERT(min(row,col) >= 1 .and. max(row,col) <= this%nrow)
    ASSERT(size(value,1) == this%bsize)
    ASSERT(size(value,2) == this%bsize)
    if (row == col) then  ! diagonal element
      this%diag(:,:,row) = this%diag(:,:,row) + value
    else  ! off-diagonal element
      n = this%graph%index(row,col)
      ASSERT(n /= 0)
      this%nonz(:,:,n) = this%nonz(:,:,n) + value
    end if
  end subroutine

  !! Increment the principle submatrix specified by the list of row indices
  !! with the corresponding values in the given matrix.
  subroutine add_to_submatrix(this, rows, matrix)
    class(bcsr_matrix), intent(inout) :: this
    integer, intent(in) :: rows(:)
    real(r8), intent(in) :: matrix(:,:,:,:)
    integer :: i, j, n
    ASSERT(size(matrix,1) == this%bsize)
    ASSERT(size(matrix,2) == this%bsize)
    ASSERT(size(matrix,3) == size(matrix,4))
    ASSERT(size(rows) == size(matrix,3))
    ASSERT(minval(rows) >= 1 .and. maxval(rows) <= this%nrow)
    do i = 1, size(rows)
      do j = 1, size(rows)
        if (i == j) then
          n = rows(i)
          this%diag(:,:,n) = this%diag(:,:,n) + matrix(:,:,i,i)
        else
          n = this%graph%index(rows(i),rows(j))
          ASSERT(n /= 0)
          this%nonz(:,:,n) = this%nonz(:,:,n) + matrix(:,:,i,j)
        end if
      end do
    end do
  end subroutine


  !! Increment the principle submatrix specified by the list of row indices
  !! with the corresponding values in the given symmetric matrix in upper
  !! packed storage format.
  subroutine add_to_submatrix_upm(this, rows, matrix)
    class(bcsr_matrix), intent(inout) :: this
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
        this%nonz(:,:,n) = this%nonz(:,:,n) + matrix(:,:,l)
        n = this%graph%index(rows(j),rows(i))
        this%nonz(:,:,n) = this%nonz(:,:,n) + matrix(:,:,l)
      end do
      l = l + 1
      n = rows(j)
      this%diag(:,:,n) = this%diag(:,:,n) + matrix(:,:,l)
    end do
  end subroutine

  !! Zero-out the values of the specified row and column elements.
  !! NB: Assumes symmetric structure
  subroutine project_out(this, index)
    class(bcsr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    m = index
    this%diag(:,:,m) = 0.0_r8
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%nonz(:,:,lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      lnm = this%graph%index(n, m)
      ASSERT(lnm > 0)
      this%nonz(:,:,lnm) = 0.0_r8
    end do
  end subroutine

  !! Returns the matrix-vector product with vector X.
  function matvec(this, x) result(y)
    class(bcsr_matrix), intent(in) :: this
    real(r8), intent(in) :: x(:,:)
    real(r8) :: y(size(x,1),size(x,2))
    integer :: i, k
    real(r8) :: s(this%bsize)
    ASSERT(size(x,2) == this%nrow)
    do i = 1, this%nrow
      s = matmul(this%diag(:,:,i), x(:,i))
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        s = s + matmul(this%nonz(:,:,k), x(:,this%graph%adjncy(k)))
      end do
      y(:,i) = s
    end do
  end function

!  !! Returns true if the matrix is exactly symmetric (bit-for-bit); otherwise false.
!  logical function is_symmetric(this)
!    class(bcsr_matrix), intent(in) :: this
!    logical :: checked(size(this%nonz))
!    integer :: m, n, lmn, lnm
!    is_symmetric = .false.
!    checked = .false.
!    do m = 1, this%nrow
!      do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
!        if (checked(lmn)) cycle
!        n = this%graph%adjncy(lmn)
!        lnm = this%graph%index(n, m)
!        if (lnm == 0) return ! non-symmetric structure
!        if (this%nonz(lnm) /= this%nonz(lmn)) return
!        checked(lmn) = .true.
!        checked(lnm) = .true.
!      end do
!    end do
!    is_symmetric = .true.
!  end function

  !!
  !! Gauss-Seidel relaxation
  !!

  subroutine gs_relaxation_msr(a, f, u, pattern)

    type(bcsr_matrix), intent(in) :: a
    real(r8), intent(in) :: f(:,:)
    real(r8), intent(inout) :: u(:,:)
    character(len=*), intent(in) :: pattern

    integer :: i, j, k
    real(r8) :: s(a%bsize)

    ASSERT(a%nrow == a%ncol)
    ASSERT(size(f,1) == a%bsize)
    ASSERT(size(f,2) <= a%nrow)
    ASSERT(size(u,1) == a%bsize)
    ASSERT(size(u,2) >= a%ncol)

    do j = 1, len(pattern)
      select case (pattern(j:j))
      case ('f', 'F') ! FORWARD SWEEP
        do i = 1, size(f,2)
          s = f(:,i)
          do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
            s = s - matmul(a%nonz(:,:,k), u(:,a%graph%adjncy(k)))
          end do
          call dense_fgs(a%diag(:,:,i), s, u(:,i))
        end do

      case ('b', 'B') ! BACKWARD SWEEP
        do i = size(f,2), 1, -1
          s = f(:,i)
          do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
            s = s - matmul(a%nonz(:,:,k), u(:,a%graph%adjncy(k)))
          end do
          call dense_bgs(a%diag(:,:,i), s, u(:,i))
        end do
      end select
    end do

  contains
  
    subroutine dense_fgs(a, f, u)
      real(r8), intent(in) :: a(:,:), f(:)
      real(r8), intent(inout) :: u(:)
      integer :: i, j
      real(r8) :: s
      do i = 1, size(a,1)
        s = f(i)
        do j = 1, size(a,2)
          if (j == i) cycle
          s = s - a(i,j)*u(j)
        end do
        u(i) = s / a(i,i)
      end do
    end subroutine
  
    subroutine dense_bgs(a, f, u)
      real(r8), intent(in) :: a(:,:), f(:)
      real(r8), intent(inout) :: u(:)
      integer :: i, j
      real(r8) :: s
      do i = size(a,1), 1, -1
        s = f(i)
        do j = 1, size(a,2)
          if (j == i) cycle
          s = s - a(i,j)*u(j)
        end do
        u(i) = s / a(i,i)
      end do
    end subroutine

  end subroutine gs_relaxation_msr

end module bcsr_matrix_type
