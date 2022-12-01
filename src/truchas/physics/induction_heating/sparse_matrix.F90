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
    generic :: add_clique => msr_graph_add_clique_one, msr_graph_add_clique_many
    procedure, private :: msr_graph_add_clique_one, msr_graph_add_clique_many
    procedure :: add_complete => msr_graph_add_complete
  end type
  
  type, public :: msr_matrix
    integer :: nrow, ncol
    real(r8), allocatable :: diag(:), nonz(:)
    type(msr_graph), pointer :: graph => null()
    logical, private :: graph_is_owned = .false.
  contains
    procedure :: init => create_msr_matrix
    procedure :: add_to => increment_msr_matrix_value
    procedure :: matvec => matmul_msr
    procedure :: project_out => msr_matrix_project_out
    final :: msr_matrix_delete
  end type msr_matrix
  
  public :: gs_relaxation
  
  !! Debugging crap
  public :: is_symmetric, matmul_transpose
  
  interface gs_relaxation
    module procedure gs_relaxation_msr
  end interface
  
contains

  function matmul_msr(this, u, mask) result(v)
  
    class(msr_matrix), intent(in) :: this
    real(r8), intent(in) :: u(:)
    logical, intent(in), optional :: mask(:)
    real(r8) :: v(this%nrow)
    
    integer :: i, k
    real(r8) :: value
    
    ASSERT( size(u) == this%ncol )
    
    if (present(mask)) then
    
      ASSERT( size(mask) == size(v) )
      do i = 1, this%nrow
        if (mask(i)) then
          value = this%diag(i) * u(i)
          do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
            value = value + this%nonz(k) * u(this%graph%adjncy(k))
          end do
          v(i) = value
        else
          v(i) = 0.0_r8
        end if
      end do
    
    else
    
      do i = 1, this%nrow
        value = this%diag(i) * u(i)
        do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
          value = value + this%nonz(k) * u(this%graph%adjncy(k))
        end do
        v(i) = value
      end do
      
    end if
    
  end function matmul_msr
  
  !!
  !! Gauss-Seidel relaxation
  !!
  
  subroutine gs_relaxation_msr (a, f, u, pattern)
  
    type(msr_matrix), intent(in)    :: a
    real(r8),    intent(in)    :: f(:)
    real(r8),    intent(inout) :: u(:)
    character(len=*), intent(in)    :: pattern
    
    integer :: i, j, k
    real(r8) :: s
    
    ASSERT( a%nrow == a%ncol )
    ASSERT( size(f) <= a%nrow )
    ASSERT( size(u) >= a%ncol )
    
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
  
  subroutine msr_graph_init(this, n)
    class(msr_graph), intent(out) :: this
    integer, intent(in) :: n
    this%nrow = n
    allocate(this%g)
    call this%g%init(this%nrow)
  end subroutine
  
  subroutine msr_graph_add_clique_one(this, indices)
    class(msr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine msr_graph_add_clique_many(this, indices)
    class(msr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:,:)
    ASSERT(allocated(this%g))
    call this%g%add_clique(indices)
  end subroutine

  subroutine msr_graph_add_complete(this)
    class(msr_graph), intent(inout) :: this
    ASSERT(allocated(this%g))
    call this%g%get_adjacency(this%xadj, this%adjncy)
    deallocate(this%g)
  end subroutine

  !! Final subroutine for MSR_MATRIX type objects.
  subroutine msr_matrix_delete(this)
    type(msr_matrix), intent(inout) :: this
    if (this%graph_is_owned) then
      if (associated(this%graph)) deallocate(this%graph)
    end if
  end subroutine

  subroutine create_msr_matrix(this, n, clique_array)
  
    class(msr_matrix), intent(out) :: this
    integer, intent(in) :: n  ! number of rows/cols
    integer, intent(in) :: clique_array(:,:)
    
    ASSERT(size(clique_array)==0.or.minval(clique_array)==1)
    !ASSERT(size(clique_array)==0.or.maxval(clique_array)==n)
    
    this%nrow = n
    this%ncol = n
    
    allocate(this%graph)
    this%graph_is_owned = .true.
    call this%graph%init(n)
    call this%graph%add_clique(clique_array)
    call this%graph%add_complete
    
    !! Allocate storage for the matrix elements and initialize.
    allocate(this%diag(n), this%nonz(size(this%graph%adjncy)))
    this%diag = 0.0_r8
    this%nonz = 0.0_r8
    
    !! NB: In graph theory, a clique is a set of nodes such that there is
    !! an edge joining any two in the set.  In a node-based FE matrix, an
    !! element is nonzero only if the support of the corresponding basis
    !! functions overlap (row and col indices are the nodes associated with
    !! the basis functions).  Thus for typical FE methods the list of nodes
    !! of a cell are a clique in the associated matrix graph, and the entire
    !! cell node list (the clique_array argument) is sufficient to exactly
    !! construct the matrix graph.  Similarly, in an edge-based FE matrix
    !! the cell edge list will completely describe associated matrix graph
    !! whose nodes correspond to mesh edges, and whose edges correspond to
    !! coupling between mesh edges which occurs only through cells.                 
    
  end subroutine create_msr_matrix
  
  
  subroutine increment_msr_matrix_value(this, row, col, value)
  
    class(msr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    
    integer :: j, k

    if (row > this%nrow) return
    if (col > this%ncol) return

    ASSERT( row > 0 .or. row <= this%nrow )
    ASSERT( col > 0 .or. col <= this%ncol )
    
    if (row == col) then  ! diagonal element
      this%diag(row) = this%diag(row) + value
      
    else  ! off-diagonal element
      j = 0   ! scan for location

      do k = this%graph%xadj(row), this%graph%xadj(row+1)-1
        if (this%graph%adjncy(k) == col) then
          j = k
          exit
        end if
      end do
      if (j /= 0) then
        this%nonz(j) = this%nonz(j) + value
      else
        print *, 'increment_msr_matrix_value: PANIC! no place for element ', row, col
        stop
      end if
    end if
    
    !! May be faster to use this version that assembled the values from this local cell
    !! matrix, rather than calling this one for each and every element.
    
  end subroutine increment_msr_matrix_value

  !! Zero-out the values of the specified row and column elements.
  !! NB: Assumes symmetric structure
  subroutine msr_matrix_project_out(this, index)
    class(msr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    m = index
    this%diag(m) = 0.0_r8
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%nonz(lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      do lnm = this%graph%xadj(n), this%graph%xadj(n+1)-1
        if (this%graph%adjncy(lnm) == m) exit
      end do
      ASSERT(lnm < this%graph%xadj(n+1))
      this%nonz(lnm) = 0.0_r8
    end do
  end subroutine
    
  !!
  !! Debugging crap
  !!
  
  logical function is_symmetric (a)
    type(msr_matrix), intent(in) :: a
    logical :: checked(size(a%nonz))
    integer :: i, k, l
    checked = .false.
    do i = 1, a%nrow
      list: do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
        if (checked(k)) cycle
        do l = a%graph%xadj(a%graph%adjncy(k)), a%graph%xadj(a%graph%adjncy(k)+1)-1
          if (a%graph%adjncy(l) == i) then
            if (a%nonz(l) == a%nonz(k)) then
              checked(k) = .true.
              checked(l) = .true.
              cycle list
            else
              is_symmetric = .false.
              return
            end if
          end if
        end do
        is_symmetric = .false.
        return
      end do list
    end do
    is_symmetric = .true.
    if (any(.not.checked)) print *, 'some elements not checked !?!?'
  end function is_symmetric
  
  function matmul_transpose (a, u) result (v)
    type(msr_matrix), intent(in) :: a
    real(r8), intent(in) :: u(:)
    real(r8) :: v(a%ncol)
    integer :: i, k
    v = a%diag * u
    do i = 1, a%nrow
      do k = a%graph%xadj(i), a%graph%xadj(i+1)-1
        v(a%graph%adjncy(k)) = v(a%graph%adjncy(k)) + u(i) * a%nonz(k)
      end do
    end do
  end function matmul_transpose
  
end module sparse_matrix
