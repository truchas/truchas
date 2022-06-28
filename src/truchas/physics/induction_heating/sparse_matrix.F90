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
  implicit none
  private
  
  type, public :: msr_matrix
    integer :: nrow, ncol
    integer, allocatable :: xadj(:), adjncy(:)
    real(kind=r8), allocatable :: diag(:), nonz(:)
  end type msr_matrix
  
  public :: create_msr_matrix, destroy_msr_matrix, increment_msr_matrix_value
  public :: matmul, matmul_msr, gs_relaxation
  
  !! Debugging crap
  public :: is_symmetric, matmul_transpose
  
  interface matmul
    module procedure matmul_msr
  end interface
  
  interface gs_relaxation
    module procedure gs_relaxation_msr
  end interface
  
contains

  function matmul_msr (a, u, mask) result (v)
  
    type(msr_matrix), intent(in) :: a
    real(kind=r8), intent(in) :: u(:)
    logical, intent(in), optional :: mask(:)
    real(kind=r8) :: v(a%nrow)
    
    integer :: i, k
    real(kind=r8) :: value
    
    ASSERT( size(u) == a%ncol )
    
    if (present(mask)) then
    
      ASSERT( size(mask) == size(v) )
      do i = 1, a%nrow
        if (mask(i)) then
          value = a%diag(i) * u(i)
          do k = a%xadj(i), a%xadj(i+1)-1
            value = value + a%nonz(k) * u(a%adjncy(k))
          end do
          v(i) = value
        else
          v(i) = 0.0_r8
        end if
      end do
    
    else
    
      do i = 1, a%nrow
        value = a%diag(i) * u(i)
        do k = a%xadj(i), a%xadj(i+1)-1
          value = value + a%nonz(k) * u(a%adjncy(k))
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
    real(kind=r8),    intent(in)    :: f(:)
    real(kind=r8),    intent(inout) :: u(:)
    character(len=*), intent(in)    :: pattern
    
    integer :: i, j, k
    real(kind=r8) :: s
    
    ASSERT( a%nrow == a%ncol )
    ASSERT( size(f) <= a%nrow )
    ASSERT( size(u) >= a%ncol )
    
    do j = 1, len(pattern)
      select case (pattern(j:j))
      case ('f', 'F') ! FORWARD SWEEP
        do i = 1, size(f)
          s = f(i)
          do k = a%xadj(i), a%xadj(i+1)-1
            s = s - a%nonz(k) * u(a%adjncy(k))
          end do
          u(i) = s / a%diag(i)
        end do
        
      case ('b', 'B') ! BACKWARD SWEEP
        do i = size(f), 1, -1
          s = f(i)
          do k = a%xadj(i), a%xadj(i+1)-1
            s = s - a%nonz(k) * u(a%adjncy(k))
          end do
          u(i) = s / a%diag(i)
        end do
      end select
    end do

  end subroutine gs_relaxation_msr
  
  subroutine create_msr_matrix (n, clique_array, a)
  
    use graph_type
    
    integer, intent(in) :: n  ! number of rows/cols
    integer, intent(in) :: clique_array(:,:)
    type(msr_matrix), intent(out) :: a
    
    integer :: j
    type(graph), allocatable :: g
    
    ASSERT(size(clique_array)==0.or.minval(clique_array)==1)
    !ASSERT(size(clique_array)==0.or.maxval(clique_array)==n)
    
    a%nrow = n
    a%ncol = n
    
    !! Create a graph of the matrix nonzero structure
    allocate(g)
    call g%init (n)
    do j = 1, size(clique_array,dim=2)
      call g%add_clique (clique_array(:,j))
    end do
    
    !! Extract the adjacency structure and toss the graph.
    call g%get_adjacency (a%xadj, a%adjncy)
    deallocate(g)
    
    !! Allocate storage for the matrix elements and initialize.
    allocate(a%diag(n), a%nonz(size(a%adjncy)))
    a%diag = 0.0_r8
    a%nonz = 0.0_r8
    
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
  
  
  subroutine destroy_msr_matrix (a)
    type(msr_matrix), intent(inout) :: a
    deallocate(a%xadj, a%adjncy, a%diag, a%nonz)
    a%nrow = 0
    a%ncol = 0
  end subroutine destroy_msr_matrix
    

  subroutine increment_msr_matrix_value (a, row, col, value)
  
    type(msr_matrix), intent(inout) :: a
    integer, intent(in) :: row, col
    real(kind=r8), intent(in) :: value
    
    integer :: j, k

    if (row > a%nrow) return
    if (col > a%ncol) return

    ASSERT( row > 0 .or. row <= a%nrow )
    ASSERT( col > 0 .or. col <= a%ncol )
    
    if (row == col) then  ! diagonal element
      a%diag(row) = a%diag(row) + value
      
    else  ! off-diagonal element
      j = 0   ! scan for location
      do k = a%xadj(row), a%xadj(row+1)-1
        if (a%adjncy(k) == col) then
          j = k
          exit
        end if
      end do
      if (j /= 0) then
        a%nonz(j) = a%nonz(j) + value
      else
        print *, 'increment_msr_matrix_value: PANIC! no place for element ', row, col
        stop
      end if
    end if
    
    !! May be faster to use a version that assembled the values from a local cell
    !! matrix, rather than calling this one for each and every element.
    
  end subroutine increment_msr_matrix_value
    
  !!
  !! Debugging crap
  !!
  
  logical function is_symmetric (a)
    type(msr_matrix), intent(in) :: a
    logical :: checked(size(a%nonz))
    integer :: i, k, l
    checked = .false.
    do i = 1, a%nrow
      list: do k = a%xadj(i), a%xadj(i+1)-1
        if (checked(k)) cycle
        do l = a%xadj(a%adjncy(k)), a%xadj(a%adjncy(k)+1)-1
          if (a%adjncy(l) == i) then
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
    real(kind=r8), intent(in) :: u(:)
    real(kind=r8) :: v(a%ncol)
    integer :: i, k
    v = a%diag * u
    do i = 1, a%nrow
      do k = a%xadj(i), a%xadj(i+1)-1
        v(a%adjncy(k)) = v(a%adjncy(k)) + u(i) * a%nonz(k)
      end do
    end do
  end function matmul_transpose
  
end module sparse_matrix
