!!
!! CGSolver
!!
!! Module provides a basic PCG iterative solver for use by the electromagnetic
!! solver only.  It is intended to be a temporary solution, until UbikSolve is
!! extended to handle the EM system and its specialized preconditioner.
!!
!! Neil N. Carlson <nnc@newmexico.com>
!!
!! PARALLELIZATION NOTES
!!
!! o All vectors are extended; that is, off-PE data follows the on-PE data.
!! o Procedures require that the off-PE data in input vectors is correct.
!! o Procedures must ensure that the off-PE data in output vectors and in
!!   vectors passed to other procedures is correct.
!! o Although the actual system does not include the off-PE data, the CG solver
!!   carries along this data because the procedures it interfaces with need it.
!!   This is not an issue except when computing the global inner product, where
!!   the off-PE data must be excluded.   The argument DIM specifies the true
!!   dimension of the system; that is, the size of the on-PE data which form
!!   the first part of the extended vectors.
!!

#include "f90_assert.fpp"

module CGSolver

  use kinds, only: r8
  use parallel_communication
  use EM_utilities
  !use index_partitioning
  implicit none
  private
  
  public :: SolveCG
    
  interface SolveCG
    module procedure solve_cg_1
  end interface
  
  type, public :: cg_desc
    integer :: minitr, maxitr ! min and max number of CG iterations
    integer :: output_level
    real(kind=r8) :: tol  ! CG stopping tolerance on the residual error
    real(kind=r8) :: red  ! CG stopping tolerance on the residual error reduction
    !! Pointers to the system matvec and PC procedures belong here...
  end type cg_desc
  
contains

  subroutine solve_cg_1 (cg, ax, pc, b, x, dim, stat)!, gsd)
  
    type(cg_desc), intent(in) :: cg   ! CG control descriptor
    real(kind=r8), intent(in) :: b(:) ! Equation RHS
    real(kind=r8), intent(inout) :: x(:) ! Solution (initial guess/final)
    integer, intent(in)  :: dim   ! True dimension of the system
    integer, intent(out) :: stat
    !type(gs_desc), intent(in) :: gsd
    
    interface
      subroutine ax (x, y)
        use kinds, only: r8
        real(kind=r8), intent(in)  :: x(:)
        real(kind=r8), intent(out) :: y(:)
      end subroutine ax
      subroutine pc (x, y)
        use kinds, only: r8
        real(kind=r8), intent(in)  :: x(:)
        real(kind=r8), intent(out) :: y(:)
      end subroutine pc
    end interface
    
    real(kind=r8) :: r(size(x)), p(size(x)), q(size(x))
    real(kind=r8) :: rho, rho_init, rho_last, s
    integer :: itr, n, j
    character(len=80) :: message
    
    ASSERT( size(x) == size(b) )
    ASSERT( dim <= size(x) )
    ASSERT( size(x) == 0 .or. dim >= 1)
    
    n = size(b)
    
    !r = b - ax(x)
    call ax (x, r)
    do j = 1, n
      r(j) = b(j) - r(j)
    end do
    
    !p = pc(r)
    call pc (r, p)
    
    rho = global_dot_product(r(:dim), p(:dim))
    rho_init = rho
    INSIST( rho_init >= 0.0_r8 )
    stat = -1
    
    if (rho == 0.0_r8) then
      stat = 0
      return
    end if
    
    do itr = 1, cg%maxitr
    
      !q = ax(p)
      call ax (p, q)
      s = rho / global_dot_product(p(:dim), q(:dim))
      INSIST( s >= 0.0_r8 )
      do j = 1, n
        x(j) = x(j) + s * p(j)
        r(j) = r(j) - s * q(j)
      end do
      
      !q = pc(r)
      call pc (r, q)
      rho_last = rho
      rho = global_dot_product(r(:dim), q(:dim))
      INSIST( rho >= 0.0_r8 )
      s = rho / rho_last
      do j = 1, n
        p(j) = q(j) + s * p(j)
      end do
      
      if (cg%output_level >= 4) then
        write(message,fmt="(t8,a,i4,2(a,es10.3))") 'cg iterate', itr, &
            ': |r|=', sqrt(rho), ', |r|/|r_prev|=', sqrt(s)
        call EM_info (message)
      end if
      
      if (itr >= cg%minitr .and. (rho < cg%tol**2 .or. rho < cg%red**2 * rho_init)) then
        stat = 0
        exit
      end if
      
    end do
    
    if (cg%output_level >= 2) then
      write(message,fmt="(t6,a,i4,2(a,es10.3))") 'step cg summary:', itr, &
          ' iterations, |r|=', sqrt(rho), ', |r0|=', sqrt(rho_init)
      call EM_info (message)
    end if
    
  end subroutine solve_cg_1
  
end module CGSolver
