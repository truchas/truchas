#include "f90_assert.fpp"

module block_linear_solver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: fct, slv

contains

  !! LU factorization of a square matrix.  No pivoting (intentionally).
  !! Unit lower triangular factor.

  subroutine fct(a)

    real(r8), intent(inout) :: a(:,:)

    integer :: i, j, k
    real(r8) :: lkk, lkj, ujk

    ASSERT(size(a,1) == size(a,2))

    do k = 2, size(a,1)
      lkk = a(k,k)
      do j = 1, k - 1
        lkj = a(k,j)
        ujk = a(j,k)
        do i = 1, j - 1
          lkj = lkj - a(k,i)*a(i,j)
          ujk = ujk - a(j,i)*a(i,k)
        end do
        lkj = lkj / a(j,j)
        lkk = lkk - lkj*ujk
        a(k,j) = lkj
        a(j,k) = ujk
      end do
      a(k,k) = lkk
    end do

  end subroutine fct


  subroutine slv(a, b)

    real(r8), intent(in) :: a(:,:)
    real(r8), intent(inout) :: b(:)

    integer  :: n, i, j
    real(r8) :: bj

    ASSERT(size(a,1) == size(a,2))
    ASSERT(size(b) == size(a,1))

    n = size(a,1)
    do j = 2, n
      bj = b(j)
      do i = 1, j-1
        bj = bj - a(j,i)*b(i)
      end do
      b(j) = bj
    end do
    b(n) = b(n) / a(n,n)
    do j = n-1, 1, -1
      bj = b(j)
      do i = j+1, n
        bj = bj - a(j,i)*b(i)
      end do
      b(j) = bj / a(j,j)
    end do

  end subroutine slv

end module block_linear_solver
