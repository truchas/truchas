!!
!! UPPER_PACKED_MATRIX_PROCS
!!
!! A collection of procedures that work with symmetric matrices stored in the
!! upper packed storage format, where the elements of the upper triangle are
!! stored contiguously by column in a rank-1 array.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 10 Jan 2006
!!
!! NB: Not all rank-1 arrays can be viewed as a symmetric matrix. Only arrays
!! of size 1, 3, 6, 10, 15, ... (for matrices of order 1, 2, 3, 4, 5, ...) are
!! valid arguments for a symmetric matrix argument. THIS IS NOT CHECKED. This
!! approach is a lightweight alternative to introducing a dedicated "upper
!! packed matrix" derived type that could ensure this constraint.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module upper_packed_matrix_procs

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: upm_factor, upm_solve, upm_invert, sym_matmul, upm_col_sum, upm_quad_form, upm_cong_prod

  interface sym_matmul
    procedure sym_matmul_vec_f77, sym_matmul_vec_f90
    procedure sym_matmul_mat_f77, sym_matmul_mat_f90
  end interface

  interface upm_quad_form
    procedure upm_quad_form_f77, upm_quad_form_f90, upm_zquad_form_f90
  end interface

  interface upm_inner_prod
    procedure :: upm_inner_prod_f77, upm_inner_prod_f90
  end interface

  interface upm_cong_prod
    procedure upm_cong_prod_f77, upm_cong_prod_f90
  end interface

contains

  !! Compute the Choleski factorization of a symmetric positive definite matrix.
  !! The matrix is provided in upper packed storage format by the argument A,
  !! and the upper Choleski factor (also in upper packed storage) is returned
  !! in A, overwriting the original matrix.

  pure subroutine upm_factor(a)

    use ieee_exceptions

    real(r8), intent(inout) :: a(:)

    integer :: l, m, p, q
    real(r8) :: s, t

    if (size(a) == 0) return

    !! Compute the upper Choleski factor.
    a(1) = sqrt(a(1))
    p = 2
    do while (p <= size(a))
      l = 1
      q = p
      t = 0.0_r8
      do while (l < p)
        s = a(q)
        m = p
        do while (m < q)
          s = s - a(l)*a(m)
          l = l + 1
          m = m + 1
        end do
        s = s / a(l)
        a(q) = s
        t = t + s**2
        l = l + 1
        q = q + 1
        if (q > size(a)) call ieee_set_flag(ieee_invalid, .true.)
      end do
      s = a(q) - t
      a(q) = sqrt(s)
      p = q + 1
    end do

  end subroutine upm_factor

  !! Solve the linear system A*x = b, where A is a symmetric matrix. The
  !! argument A holds the upper Choleski factor (as computed by UPM_FACTOR)
  !! in upper packed storage format. The RHS vector is provided in X, and
  !! is overwritten with the solution.

  pure subroutine upm_solve(a, x)

    real(r8), intent(in) :: a(:)
    real(r8), intent(inout) :: x(:)

    integer :: i, j, l
    real(r8) :: s

    if (size(x) == 0) return

    !! Forward substitution.
    x(1) = x(1)/a(1)
    l = 2
    do i = 2, size(x)
      s = x(i)
      do j = 1, i-1
        s = s - a(l)*x(j)
        l = l + 1
      end do
      x(i) = s/a(l)
      l = l + 1
    end do

    !! Backward substitution.
    l = l - 1
    do i = size(x), 2, -1
      x(i) = x(i)/a(l)
      l = l - 1
      do j = i-1, 1, -1
        x(j) = x(j) - a(l)*x(i)
        l = l - 1
      end do
    end do
    x(1) = x(1)/a(1)

  end subroutine upm_solve

  !! Computes the inverse of the symmetric, positive-definite matrix A stored
  !! in upper packed storage format. A is overwritten with its (symmetric)
  !! inverse, also in upper packed format. NB: To solve a system Ax = b with
  !! one or more RHS b, it is more efficient to use UPM_FACTOR and UPM_SOLVE,
  !! rather that computing the inverse of A and forming its product with b.

  pure subroutine upm_invert(a)

    real(r8), intent(inout) :: a(:)

    integer :: l, m, p, q
    real(r8) :: s

    if (size(a) == 0) return

    !! Compute the upper Choleski factor.
    call upm_factor(a)

    !! Compute the inverse of the Choleski factor.
    a(1) = 1.0_r8 / a(1)
    p = 2
    do while (p <= size(a))
      l = 1
      q = p
      do while (l < p)
        s = a(q)
        m = p
        do while (m < q)
          a(m) = a(m) + s*a(l)
          l = l + 1
          m = m + 1
        end do
        a(m) = s*a(l)
        l = l + 1
        q = q + 1
      end do
      s = 1.0_r8 / a(q)
      a(q) = s
      s = -s
      m = p
      do while (m < q)
        a(m) = s*a(m)
        m = m + 1
      end do
      p = q + 1
    end do

    !! Compute the product of the inverse Choleski factors.
    a(1) = a(1)*a(1)
    p = 2
    do while (p <= size(a))
      l = 1
      q = p
      do while (l < p)
        s = a(q)
        m = p
        do while (m <= q)
          a(l) = a(l) + s*a(m)
          l = l + 1
          m = m + 1
        end do
        q = q + 1
      end do
      s = a(q)
      do while (l <= q)
        a(l) = s*a(l)
        l = l + 1
      end do
      p = q + 1
    end do

  end subroutine upm_invert

  !! These functions are analogs of the intrinsic MATMUL function that return
  !! return the matrix-vector or matrix_matrix product when the first argument
  !! (and left factor) is a symmetric matrix stored in upper-packed storage
  !! format. The result is either a rank-1 or rank-2 array as appropriate.

  pure function sym_matmul_vec_f90(a, b) result(c)

    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: c(size(b))

    integer :: i, j, l
    real(r8) :: bi, ci

    l = 1
    do i = 1, size(b)
      bi = b(i)
      ci = 0.0_r8
      do j = 1, i-1
        ci = ci + a(l)*b(j)
        c(j) = c(j) + a(l) * bi
        l = l + 1
      end do
      c(i) = ci + a(l)*bi
      l = l + 1
    end do

  end function sym_matmul_vec_f90

  pure function sym_matmul_mat_f90(a, b) result(c)

    real(r8), intent(in) :: a(:), b(:,:)
    real(r8) :: c(size(b,1),size(b,2))

    integer :: i, j, k, l
    real(r8) :: bi, ci

    do k = 1, size(b,2)
      l = 1
      do i = 1, size(b,1)
        bi = b(i,k)
        ci = 0.0_r8
        do j = 1, i-1
          ci = ci + a(l)*b(j,k)
          c(j,k) = c(j,k) + a(l) * bi
          l = l + 1
        end do
        c(i,k) = ci + a(l)*bi
        l = l + 1
      end do
    end do

  end function sym_matmul_mat_f90

  pure function sym_matmul_vec_f77(n, a, b) result(c)

    integer, intent(in) :: n
    real(r8), intent(in) :: a(*), b(*)
    real(r8) :: c(n)

    integer :: i, j, l
    real(r8) :: bi, ci

    l = 1
    do i = 1, n
      bi = b(i)
      ci = 0.0_r8
      do j = 1, i-1
        ci = ci + a(l)*b(j)
        c(j) = c(j) + a(l) * bi
        l = l + 1
      end do
      c(i) = ci + a(l)*bi
      l = l + 1
    end do

  end function sym_matmul_vec_f77

  pure function sym_matmul_mat_f77(n, m, a, b) result(c)

    integer, intent(in) :: n, m
    real(r8), intent(in) :: a(*), b(n,*)
    real(r8) :: c(n,m)

    integer :: i, j, k, l
    real(r8) :: bi, ci

    do k = 1, m
      l = 1
      do i = 1, n
        bi = b(i,k)
        ci = 0.0_r8
        do j = 1, i-1
          ci = ci + a(l)*b(j,k)
          c(j,k) = c(j,k) + a(l) * bi
          l = l + 1
        end do
        c(i,k) = ci + a(l)*bi
        l = l + 1
      end do
    end do

  end function sym_matmul_mat_f77

  !! Compute the quadratic form x^T A x where A is a symmetric matrix
  !! in upper packed storage format.

  pure function upm_quad_form_f90(a, x) result(xtax)
    !use,intrinsic :: iso_fortran_env, only: r16 => real128
    real(r8), intent(in) :: a(:), x(:)
    real(r8) :: xtax
    integer :: i, j, l
    !real(r16) :: s ! quad precision accumulator
    real(r8) :: s
    s = 0.0_r8
    l = 1
    do i = 1, size(x)
      do j = 1, i-1
        s = s + 2*a(l)*x(i)*x(j)
        l = l + 1
      end do
      s = s + a(l)*x(i)*x(i)
      l = l + 1
    end do
    xtax = s
  end function

  pure function upm_quad_form_f77(n, a, x) result(xtax)
    !use,intrinsic :: iso_fortran_env, only: r16 => real128
    integer, intent(in) :: n
    real(r8), intent(in) :: a(*), x(*)
    real(r8) :: xtax
    integer :: i, j, l
    !real(r16) :: s ! quad precision accumulator
    real(r8) :: s
    s = 0.0_r8
    l = 1
    do i = 1, n
      do j = 1, i-1
        s = s + 2*a(l)*x(i)*x(j)
        l = l + 1
      end do
      s = s + a(l)*x(i)*x(i)
      l = l + 1
    end do
    xtax = s
  end function

  !! Compute the quadratic form x^H A x where A is a real symmetric matrix
  !! in upper packed storage format and x is a complex vector.

  pure function upm_zquad_form_f90(a, x) result(xtax)
    !use,intrinsic :: iso_fortran_env, only: r16 => real128
    real(r8), intent(in) :: a(:)
    complex(r8), intent(in) :: x(:)
    real(r8) :: xtax
    integer :: i, j, l
    !real(r16) :: s ! quad precision accumulator
    real(r8) :: s
    s = 0.0_r8
    l = 1
    do i = 1, size(x)
      do j = 1, i-1
        s = s + 2*a(l)*(x(i)%re*x(j)%re + x(i)%im*x(j)%im)
        l = l + 1
      end do
      s = s + a(l)*(x(i)%re*x(i)%re + x(i)%im*x(i)%im)
      l = l + 1
    end do
    xtax = s
  end function

  pure function upm_inner_prod_f90(a, x, y) result(xtay)
    real(r8), intent(in) :: a(:), x(:), y(:)
    real(r8) :: xtay
    integer :: i, j, l
    real(r8) :: s
    s = 0.0_r8
    l = 1
    do i = 1, size(x)
      do j = 1, i-1
        s = s + a(l)*(x(i)*y(j) + x(j)*y(i))
        l = l + 1
      end do
      s = s + a(l)*x(i)*y(i)
      l = l + 1
    end do
    xtay = s
  end function

  pure function upm_inner_prod_f77(n, a, x, y) result(xtay)
    integer, intent(in) :: n
    real(r8), intent(in) :: a(*), x(*), y(*)
    real(r8) :: xtay
    integer :: i, j, l
    real(r8) :: s
    s = 0.0_r8
    l = 1
    do i = 1, n
      do j = 1, i-1
        s = s + a(l)*(x(i)*y(j) + x(j)*y(i))
        l = l + 1
      end do
      s = s + a(l)*x(i)*y(i)
      l = l + 1
    end do
    xtay = s
  end function

  pure function upm_cong_prod_f90(a, b) result(c)
    real(r8), intent(in) :: a(:), b(:,:)
    real(r8) :: c((size(b,2)*(size(b,2)+1))/2)
    integer :: i, j, l
    l = 1
    do i = 1, size(b,2)
      do j = 1, i-1
        c(l) = upm_inner_prod(a, b(:,i), b(:,j))
        l = l + 1
      end do
      c(l) = upm_quad_form(a, b(:,i))
      l = l + 1
    end do
  end function

  pure function upm_cong_prod_f77(n, m, a, b) result(c)
    integer, intent(in) :: n, m
    real(r8), intent(in) :: a(*), b(n,*)
    real(r8) :: c((m*(m+1))/2)
    integer :: i, j, l
    l = 1
    do i = 1, m
      do j = 1, i-1
        c(l) = upm_inner_prod(n, a, b(:,i), b(:,j))
        l = l + 1
      end do
      c(l) = upm_quad_form(n, a, b(:,i))
      l = l + 1
    end do
  end function

  !! Compute the column sums of the elements of the symmetric matrix A stored
  !! in upper packed format.

  subroutine upm_col_sum(a, csum)

    real(r8), intent(in)  :: a(:)
    real(r8), intent(out) :: csum(:)

    integer :: i, j, l
    real(r8) :: s

    ASSERT(size(a) == (size(csum)*(size(csum)+1))/2)

    l = 1
    do i = 1, size(csum)
      s = 0.0_r8
      do j = 1, i-1
        s = s + a(l)
        csum(j) = csum(j) + a(l)
        l = l + 1
      end do
      csum(i) = s + a(l)
      l = l + 1
    end do

  end subroutine upm_col_sum

end module upper_packed_matrix_procs
