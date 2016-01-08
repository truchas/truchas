!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE LU_SOLVE_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Solution of equations A*x=b and/or the factorization of A using
  !   column-wise LU-factorization with or without partial pivoting.
  !
  !   LU-Factorization: A*x=L*U*x=b
  !     Upper Triangular Matrix (U) - Full Matrix
  !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
  !   Back Substitution: L*y=b and U*x=y
  !
  ! Contains: LU_SOLVE
  !
  !           LUINV
  !           LUINV_FACTORIZE
  !           LUINV_BACKSOLVE
  !
  !           LUINP
  !           LUINP_FACTORIZE
  !           LUINP_BACKSOLVE
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  ! Public Subroutines
  public :: LU_SOLVE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Factorization/Solution Parameters
  integer, parameter, public :: factorize = 1 ! Factorization Only
  integer, parameter, public :: backsolve = 2 ! Back Substitution Only
  integer, parameter, public :: factsolve = 3 ! Factorization
                                                               ! and Back Substitution

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE LU_SOLVE (A, B, Indx, solve, pivot)
    !=======================================================================
    ! Purpose(s):
    !   Solution of equations A*x=b and/or the factorization of A using
    !   column-wise LU-factorization with or without partial pivoting.
    !   Solution returned in vector b.
    !
    !   Note: The objective of these routines was to have a LU solution
    !   capability with an 'LU re-use' option, the repeated use of a
    !   factorized matrix, for both the pivot and non-pivot options. In
    !   order to also make these routines simple, if the pivot option is
    !   selected the user is required to supply a pivot indexing array.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !
    !   solve = factorize : LU-Factorization Only
    !   solve = backsolve : Back Substitution Only
    !   solve = factsolve : LU-Factorization and Back Substitution
    !
    !   pivot = .true.  : LU-Factorization with Partial Pivoting
    !   pivot = .false. : LU-Factorization without Partial Pivoting
    !=======================================================================

    ! Argument List
    logical, intent(IN), optional :: pivot
    integer, intent(IN), optional :: solve

    real(r8), dimension(:,:), intent(INOUT) :: A
    real(r8), dimension(:),   intent(INOUT) :: B
    integer,  dimension(:),   intent(INOUT), optional :: Indx

    ! Local Variables
    logical :: pivot_flag
    logical :: pivot_default = .false.

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Pivot Flag
    if (PRESENT(pivot)) then
       pivot_flag = pivot         ! Input Pivot
    else
       pivot_flag = pivot_default ! Default Pivot
    end if

    ! Pivot Index Array
    INSIST(PRESENT(Indx) .or. .not.pivot_flag)

    ! Factorize/Solve Equations
    select case (pivot_flag)
    case (.false.) ! Non-Pivot Solution
       call LUINV (A, B,       solve)
    case (.true.)  ! Pivot Solution
       call LUINP (A, B, Indx, solve)
    end select

  END SUBROUTINE LU_SOLVE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LUINV (A, B, solve)
    !=======================================================================
    ! Purpose(s):
    !   Solution of equations A*x=b and/or the factorization of A using
    !   column-wise LU-factorization. Solution returned in vector b.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !
    !   solve = factorize : LU-Factorization Only
    !   solve = backsolve : Back Substitution Only
    !   solve = factsolve : LU-Factorization and Back Substitution
    !=======================================================================

    ! Argument List
    integer, intent(IN), optional :: solve

    real(r8), dimension(:,:), intent(INOUT) :: A
    real(r8), dimension(:),   intent(INOUT) :: B

    ! Local Variables
    integer :: solve_flag
    integer :: solve_default = factsolve

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Factorize/Solve Flag
    if (PRESENT(solve)) then
       solve_flag = solve         ! Input Flag
    else
       solve_flag = solve_default ! Default Flag
    end if

    ! Factorize/Solve Equations
    select case (solve_flag)
    case (factorize) ! LU-Factorization Only
       call LUINV_FACTORIZE (A)
    case (backsolve) ! Back Substitution Only
       call LUINV_BACKSOLVE (A, B)
    case (factsolve) ! LU-Factorization and Back Substitution
       call LUINV_FACTORIZE (A)
       call LUINV_BACKSOLVE (A, B)
    end select

  END SUBROUTINE LUINV

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LUINV_FACTORIZE (A)
    !=======================================================================
    ! Purpose(s):
    !   Factorization of an n*n matrix using column-wise LU-factorization.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !=======================================================================

    ! Argument List
    real(r8), dimension(:,:), intent(INOUT) :: A

    ! Local Variables
    integer :: n
    integer :: i, j, k

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Matrix/Vector Size
    n = SIZE(A,1)

    ! LU-Factorization: A*x=L*U*x=b
    do j = 1, n
       ! Upper Triangular Matrix - Full Matrix
       do i = 1, j
          do k = 1, i-1
             A(i,j) = A(i,j) - A(i,k)*A(k,j)
          end do
       end do

       ! Lower Triangular Matrix - Full Matrix with Identity Diagonal
       do i = j+1, n
          do k = 1, j-1
             A(i,j) = A(i,j) - A(i,k)*A(k,j)
          end do
          A(i,j) = A(i,j)/A(j,j)
       end do
    end do

  END SUBROUTINE LUINV_FACTORIZE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LUINV_BACKSOLVE (A, B)
    !=======================================================================
    ! Purpose(s):
    !   Solution of equations A*x=b where A is column-wise LU-factorized.
    !   Solution returned in vector b.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !=======================================================================

    ! Argument List
    real(r8), dimension(:,:), intent(IN)    :: A
    real(r8), dimension(:),   intent(INOUT) :: B

    ! Local Variables
    integer :: n
    integer :: i, j

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Matrix/Vector Size
    n = SIZE(A,1)

    ! Back Substitution: L*y=b and U*x=y
    ! Lower Triangular Matrix
    do i = 2, n
       do j = 1, i-1
          B(i) = B(i) - A(i,j)*B(j)
       end do
    end do

    ! Upper Triangular Matrix
    B(n) = B(n)/A(n,n)
    do i = n-1, 1, -1
       do j = n, i+1, -1
          B(i) = B(i) - A(i,j)*B(j)
       end do
       B(i) = B(i)/A(i,i)
    end do

  END SUBROUTINE LUINV_BACKSOLVE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LUINP (A, B, Indx, solve)
    !=======================================================================
    ! Purpose(s):
    !   Solution of equations A*x=b and/or the factorization of A using
    !   column-wise LU-factorization with partial pivoting. Solution
    !   returned in vector b.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !
    !   solve = factorize : LU-Factorization Only
    !   solve = backsolve : Back Substitution Only
    !   solve = factsolve : LU-Factorization and Back Substitution
    !=======================================================================

    ! Argument List
    integer, intent(IN), optional :: solve

    real(r8), dimension(:,:), intent(INOUT) :: A
    real(r8), dimension(:),   intent(INOUT) :: B
    integer,  dimension(:),   intent(INOUT) :: Indx

    ! Local Variables
    integer :: solve_flag
    integer :: solve_default = factsolve

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Factorize/Solve Flag
    if (PRESENT(solve)) then
       solve_flag = solve         ! Input Flag
    else
       solve_flag = solve_default ! Default Flag
    end if

    ! Factorize/Solve Equations
    select case (solve_flag)
    case (factorize) ! LU-Factorization Only
       call LUINP_FACTORIZE (A, Indx)
    case (backsolve) ! Back Substitution Only
       call LUINP_BACKSOLVE (A, B, Indx)
    case (factsolve) ! LU-Factorization and Back Substitution
       call LUINP_FACTORIZE (A, Indx)
       call LUINP_BACKSOLVE (A, B, Indx)
    end select

  END SUBROUTINE LUINP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LUINP_FACTORIZE (A, Indx)
    !=======================================================================
    ! Purpose(s):
    !   Factorization of an n*n matrix using column-wise LU-factorization
    !   with partial pivoting.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !=======================================================================
    ! Argument List
    real(r8), dimension(:,:), intent(INOUT) :: A
    integer, dimension(SIZE(A,1)), intent(INOUT) :: Indx

    ! Local Variables
    integer :: n
    integer :: imax
    integer :: i, j, k

    real(r8) :: dum
    real(r8) :: smax
    real(r8), dimension(SIZE(A,1)) :: T

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Matrix/Vector Size
    n = SIZE(A,1)

    ! LU-Factorization: A*x=L*U*x=b
    do i = 1, n
       T(i) = 0.0_r8
       do j = 1, n
          if (ABS(A(i,j)) > T(i)) T(i) = ABS(A(i,j))
       end do
       T(i) = 1.0_r8/T(i)
    end do

    do j = 1, n
       do i = 1, j-1
          do k = 1, i-1
             A(i,j) = A(i,j) - A(i,k)*A(k,j)
          end do
       end do

       imax = j
       smax = 0.0_r8
       do i = j, n
          do k = 1, j-1
             A(i,j) = A(i,j) - A(i,k)*A(k,j)
          end do
          if (T(i)*ABS(A(i,j)) >= smax) then
             imax = i
             smax = T(i)*ABS(A(i,j))
          end if
       end do

       if (imax /= j) then
          do k = 1, n
             dum = A(imax,k)
             A(imax,k) = A(j,k)
             A(j,k) = dum
          end do
          T(imax) = T(j)
       end if
       Indx(j) = imax

       do i = j+1, n
          A(i,j) = A(i,j)/A(j,j)
       end do
    end do

  END SUBROUTINE LUINP_FACTORIZE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LUINP_BACKSOLVE (A, B, Indx)
    !=======================================================================
    ! Purpose(s):
    !   Solution of equations A*x=b where A is column-wise LU-factorized
    !   with partial pivoting. Solution returned in vector b.
    !
    !   LU-Factorization: A*x=L*U*x=b
    !     Upper Triangular Matrix (U) - Full Matrix
    !     Lower Triangular Matrix (L) - Full Matrix with Identity Diagonal
    !   Back Substitution: L*y=b and U*x=y
    !=======================================================================

    ! Argument List
    real(r8), dimension(:,:), intent(IN)    :: A
    real(r8), dimension(:),   intent(INOUT) :: B
    integer, dimension(SIZE(A,1)), intent(IN) :: Indx

    ! Local Variables
    integer :: n
    integer :: i, j
    real(r8), dimension(SIZE(A,1)) :: T

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Matrix/Vector Size
    n = SIZE(A,1)

    ! Back Substitution: L*y=b and U*x=y
    ! Lower Triangular Matrix
    do i = 1, n
       T(i) = B(Indx(i))
       B(Indx(i)) = B(i)
       do j = 1, i-1
          T(i) = T(i) - A(i,j)*B(j)
       end do
       B(i) = T(i)
    end do

    ! Upper Triangular Matrix
    B(n) = B(n)/A(n,n)
    do i = n-1, 1, -1
       do j = n, i+1, -1
          B(i) = B(i) - A(i,j)*B(j)
       end do
       B(i) = B(i)/A(i,i)
    end do

  END SUBROUTINE LUINP_BACKSOLVE

END MODULE LU_SOLVE_MODULE
