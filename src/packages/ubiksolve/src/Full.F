! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: Full.F,v 1.15 2003/10/19 22:12:48 turner Exp $
! 
!**********************************************************************
module Ubik_Full_module
  use ubik_types_module
  implicit none

  ! Intrinsics.

  interface MatVec
     module procedure MatVec_Full
  end interface

  ! Functions.

  interface Diag
     module procedure Diag_Full
  end interface

  interface Norm
     module procedure Norm_Full
  end interface

  interface Ubik_Sum
     module procedure Sum_rank2
  end interface

  ! Subroutines.

  interface Dump
     module procedure DumpMatrix_Full
     module procedure DumpSystem_Full
  end interface

  interface EnsureFormat
     module procedure EnsureFormat_Full
  end interface

  interface IC
     module procedure IC_Full
  end interface

  interface ILU
     module procedure ILU_Full
  end interface

  interface Jacobi
     module procedure Jacobi_Full
  end interface

  interface Load
     module procedure LoadMatrix_Full
     module procedure LoadSystem_Full
  end interface

  interface LUdecomp
     module procedure LUdecomp_Full
  end interface

  interface ScaleMatrixByDiagonal
     module procedure ScaleMatrixByDiagonal_Full
  end interface

  interface ScaleMatrixByRow
     module procedure ScaleMatrixByRow_Full
  end interface

  interface ScaleSystemByDiagonal
     module procedure ScaleSystemByDiagonal_Full
  end interface

  interface ScaleSystemByRow
     module procedure ScaleSystemByRow_Full
  end interface

  interface SolveLTriang
     module procedure SolveLTriang_Full
  end interface

  interface SolveUTriang
     module procedure SolveUTriang_Full
  end interface

  interface SOR
     module procedure SOR_Full
  end interface

  interface Spew
     module procedure WriteMatrixFloat_Full
     module procedure WriteMatrixInt_Full
     module procedure WriteSystem_Full
  end interface

  private
  public:: &
       Diag, &
       Norm, &
       Ubik_Sum, &
       Dump, &
       EnsureFormat, &
       IC, &
       ILU, &
       Jacobi, &
       Load, &
       LUdecomp, &
       MatVec, &
       ScaleMatrixByDiagonal, &
       ScaleMatrixByRow, &
       ScaleSystemByDiagonal, &
       ScaleSystemByRow, &
       SolveLTriang, &
       SolveUTriang, &
       SOR, &
       Spew

contains

  function MatVec_Full(a, x)

    ! <arguments>
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type), intent(in), dimension(:) :: x
    real(real_type), dimension(SIZE(a,dim=1)) :: MatVec_Full
    ! </arguments>

!!$    ! Local:
!!$    integer(int_type) :: i, j, m, n
!!$
!!$    m = SIZE(a,dim=1)
!!$    n = SIZE(a,dim=2)

    ! just use MatMul intrinsic
    MatVec_Full = MatMul(a, x)

!!$    ! use array syntax
!!$    MatVec_Full = 0.0
!!$    do j=1,n
!!$       MatVec_Full = MatVec_Full + a(:,j)*x(j)
!!$    end do
!!$
!!$    ! fully-expanded, ji ordering
!!$    MatVec_Full = 0.0
!!$    do j=1,n
!!$       do i=1,m
!!$          MatVec_Full(i) = MatVec_Full(i) + a(i,j)*x(j)
!!$       end do
!!$    end do
!!$
!!$    ! fully-expanded, ij ordering
!!$    MatVec_Full = 0.0
!!$    do i=1,m
!!$       do j=1,n
!!$          MatVec_Full(i) = MatVec_Full(i) + a(i,j)*x(j)
!!$       end do
!!$    end do

    return
  end function MatVec_Full

  !*********************************************************************
  !
  ! <para>
  ! Return diagonal elements of full-storage matrix as a vector.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In:
  !   a - matrix in full-storage format
  ! 
  !  Out:
  !   diag - vector of diagonal elements
  !
  !*********************************************************************
  function Diag_Full(a)

    ! <arguments>
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type), dimension(MIN( SIZE(a, dim=1), SIZE(a, dim=2) )) :: Diag_Full
    ! </arguments>

    ! Local:
    integer(int_type) :: i, n

    n = MIN( SIZE(a, dim=1), SIZE(a, dim=2) )
    do i=1,n
       Diag_Full(i) = a(i,i)
    end do

    return
  end function Diag_Full
  
  !*****************************************************************
  !
  ! <para>
  ! Compute the norm of a rectangular matrix in full-storage format.
  ! </para>
  !
  !                                     n
  !                                    ---
  !                                    \  
  !        infinity norm  ==>    max   /   | a  |
  !                            1<=i<=m ---    ij
  !                                    j=1
  !                              m 
  !                             ---
  !                             \  
  !        1-norm  ==>    max   /   | a  |
  !                     1<=j<=n ---    ij
  !                             i=1
  !                                   m   n
  !                                  --- ---
  !                                  \   \    2
  !        Frobenius norm  ==>  sqrt /   /   a
  !                                  --- ---  ij
  !                                  i=1 j=1
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   inorm - determines which vector norm to use
  !      0  ==>  infinity norm
  !      1  ==>  1-norm
  !      2  ==>  Frobenius norm
  !   a - matrix
  ! 
  !  Output:
  !   MatrixNorm_Full
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  MaxVal
  !  Sum
  !
  !*****************************************************************
  function Norm_Full(inorm, a)
    use Ubik_Vector_module, only: Ubik_MaxVal, Ubik_Sum

    ! <arguments>
    integer(int_type), intent(in) :: inorm
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type) :: Norm_Full
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, nrows, ncols
    real(real_type), dimension(SIZE(a, dim=1)) :: rowsum
    real(real_type), dimension(SIZE(a, dim=2)) :: colsum

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    select case (inorm)
    case (0)
       rowsum = 0.0
       do i=1,nrows
          rowsum(i) = Ubik_Sum( ABS( a(i,:) ) )
       end do
       Norm_Full = Ubik_MaxVal(rowsum)

    case (1)
       colsum = 0.0
       do j=1,ncols
          colsum(j) = Ubik_Sum( ABS( a(:,j) ) )
       end do
       Norm_Full = Ubik_MaxVal(colsum)

    case (2)
       Norm_Full = SQRT ( Ubik_Sum( a**2 ) )

    end select

    return
  end function Norm_Full
  
  !*********************************************************************
  !
  ! <para>
  ! Compute the sum of the elements of an array.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   x - array to be summed
  !
  !  Output:
  !   Sum_rank2
  !
  !*********************************************************************
  function Sum_rank2(x)
#ifdef use_PGSLib
    use Ubik_ParallelInfo_module, only: PGSLib_scope_flag
    use PGSLib_module
#endif

    ! <arguments>
    real(real_type), dimension(:,:), intent(in) :: x
    real(real_type) :: Sum_rank2
    ! </arguments>

#ifdef use_PGSLib
    Sum_rank2 = PGSLib_Global_SUM(x, scope=PGSLib_scope_flag)
#else
    Sum_rank2 = SUM(x)
#endif

    return
  end function Sum_rank2

  !*********************************************************************
  !
  ! <para>
  ! Dump the nonzero elements of a coefficient matrix in full storage format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit
  !   a - coefficient matrix
  !
  !*********************************************************************
  subroutine DumpMatrix_Full (lu, a)
    use Ubik_WriteString_module

    ! <arguments>
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    character(len=200) :: out_string
    integer(int_type) :: i, j, nrows, ncols

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    ! Write out sizes.
    write(out_string,*) '#', nrows, ncols
    call Spew (lu, out_string)

    ! Write out non-zero elements.
    do i=1,nrows
       do j=1,ncols
          if (ABS(a(i,j)) > TINY(1.0d0)) then
             write(out_string,100) i, j, a(i,j)
             call Spew (lu, out_string)
          end if
       end do
    end do

    return
100 format (2(i8, 1x), 1pg15.8)
  end subroutine DumpMatrix_Full

  !*********************************************************************
  !
  ! <para>
  ! Dump an entire linear system with coefficient in full-storage format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit to which to write
  !   a - coefficient matrix
  !   x - vector of unknowns
  !   b - vector containing right-hand side
  ! 
  ! <SUBROUTINES REQUIRED>
  ! 
  !  DumpCoeff
  !
  !*********************************************************************
  subroutine DumpSystem_Full (lu, a, x, b)
    use Ubik_Vector_module, only: Dump
    use Ubik_WriteString_module

    ! <arguments>
    real(real_type), intent(in), dimension(:,:) :: a

#include "DumpSystem-guts.F"

  end subroutine DumpSystem_Full

  !*********************************************************************
  !
  ! <para>
  ! Check a matrix in full-storage format.
  ! </para>
  !
  !*********************************************************************
  subroutine EnsureFormat_Full (a)

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a
    ! </arguments>

    return
  end subroutine EnsureFormat_Full
  
  !***************************************************************
  !
  ! <para>
  ! Overwrite a matrix in full-storage format by a level zero
  ! relaxed incomplete Cholesky factorization (fill elements
  ! are lumped into diagonal rather than simply being dropped).
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   omega - relaxation parameter
  ! 
  !  In/Out:
  !   a - matrix to be factored
  !
  !***************************************************************
  subroutine IC_Full (omega, a)

    ! <arguments>
    real(real_type), intent(in) :: omega
    real(real_type), intent(inout), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, k, n

    n = MIN(SIZE(a, dim=1), SIZE(a, dim=2))

    do k=1,n-1
       a(k,k) = SQRT(a(k,k))
       do i=k+1,n
          if (a(i,k) /= 0.0) then
             a(i,k) = a(i,k) / a(k,k)
             do j=i,n
                if (a(k,j) /= 0.0) then
                   if (a(i,j) /= 0.0) then
                      a(i,j) = a(i,j) - a(i,k)*a(k,j)
                   else
                      a(i,i) = a(i,i) - omega*a(i,k)*a(k,j)
                   end if
                end if
             end do
          end if
       end do
    end do
    a(n,n) = SQRT(a(n,n))

    return
  end subroutine IC_Full

  !***************************************************************
  !
  ! <para>
  ! Overwrite a matrix in full-storage format by a level zero
  ! relaxed incomplete LU-factorization (fill elements are
  ! lumped into diagonal rather than simply being dropped).
  ! 
  ! NOTE: No pivoting is performed.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   omega - relaxation parameter
  ! 
  !  In/Out:
  !   a - matrix to be factored
  !
  !***************************************************************
  subroutine ILU_Full (omega, a)

    ! <arguments>
    real(real_type), intent(in) :: omega
    real(real_type), intent(inout), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, k, n
    real(real_type) :: amult

    n = MIN(SIZE(a, dim=1), SIZE(a, dim=2))

    do k=1,n-1
       amult = 1.0/a(k,k)
       do i=k+1,n
          if (a(i,k) /= 0.0) then
             a(i,k) = amult*a(i,k)
             do j=k+1,n
                if (a(k,j) /= 0.0) then
                   if (a(i,j) /= 0.0) then
                      a(i,j) = a(i,j) - a(i,k)*a(k,j)
                   else
                      a(i,i) = a(i,i) - omega*a(i,k)*a(k,j)
                   end if
                end if
             end do
          end if
       end do
    end do

    return
  end subroutine ILU_Full

  !*******************************************************************
  !
  ! <para>
  ! Perform a single Jacobi iteration using a matrix in full-storage format.
  ! </para>
  !
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   omega - relaxation parameter
  !   b - source vector
  !   a - coefficient matrix
  ! 
  !  In/Out:
  !   x - solution vector
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  Dot_Product
  !
  !**********************************************************************
  subroutine Jacobi_Full (omega, b, a, x)
    ! <arguments>
    real(real_type), intent(in), dimension(:,:) :: a

#include "Jacobi-guts.F"

  end subroutine Jacobi_Full

  !*********************************************************************
  !
  ! <para>
  ! Load the nonzero elements of a matrix in full storage format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit from which to read
  ! 
  !  Output:
  !   a - coefficient matrix
  !
  !*********************************************************************
  subroutine LoadMatrix_Full (lu, a)

    ! <arguments>
    integer(int_type), intent(in) :: lu
    real(real_type), intent(out), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, ij, nrows, ncols
    real(real_type) :: element

    ! Very bad things will happen if these do not match what was dumped.
    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    ! Read non-zero elements.
    do ij=1,nrows*ncols
       read(lu,*,end=9999,err=9999) i, j, element
       if (i == 0) exit
       a(i,j) = element
    end do

9999 continue
    return
  end subroutine LoadMatrix_Full

  !*********************************************************************
  !
  ! <para>
  ! Load an entire linear system for a matrix in full-storage format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit from which to read
  !   a - coefficient matrix
  !   x - vector of unknowns
  !   b - vector containing right-hand side
  !
  !*********************************************************************
  subroutine LoadSystem_Full (lu, a, x, b)
    use Ubik_Vector_module, only: Load

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a

#include "LoadSystem-guts.F"

  end subroutine LoadSystem_Full
  
  !*****************************************************************
  !
  ! <para>
  ! Overwrite a square n by n matrix, a, stored in full format, by its 
  ! LU decomposition.  KJI ordering of the loops is used.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In/Out:
  !   a - matrix to be factored on entry, LU decomposition on exit
  !
  !*****************************************************************
  subroutine LUdecomp_Full (a)

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: j, k, n
    real(real_type) :: amult

    n = MIN(SIZE(a,dim=1),SIZE(a,dim=2))

    k_loop: do k=1,n-1
       amult = 1.0 / a(k,k)
       a(k+1:n,k) = amult*a(k+1:n,k)
       j_loop: do j=k+1,n
          a(k+1:n,j) = a(k+1:n,j) - a(k,j)*a(k+1:n,k)
       end do j_loop
    end do k_loop

    return
  end subroutine LUdecomp_Full

  !*********************************************************************
  !
  ! <para>
  ! Solve lower triangular system stored in full format.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   idiag - indicates whether diagonal is unity
  !     0 => general case
  !     non-zero => unit diagonal
  !   itrans - indicates whether to use transpose
  !     0 => no
  !     non-zero => use transpose
  !   a - matrix
  !
  !  In/Out:
  !   x - source vector on input, solution on output
  !
  !*********************************************************************
  subroutine SolveLTriang_Full (unit_diagonal, use_transpose, a, x)

    ! <arguments>
    logical(logical_type), intent(in) :: unit_diagonal, use_transpose
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, n

    n = MIN(SIZE(a,dim=1),SIZE(a,dim=2))

    if (use_transpose) then

       ! Use transpose.
       if (unit_diagonal) then

          ! Unit diagonal (UNUSED, UNTESTED).
          do i=2,n
             do j=1,i-1
                x(i) = x(i) - a(j,i)*x(j)
             end do
          end do
       else

          ! Non-unit diagonal (UNUSED, UNTESTED).
          x(1) = x(1) / a(1,1)
          do i=2,n
             do j=1,i-1
                x(i) = x(i) - a(j,i)*x(j)
             end do
             x(i) = x(i) / a(i,i)
          end do
       end if
    else

       ! Use coeff. as is.
       if (unit_diagonal) then

          ! Unit diagonal (used in ILU preconditioning).
          do j=1,n-1
             x(j+1:n) = x(j+1:n) - x(j)*a(j+1:n,j)
          end do
       else

          ! Non-unit diagonal (used in IC preconditioning).
          do j=1,n-1
             x(j) = x(j) / a(j,j)
             x(j+1:n) = x(j+1:n) - x(j)*a(j+1:n,j)
          end do
          x(n) = x(n) / a(n,n)
       end if
    end if

    return
  end subroutine SolveLTriang_Full

  !*********************************************************************
  !
  ! <para>
  ! Solve upper triangular system in full format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   idiag - indicates whether diagonal is unity
  !     0 => general case
  !     non-zero => unit diagonal
  !   itrans - indicates whether to use transpose
  !     0 => no
  !     non-zero => use transpose
  !   a - matrix
  ! 
  !  In/Out:
  !   x - source vector on input, solution on output
  !
  !*********************************************************************
  subroutine SolveUTriang_Full (unit_diagonal, use_transpose, a, x)

    ! <arguments>
    logical(logical_type), intent(in) :: unit_diagonal, use_transpose
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, n

    n = MIN(SIZE(a,dim=1),SIZE(a,dim=2))

    if (use_transpose) then

       ! Use transpose.
       if (unit_diagonal) then

          ! Unit diagonal (UNUSED, UNTESTED).
          do i=n-1,1,-1
             do j=i+1,n
                x(i) = x(i) - a(j,i)*x(j)
             end do
          end do
       else

          ! Non-unit diagonal (used in IC preconditioning).
          x(n) = x(n) / a(n,n)
          do i=n-1,1,-1
             do j=i+1,n
                x(i) = x(i) - a(j,i)*x(j)
             end do
             x(i) = x(i) / a(i,i)
          end do
       end if
    else

       ! Use coeff. as is.
       if (unit_diagonal) then

          ! Unit diagonal (UNUSED, UNTESTED).
          do j=n,2,-1
             x(1:j-1) = x(1:j-1) - x(j)*a(1:j-1,j)
          end do
       else

          ! Non-unit diagonal (used in ILU preconditioning and GMRES).
          do j=n,2,-1
             x(j) = x(j) / a(j,j)
             x(1:j-1) = x(1:j-1) - x(j)*a(1:j-1,j)
          end do
          x(1) = x(1) / a(1,1)
       end if
    end if

    return
  end subroutine SolveUTriang_Full
  
  !**********************************************************************
  !
  ! <para>
  ! Performs one iteration of SOR or SSOR for a system in full-storage
  ! format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   imeth - determines whether or not to perform SSOR
  !   omega - relaxation parameter
  !   b - source vector
  !   a - coefficient matrix
  ! 
  !  In/Out:
  !   x - solution vector
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  Dot_Product
  !
  !**********************************************************************
  subroutine SOR_Full (imeth, omega, b, a, x)
    use ubik_constants_module
    use Ubik_Vector_module, only: Ubik_Dot_Product

    ! <arguments>
    integer(int_type), intent(in) :: imeth
    real(real_type), intent(in) :: omega
    real(real_type), intent(in), dimension(:) :: b
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: i, n
    real(real_type) :: dot

    n = SIZE(x)

    do i=1,n
       dot = Ubik_Dot_Product(a(i,:), x) - a(i,i)*x(i)
       dot = ( b(i) - dot ) / a(i,i)
       x(i) = x(i) + omega*( dot - x(i) )
    end do

    if (imeth == method_SSOR) then
       do i=n,1,-1
          dot = Ubik_Dot_Product(a(i,:), x) - a(i,i)*x(i)
          dot = ( b(i) - dot ) / a(i,i)
          x(i) = x(i) + omega*( dot - x(i) )
       end do
    end if

    return
  end subroutine SOR_Full

  !*********************************************************************
  !
  ! <para>
  ! Write out rectangular floating-point matrix.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   lu - logical unit to which to write
  !   a - matrix
  !   title - character description
  !
  ! <SUBROUTINES REQUIRED>
  !
  !  WriteString
  !
  !*********************************************************************
  subroutine WriteMatrixFloat_Full (lu, a, title)
    use Ubik_WriteString_module

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    character(len=200) :: out_string
    integer(int_type) :: i, j, nrows, ncols

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    call Spew (lu, ' ')
    call Spew (lu, title)
    call Spew (lu, ' ')
    if (ncols <= 15) then
       do i=1,nrows
          write(out_string,500) i, (a(i,j),j=1,ncols)
          call Spew (lu, out_string)
       end do
    else
       do j=1,ncols
          call Spew (lu, ' ')
          write(out_string,*) 'Column: ', j
          call Spew (lu, out_string)
          call Spew (lu, ' ')
          do i=1,nrows
             write(out_string,500) i, a(i,j)
             call Spew (lu, out_string)
          end do
       end do
    end if

    return

500 format (t2, i7, 15(2x,f5.1))
  end subroutine WriteMatrixFloat_Full
  
  !*********************************************************************
  !
  ! <para>
  ! Write out rectangular integer matrix.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   lu - logical unit to which to write
  !   a - matrix
  !   title - character description
  !
  !*********************************************************************
  subroutine WriteMatrixInt_Full (lu, a, title)
    use Ubik_WriteString_module

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    integer(int_type), intent(in), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    character(len=132) :: out_string
    integer(int_type) :: i, j, nrows, ncols

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    call Spew (lu, ' ')
    call Spew (lu, title)
    call Spew (lu, ' ')
    if (ncols <= 10) then
       do i=1,nrows
          write(out_string,501) i, (a(i,j),j=1,ncols)
          call Spew (lu, out_string)
       end do
    else
       do j=1,ncols
          call Spew (lu, ' ')
          write(out_string,*) 'Column: ', j
          call Spew (lu, out_string)
          call Spew (lu, ' ')
          do i=1,nrows
             write(out_string,501) i, a(i,j)
             call Spew (lu, out_string)
          end do
       end do
    end if

    return

501 format (t2, i7, 2x, 10(2x,i7))
  end subroutine WriteMatrixInt_Full
  
  !*********************************************************************
  !
  ! <para>
  ! Write out entire linear system, Ax = b, for A in full-storage format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit to which to write
  !   a - coefficient matrix
  !   x - unknown vector
  !   b - source vector
  !   title - character description
  ! 
  ! <SUBROUTINES REQUIRED>
  ! 
  !  WriteString
  !
  !*********************************************************************
  subroutine WriteSystem_Full (lu, a, x, b, title)
    use Ubik_WriteString_module

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:,:) :: a
    real(real_type), intent(in), dimension(:) :: b, x
    ! </arguments>

    ! Local:
    character(len=200) :: out_string
    integer(int_type) :: i, j, nrows, ncols

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    call Spew (lu, ' ')
    call Spew (lu, title)
    call Spew (lu, ' ')
    if (ncols <= 7) then
       do i=1,nrows
          write(out_string,500) i, (a(i,j),j=1,ncols), x(i), b(i)
          call Spew (lu, out_string)
       end do
    else
       do j=1,ncols
          call Spew (lu, ' ')
          write(out_string,*) 'Coefficient, Column: ', j
          call Spew (lu, out_string)
          call Spew (lu, ' ')
          do i=1,nrows
             write(out_string,500) i, a(i,j)
             call Spew (lu, out_string)
          end do
       end do
       call Spew (lu, ' ')
       out_string = 'Unknown and Source vectors'
       call Spew (lu, out_string)
       call Spew (lu, ' ')
       do i=1,nrows
          write(out_string,500) i, x(i), b(i)
          call Spew (lu, out_string)
       end do
    end if

    return

500 format (t2, i7, 1p, 9(2x,e13.6))
  end subroutine WriteSystem_Full

  !*********************************************************************
  !
  ! <para>
  ! Row-scale a matrix stored in full format.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleMatrixByRow_Full (a)

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, ncols
    real(real_type), allocatable, dimension(:) :: max_val_row

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    allocate (max_val_row(nrows))
    max_val_row = MaxVal(a, dim=2)

    a = a / SPREAD(max_val_row, dim=2, ncopies=ncols)

    deallocate (max_val_row)
    return
  end subroutine ScaleMatrixByRow_Full

  !*********************************************************************
  !
  ! <para>
  ! Row-scale a system with the coefficient stored in full format.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleSystemByRow_Full (a, b)

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a
    real(real_type), intent(inout), dimension(:) :: b
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, ncols
    real(real_type), allocatable, dimension(:) :: max_val_row

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    allocate (max_val_row(nrows))
    max_val_row = MaxVal(a, dim=2)

    a = a / SPREAD(max_val_row, dim=2, ncopies=ncols)
    b = b / max_val_row

    deallocate (max_val_row)
    return
  end subroutine ScaleSystemByRow_Full

  !*********************************************************************
  !
  ! <para>
  ! Scale a matrix stored in full format by the diagonal elements.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleMatrixByDiagonal_Full (a)

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, ncols
    real(real_type), allocatable, dimension(:) :: diagonal

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    allocate (diagonal(nrows))
    diagonal = Diag(a)

    a = a / SPREAD(diagonal, dim=2, ncopies=ncols)

    deallocate (diagonal)
    return
  end subroutine ScaleMatrixByDiagonal_Full

  !*********************************************************************
  !
  ! <para>
  ! Scale a system with coefficient stored in full format by the diagonal
  ! elements.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleSystemByDiagonal_Full (a, b)

    ! <arguments>
    real(real_type), intent(inout), dimension(:,:) :: a
    real(real_type), intent(inout), dimension(:) :: b
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, ncols
    real(real_type), allocatable, dimension(:) :: diagonal

    nrows = SIZE(a, dim=1)
    ncols = SIZE(a, dim=2)

    allocate (diagonal(nrows))
    diagonal = Diag(a)

    a = a / SPREAD(diagonal, dim=2, ncopies=ncols)
    b = b / diagonal

    deallocate (diagonal)
    return
  end subroutine ScaleSystemByDiagonal_Full

end module Ubik_Full_module
