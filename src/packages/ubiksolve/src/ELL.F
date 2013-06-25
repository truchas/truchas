! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: ELL.F,v 1.18 2003/10/19 22:12:48 turner Exp $
! 
!*********************************************************************
module Ubik_ELL_module
  use ubik_types_module
  implicit none

  ! Type definition.

  type ELL_matrix
     real(real_type), dimension(:,:), pointer :: values => NULL()
     integer(int_type), dimension(:,:), pointer :: map => NULL()
     integer(int_type) :: status
  end type ELL_matrix

  ! Operators.

  interface assignment(=)
     module procedure ELL_eq_int_type
     module procedure ELL_eq_real_type
     module procedure ELL_eq_real
  end interface

  ! Intrinsics.

  interface MatVec
     module procedure MatVec_ELL
  end interface

  interface SIZE
     module procedure SIZE_ELL
  end interface

  ! Functions.

  interface Diag
     module procedure Diag_ELL
  end interface

  interface DiagonalDominance
     module procedure DiagonalDominance_ELL
  end interface

  interface Norm
     module procedure Norm_ELL
  end interface

  interface NumberOfNonpositiveRows
     module procedure NumberOfNonpositiveRows_ELL
  end interface

  interface NumberOfNonsymmElements
     module procedure NumberOfNonsymmElements_ELL
  end interface

  ! Subroutines.

  interface create
     module procedure create_and_size
     module procedure create_and_size_and_init_int
     module procedure create_and_size_and_init_float
  end interface

  interface destroy
     module procedure destroy_ELL
  end interface

  interface set_status_alloc_failure
     module procedure set_status_alloc_failure_ELL
  end interface

  interface alloc_failure
     module procedure alloc_failure_ELL
  end interface

  interface Dump
     module procedure DumpMatrix_ELL
     module procedure DumpSystem_ELL
  end interface

  interface EnsureFormat
     module procedure EnsureFormat_ELL
  end interface

  interface IC
     module procedure IC_ELL
  end interface

  interface ILU
     module procedure ILU_ELL
  end interface

  interface Jacobi
     module procedure Jacobi_ELL
  end interface

  interface LeftJustifyMatrix
     module procedure LeftJustifyMatrix_ELL
  end interface

  interface Load
     module procedure LoadMatrix_ELL
     module procedure LoadSystem_ELL
  end interface

  interface LUdecomp
     module procedure LUdecomp_ELL
  end interface

  interface ScaleMatrixByDiagonal
     module procedure ScaleMatrixByDiagonal_ELL
  end interface

  interface ScaleMatrixByRow
     module procedure ScaleMatrixByRow_ELL
  end interface

  interface ScaleSystemByDiagonal
     module procedure ScaleSystemByDiagonal_ELL
  end interface

  interface ScaleSystemByRow
     module procedure ScaleSystemByRow_ELL
  end interface

  interface SolveLTriang
     module procedure SolveLTriang_ELL
  end interface

  interface SolveUTriang
     module procedure SolveUTriang_ELL
  end interface

  interface SOR
     module procedure SOR_ELL
  end interface

  interface Spew
     module procedure WriteMatrix_ELL
     module procedure WriteSystem_ELL
  end interface

  private
  public :: &
       ELL_matrix, &
       assignment(=), &
       create, &
       destroy, &
       set_status_alloc_failure, &
       alloc_failure, &
       MatVec, &
       SIZE, &
       Diag, &
       DiagonalDominance, &
       Norm, &
       NumberOfNonpositiveRows, &
       NumberOfNonsymmElements, &
       Dump, &
       EnsureFormat, &
       IC, &
       ILU, &
       Jacobi, &
       LeftJustifyMatrix, &
       Load, &
       LUdecomp, &
       ScaleMatrixByDiagonal, &
       ScaleMatrixByRow, &
       ScaleSystemByDiagonal, &
       ScaleSystemByRow, &
       SolveLTriang, &
       SolveUTriang, &
       SOR, &
       Spew

contains

  !*****************************************************************
  !
  ! Set routines.
  ! 
  !*****************************************************************
  ! <para>
  ! Set ELL matrix elements to a scalar value given an integer.
  ! </para>
  subroutine ELL_eq_int_type (a, i)

    ! <arguments>
    integer(int_type), intent(in) :: i
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    a%values = DBLE(i)
    a%map = i

    return
  end subroutine ELL_eq_int_type

  ! <para>
  ! Set ELL matrix elements to a scalar value given a double-precision real.
  ! </para>
  subroutine ELL_eq_real_type (a, s)

    ! <arguments>
    real(real_type), intent(in) :: s
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    a%values = s
    a%map = INT(s)

    return
  end subroutine ELL_eq_real_type

  ! <para>
  ! Set ELL matrix elements to a scalar value given a single-precision real.
  ! </para>
  subroutine ELL_eq_real (a, r)
    real, intent(in) :: r
    type(ELL_matrix), intent(inout) :: a

    ! <arguments>
    a%values = DBLE(r)
    a%map = INT(r)
    ! </arguments>

    return
  end subroutine ELL_eq_real

  !*******************************************************************
  !
  ! constructors and destructors
  !
  !**********************************************************************
  ! <para>
  ! Allocate an ELL-format matrix.
  ! </para>
  subroutine create_and_size (a, nrows, maxnz)

    ! <arguments>
    integer(int_type), intent(in) :: maxnz, nrows
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: status

    ! if resizing, need to deallocate first - note that must check
    ! association status first, since SIZE requires a status of
    ! associated
    if (ASSOCIATED(a%values) .and. ASSOCIATED(a%map)) then
       if (SIZE(a%values, dim=1) /= nrows .or. &
           SIZE(a%values, dim=2) /= maxnz) then
          call destroy (a)
          if (alloc_failure(a)) then
             return
          end if
       end if
    end if

    ! allocate
    if (.not.ASSOCIATED(a%values) .and. .not.ASSOCIATED(a%map)) then
       allocate (a%values(nrows,maxnz), a%map(nrows,maxnz), stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(a)
       end if
    end if

    return
  end subroutine create_and_size

  ! <para>
  ! Allocate an ELL-format matrix and initialize to a scalar value given an integer.
  ! </para>
  subroutine create_and_size_and_init_int (a, nrows, maxnz, i)
    ! <arguments>
    integer(int_type), intent(in) :: maxnz, nrows, i
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    call create (a, nrows, maxnz)
    if (alloc_failure(a)) then
       return
    end if

    a = i

    return
  end subroutine create_and_size_and_init_int

  ! <para>
  ! Allocate an ELL-format matrix and initialize to a scalar value given a
  ! floating-point number.
  ! </para>
  subroutine create_and_size_and_init_float (a, nrows, maxnz, x)
    ! <arguments>
    integer(int_type), intent(in) :: maxnz, nrows
    real(real_type), intent(in) :: x
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    call create (a, nrows, maxnz)
    if (alloc_failure(a)) then
       return
    end if

    a = x

    return
  end subroutine create_and_size_and_init_float

  ! <para>
  ! Safely deallocate an ELL-format matrix.
  ! </para>
  subroutine destroy_ELL (a)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>
    integer(int_type) :: status

    if (ASSOCIATED(a%values)) then
       deallocate(a%values, stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(a)
       end if
    end if

    if (ASSOCIATED(a%map)) then
       deallocate(a%map, stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(a)
       end if
    end if

    return
  end subroutine destroy_ELL

  !*****************************************************************
  !
  ! status setting and accessing
  !
  !*****************************************************************

  ! <para>
  ! Set this ELL object's status flag to indicate allocation or
  ! deallocation failure.
  ! </para>
  subroutine set_status_alloc_failure_ELL(a)
    use ubik_constants_module

    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    a%status = alloc_failure__
    return
  end subroutine set_status_alloc_failure_ELL

  ! <para>
  ! Return this ELL object's status flag.
  ! </para>
  function alloc_failure_ELL(a)
    use ubik_constants_module

    ! <arguments>
    logical(logical_type) :: alloc_failure_ELL
    type(ELL_matrix), intent(in) :: a
    ! </arguments>

    alloc_failure_ELL = (a%status == alloc_failure__)
    return
  end function alloc_failure_ELL

  !*****************************************************************
  !
  ! <para>
  ! Compute y=Ax where A is a matrix stored in ELL format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   a - ELL matrix
  !   x - x-vector
  ! 
  !  Output:
  !   y - y-vector
  !
  !*****************************************************************
  function MatVec_ELL(a, x)

    ! <arguments>
    type(ELL_matrix), intent(in) :: a
    real(real_type), intent(in), dimension(:) :: x
    real(real_type), dimension(SIZE(a%values, dim=1)) :: MatVec_ELL
    ! </arguments>

    ! Local:
    integer(int_type) :: j
    real(real_type), dimension(0:SIZE(x)) :: xtmp

    xtmp(1:) = x
    xtmp(0) = 0.0

    MatVec_ELL = 0.0
    do j=1,SIZE(a%values, dim=2)
!       where (a%map(:,j) /= 0) y_eq_Ax = y_eq_Ax + a%values(:,j)*x(a%map(:,j))
       MatVec_ELL = MatVec_ELL + a%values(:,j)*xtmp(a%map(:,j))
    end do

    return
  end function MatVec_ELL

  !*****************************************************************
  !
  ! <para>
  ! Return the size of an ELL-format matrix.
  ! </para>
  !
  !*****************************************************************
  function SIZE_ELL(a)
    ! <arguments>
    type(ELL_matrix), intent(in) :: a
    integer(int_type) :: SIZE_ELL
    ! </arguments>

    SIZE_ELL = SIZE(a%values)

    return
  end function SIZE_ELL

  !*****************************************************************
  !
  ! <para>
  ! Return diagonal elements of a matrix stored in ELL format as a vector.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   a - ELL matrix
  ! 
  !  Output:
  !   diag - vector of diagonal elements
  !
  !*****************************************************************
  function Diag_ELL(a)
    ! <arguments>
    type(ELL_matrix), intent(in) :: a
    real(real_type), dimension(SIZE(a%values, dim=1)) :: Diag_ELL
    ! </arguments>

    Diag_ELL = a%values(:,1)

    return
  end function Diag_ELL

  !*****************************************************************
  !
  ! <para>
  ! Compute diagonal dominance of a matrix stored in ELL format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   a - ELL matrix
  ! 
  !  Output:
  !   DiagonalDominance
  !
  !*****************************************************************
  function DiagonalDominance_ELL (a, dd)

    ! <arguments>
    type(ELL_matrix), intent(in) :: a
    real(real_type), dimension(:), intent(out) :: dd
    real(real_type) :: DiagonalDominance_ELL
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, maxnz
    real(real_type), dimension(:), allocatable :: off_diag

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

!!$    if (.not.allocated(dd)) then
!!$       allocate (dd(nrows))
!!$    end if

    allocate (off_diag(nrows))
    off_diag = SUM(ABS(a%values(:,2:maxnz)), dim=2)
    where(off_diag > TINY(1.0))
       dd = a%values(:,1) / off_diag
    end where

    DiagonalDominance_ELL = MinVal(dd, MASK=(dd > TINY(1.0)))

    return
  end function DiagonalDominance_ELL

  !*****************************************************************
  !
  ! <para>
  ! Compute positivity (number of non-positive rows) of a matrix stored
  ! in ELL format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   a - ELL matrix
  ! 
  !  Output:
  !   DiagonalDominance
  !
  !*****************************************************************
  function NumberOfNonpositiveRows_ELL (a)

    ! <arguments>
    type(ELL_matrix), intent(in) :: a
    integer(int_type) :: NumberOfNonpositiveRows_ELL
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, nrows, maxnz
    real(real_type) :: diagonal

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    NumberOfNonpositiveRows_ELL = 0
    do i = 1, nrows
       diagonal = a%values(i,1)
       if (diagonal > 0.0) then
          do j = 2, maxnz
             if (a%values(i,j) > 0.0) then
                NumberOfNonpositiveRows_ELL = NumberOfNonpositiveRows_ELL + 1
             end if
          end do
       else
          do j = 2, maxnz
             if (a%values(i,j) < 0.0) then
                NumberOfNonpositiveRows_ELL = NumberOfNonpositiveRows_ELL + 1
             end if
          end do
       end if
    end do

    return
  end function NumberOfNonpositiveRows_ELL

  !*****************************************************************
  !
  ! <para>
  ! Check symmetry of a matrix stored in ELL format.  Number of nonsymmetric
  ! elements is returned.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   a - ELL matrix
  ! 
  !  Output:
  !   DiagonalDominance
  !
  !*****************************************************************
  function NumberOfNonsymmElements_ELL (a, maxdiff)

    ! <arguments>
    type(ELL_matrix), intent(in) :: a
    real(real_type), intent(out) :: maxdiff
    integer(int_type) :: NumberOfNonsymmElements_ELL
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, maxnz
    integer(int_type) :: i, j, jj, j1, j2
    real(real_type) :: diff

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    maxdiff = 0.0
    NumberOfNonsymmElements_ELL = 0
    i1_loop: &
    do i = 1, nrows
      j1_loop: &
      do j = 2, maxnz
        j1 = a%map(i,j)
        if (j1 == 0) then
           cycle i1_loop
        end if
        j2_loop: &
        do jj = 2, maxnz
           j2 = a%map(j1,jj)
           if (j2 == 0) then
              cycle j1_loop
           end if
           if (j2 == i) then
              diff = ABS(a%values(i,j) - a%values(j1,jj))
              if (diff > TINY(1.0)) then
                 NumberOfNonsymmElements_ELL = NumberOfNonsymmElements_ELL + 1
                 maxdiff = MAX(diff, maxdiff)
!!$                 write(3,*) " "
!!$                 write(3,*) i, j, j1, a%values(i,j)
!!$                 write(3,*) j1, jj, j2, a%values(j1,jj)
              end if
              cycle j1_loop
           end if
        end do j2_loop
        NumberOfNonsymmElements_ELL = NumberOfNonsymmElements_ELL + 1
        maxdiff = MAX(ABS(a%values(i,j)), maxdiff)
!!$        write(3,*) " "
!!$        write(3,*) i, j, j1, a%values(i,j)
!!$        write(3,*) j1, i, " missing"
      end do j1_loop
    end do i1_loop

    return
  end function NumberOfNonsymmElements_ELL
  
  !*****************************************************************
  !
  ! <para>
  ! Compute the norm of a matrix stored in ELL format.
  ! </para>
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
  !      1  ==>  1-norm (NOT IMPLEMENTED)
  !      2  ==>  Frobenius norm
  !   a - matrix
  ! 
  !  Output:
  !   Norm - norm of matrix A
  ! 
  !  <FUNCTIONS REQUIRED>
  ! 
  !   MaxVal
  !   Sum
  !
  !*********************************************************************
  function Norm_ELL(inorm, a)
    use Ubik_Full_module, only: Ubik_Sum
    use Ubik_Vector_module, only: Ubik_MaxVal, Ubik_Sum

    ! <arguments>
    integer(int_type), intent(in) :: inorm
    type(ELL_matrix), intent(in) :: a
    real(real_type) :: Norm_ELL
    ! </arguments>

    ! Local:
    integer(int_type) :: i, nrows, ncols
    real(real_type), dimension(SIZE(a%values, dim=1)) :: rowsum

    nrows = SIZE(a%values, dim=1)
!    ncols = SIZE(a%values, dim=2)

    select case (inorm)
    case (0)
       rowsum = 0.0
       do i=1,nrows
          rowsum(i) = Ubik_Sum( ABS( a%values(i,:) ) )
       end do
       Norm_ELL = Ubik_MaxVal(rowsum)

    case (1)
       ! Need a chunk of code here to do 1-norm.

    case (2)
       Norm_ELL = SQRT ( Ubik_Sum( a%values**2 ) )

    end select

    return
  end function Norm_ELL

  !*********************************************************************
  !
  ! <para>
  ! Dump the nonzero elements of a matrix stored in ELL format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit
  !   a - coefficient matrix
  !
  !*********************************************************************
  subroutine DumpMatrix_ELL (lu, a)
    use Ubik_WriteString_module

    ! <arguments>
    integer(int_type), intent(in) :: lu
    type(ELL_matrix), intent(in) :: a
    ! </arguments>

    ! Local:
    character(len=200) :: out_string
    integer(int_type) :: i, j, nrows, maxnz

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    ! Write out sizes.
    write(out_string,*) '#', nrows, maxnz
    call Spew (lu, out_string)

    ! Write out non-zero elements.
    do i=1,nrows
       do j=1,maxnz
          if (ABS(a%values(i,j)) > TINY(1.0d0)) then
             write(out_string,100) i, a%map(i,j), a%values(i,j), j
             call Spew (lu, out_string)
          end if
       end do
    end do

    return
100 format (2(i8, 1x), 1pg15.8, 1x, i8)
  end subroutine DumpMatrix_ELL

  !*********************************************************************
  !
  ! <para>
  ! Dump an entire linear system with coefficient stored in ELL format.
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
  !*********************************************************************
  subroutine DumpSystem_ELL (lu, a, x, b)
    use Ubik_Vector_module, only: Dump
    use Ubik_WriteString_module

    ! <arguments>
    type(ELL_matrix), intent(in) :: a

#include "DumpSystem-guts.F"

  end subroutine DumpSystem_ELL
  
  !***************************************************************
  !
  ! <para>
  ! Overwrite a matrix in ELL storage format by a level zero
  ! relaxed incomplete Cholesky factorization.
  ! 
  ! Note that the j loop goes from i+1 to n rather than from 1 to
  ! maxnz (which it could do in a standard IC routine).  This
  ! increases cost significantly for large systems.  Is there a
  ! better way?
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
  !   a - original matrix on input, factored matrix on output
  !
  !***************************************************************
  subroutine IC_ELL (omega, a)

    ! <arguments>
    real(real_type), intent(in) :: omega
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, ij, j, k, kj, kk, n, maxnz

    n = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    k_loop: do k=1,n-1
       a%values(k,1) = SQRT(a%values(k,1))

       ! Loop over rows k+1 to n.
       i_loop: do i=k+1,n

          ! Find kk such that a%map(i,kk) = k.  Note that if one is found,
          ! we also know that a%values(i,kk) is non-zero.
          kk_loop: do kk=2,maxnz
             if (a%map(i,kk) == 0) then

                ! Out of elements, so a%values(i,kk) = 0, so bail out and get another i.
                cycle i_loop
             else if (a%map(i,kk) == k) then

                ! Found kk.
                go to 15
             end if
          end do kk_loop

          ! Didn't find kk, so a%values(i,kk) = 0, so bail out and get another i.
          cycle i_loop

15        continue
          a%values(i,kk) = a%values(i,kk) / a%values(k,1)

          ! Loop over columns k+1 to n.
          j_loop: do j=k+1,n

             ! Find kj such that a%map(k,kj) = j.  Note that if one is found,
             ! we also know that a%values(k,kj) is non-zero.
             kj_loop: do kj=2,maxnz
                if (a%map(k,kj) == 0) then

                   ! Out of elements, so a%values(k,kj) = 0, so bail out and get
                   ! another j.
                   cycle j_loop
                else if (a%map(k,kj) == j) then

                   ! Found kj.
                   go to 20
                end if
             end do kj_loop

             ! Didn't find kj, so a%values(k,kj) = 0, so bail out and get another j.
             cycle j_loop

             ! Find ij such that a%map(i,ij) = j.  Whether one is found or not
             ! will determine whether to modify the diagonal or the off-diagonal.
20           continue
             ij_loop: do ij=2,maxnz
                if (a%map(i,ij) == 0) then

                   ! Out of elements, so a%values(i,ij) = 0, so bail out and modify
                   ! diagonal.
                   exit ij_loop
                else if (a%map(i,ij) == j) then

                   ! Found ij, so modify element a%values(i,ij).
                   a%values(i,ij) = a%values(i,ij) - a%values(i,kk)*a%values(k,kj)
                   cycle j_loop
                end if
             end do ij_loop

             ! Didn't find ij, so a%values(i,ij) = 0, so modify diagonal.
             a%values(i,1) = a%values(i,1) - omega*a%values(i,kk)*a%values(k,kj)
          end do j_loop
       end do i_loop
    end do k_loop
    a%values(n,1) = SQRT(a%values(n,1))

    return
  end subroutine IC_ELL
  
  !***************************************************************
  !
  ! <para>
  ! Overwrite a matrix in ELL storage format by a level zero
  ! relaxed incomplete LU-factorization (fill elements are
  ! lumped into diagonal rather than simply being dropped).
  ! 
  ! Note that the main difference between this routine and
  ! ILU_ELL is that the j loop goes from i+1 to n rather
  ! than from 1 to maxnz.  This increases cost significantly
  ! for large systems.  Is there a better way?
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
  subroutine ILU_ELL (omega, a)

    ! <arguments>
    real(real_type), intent(in) :: omega
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, ij, j, k, kj, kk, n, maxnz
    real(real_type) :: amult

    n = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    k_loop: do k=1,n-1
       amult = 1.0/a%values(k,1)

       ! Loop over rows k+1 to n.
       i_loop: do i=k+1,n

          ! Find kk such that a%map(i,kk) = k.  Note that if one is found,
          ! we also know that a%values(i,kk) is non-zero.
          kk_loop: do kk=2,maxnz
             if (a%map(i,kk) == 0) then

                ! Out of elements, so a%values(i,kk) = 0, so bail out and get another i.
                cycle i_loop
             else if (a%map(i,kk) == k) then

                ! Found kk.
                go to 15
             end if
          end do kk_loop

          ! Didn't find kk, so a%values(i,kk) = 0, so bail out and get another i.
          cycle i_loop

15        continue
          a%values(i,kk) = amult*a%values(i,kk)

          ! Loop over columns k+1 to n.
          j_loop: do j=k+1,n

             ! Find kj such that a%map(k,kj) = j.  Note that if one is found,
             ! we also know that a%values(k,kj) is non-zero.
             kj_loop: do kj=2,maxnz
                if (a%map(k,kj) == 0) then

                   ! Out of elements, so a%values(k,kj) = 0, so bail out and get
                   ! another j.
                   cycle j_loop
                else if (a%map(k,kj) == j) then

                   ! Found kj.
                   go to 20
                end if
             end do kj_loop

             ! Didn't find kj, so a%values(k,kj) = 0, so bail out and get another j.
             cycle j_loop

             ! Find ij such that a%map(i,ij) = j.  Whether one is found or not
             ! will determine whether to modify the diagonal or the off-diagonal.
20           continue
             ij_loop: do ij=2,maxnz
                if (a%map(i,ij) == 0) then

                   ! Out of elements, so a%values(i,ij) = 0, so bail out and modify
                   ! diagonal.
                   exit ij_loop
                else if (a%map(i,ij) == j) then

                   ! Found ij, so modify element a%values(i,ij).
                   a%values(i,ij) = a%values(i,ij) - a%values(i,kk)*a%values(k,kj)
                   cycle j_loop
                end if
             end do ij_loop

             ! Didn't find ij, so a%values(i,ij) = 0, so modify diagonal.
             a%values(i,1) = a%values(i,1) - omega*a%values(i,kk)*a%values(k,kj)
          end do j_loop
       end do i_loop
    end do k_loop

    return
  end subroutine ILU_ELL

  !*******************************************************************
  !
  ! <para>
  ! Perform a single Jacobi iteration using a matrix stored in ELL format.
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
  subroutine Jacobi_ELL (omega, b, a, x)
    ! <arguments>
    type(ELL_matrix), intent(in) :: a

#include "Jacobi-guts.F"

  end subroutine Jacobi_ELL

  !*********************************************************************
  !
  ! <para>
  ! Left justifies non-zero elements in an ELL matrix.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In/Out:
  !   a - coefficient array
  !
  !*********************************************************************
  subroutine LeftJustifyMatrix_ELL (a)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, icol, j, n, maxnz, status
    integer(int_type), allocatable, dimension(:) :: jtmp
    real(real_type), allocatable, dimension(:) :: tmp

    n = SIZE(a%values,dim=1)
    maxnz = SIZE(a%values,dim=2)

    ! Allocate local arrays.
    allocate (jtmp(maxnz), tmp(maxnz), stat=status)
    if (status /= 0) then
       return
    end if

    do i=1,n
       icol = 1
       tmp = 0.0
       jtmp = 0
       do j=1,maxnz
          if (a%map(i,j) /= 0) then
             tmp(icol) = a%values(i,j)
             jtmp(icol) = a%map(i,j)
             icol = icol + 1
          end if
       end do

       a%values(i,:) = tmp
       a%map(i,:) = jtmp
    end do

    ! Deallocate local arrays.
    deallocate (jtmp, tmp)

    return
  end subroutine LeftJustifyMatrix_ELL

  !*********************************************************************
  !
  ! <para>
  ! Check a matrix stored in ELL format.
  ! </para>
  !
  !*********************************************************************
  subroutine EnsureFormat_ELL (a)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, jtmp, nrows, maxnz
    real(real_type) :: tmp

    nrows = SIZE(a%values,dim=1)
    maxnz = SIZE(a%values,dim=2)

    rows: &
    do i=1,nrows
       if (a%map(i,1) /= i) then

          ! diagonal element is not in column 1
          columns: &
          do j=2,maxnz
             if (a%map(i,j) == i) then

                ! found the diagonal element - swap it into column 1
                tmp = a%values(i,j)
                a%values(i,j) = a%values(i,1)
                a%values(i,1) = tmp

                ! don't forget the map array
                jtmp = a%map(i,j)
                a%map(i,j) = a%map(i,1)
                a%map(i,1) = jtmp
                cycle rows
             end if
          end do columns
       end if
    end do rows

    ! also left-justify
    call LeftJustifyMatrix (a)

    return
  end subroutine EnsureFormat_ELL

  !*********************************************************************
  !
  ! <para>
  ! Load the nonzero elements of a matrix stored in ELL format.
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
  subroutine LoadMatrix_ELL (lu, a)
    ! <arguments>
    integer(int_type), intent(in) :: lu
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: i, ij, j, jj, nrows, maxnz
    real(real_type) :: element

    ! Very bad things will happen if these do not match what was dumped.
    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    ! Read non-zero elements.
    do ij=1,nrows*maxnz
       read(lu,*,end=9999,err=9999) i, j, element, jj
       if (i == 0) then
          exit
       end if
       a%values(i,jj) = element
       a%map(i,jj) = j
    end do

9999 continue

    return
  end subroutine LoadMatrix_ELL

  !*********************************************************************
  !
  ! <para>
  ! Load an entire linear system for a matrix stored in ELL format.
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
  subroutine LoadSystem_ELL (lu, a, x, b)
    use Ubik_Vector_module, only: Load

    ! <arguments>
    type(ELL_matrix), intent(inout) :: a

#include "LoadSystem-guts.F"

  end subroutine LoadSystem_ELL
  
  !*****************************************************************
  !
  ! <para>
  ! Overwrite a square n by n matrix, a, stored in ELL format, by its 
  ! LU decomposition.  NOT IMPLEMENTED.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In/Out:
  !   a - matrix to be factored on entry, LU decomposition on exit
  !
  !*****************************************************************
  subroutine LUdecomp_ELL (a)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! NOT IMPLEMENTED

    return
  end subroutine LUdecomp_ELL

  !*********************************************************************
  !
  ! <para>
  ! Row-scale a matrix stored in ELL format.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleMatrixByRow_ELL (a)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, maxnz
    real(real_type), allocatable, dimension(:) :: max_val_row

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    allocate (max_val_row(nrows))
    max_val_row = MaxVal(a%values, dim=2)

    a%values = a%values / SPREAD(max_val_row, dim=2, ncopies=maxnz)

    deallocate (max_val_row)
    return
  end subroutine ScaleMatrixByRow_ELL

  !*********************************************************************
  !
  ! <para>
  ! Row-scale a system with the coefficient stored in ELL format.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleSystemByRow_ELL (a, b)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    real(real_type), dimension(:), intent(inout) :: b
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, maxnz
    real(real_type), allocatable, dimension(:) :: max_val_row

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    allocate (max_val_row(nrows))
    max_val_row = MaxVal(a%values, dim=2)

    a%values = a%values / SPREAD(max_val_row, dim=2, ncopies=maxnz)
    b = b / max_val_row

    deallocate (max_val_row)
    return
  end subroutine ScaleSystemByRow_ELL

  !*********************************************************************
  !
  ! <para>
  ! Scale a matrix stored in ELL format by the diagonal elements.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleMatrixByDiagonal_ELL (a)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, maxnz
    real(real_type), allocatable, dimension(:) :: diagonal

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    allocate (diagonal(nrows))
    diagonal = Diag(a)

    a%values(:,1) = 1.0
    a%values(:,2:maxnz) = a%values(:,2:maxnz) / SPREAD(diagonal, dim=2, ncopies=(maxnz-1))

    deallocate (diagonal)
    return
  end subroutine ScaleMatrixByDiagonal_ELL

  !*********************************************************************
  !
  ! <para>
  ! Scale a system with coefficient stored in ELL format by the diagonal
  ! elements.
  ! </para>
  !
  !*********************************************************************
  subroutine ScaleSystemByDiagonal_ELL (a, b)
    ! <arguments>
    type(ELL_matrix), intent(inout) :: a
    real(real_type), dimension(:), intent(inout) :: b
    ! </arguments>

    ! Local:
    integer(int_type) :: nrows, maxnz
    real(real_type), allocatable, dimension(:) :: diagonal

    nrows = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    allocate (diagonal(nrows))
    diagonal = Diag(a)

    a%values(:,1) = 1.0
    a%values(:,2:maxnz) = a%values(:,2:maxnz) / SPREAD(diagonal, dim=2, ncopies=(maxnz-1))
    b = b / diagonal

    deallocate (diagonal)
    return
  end subroutine ScaleSystemByDiagonal_ELL

  !*******************************************************************
  !
  ! <para>
  ! Solve unit lower triangular system stored in ELL format.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   unit_diagonal - indicates whether diagonal is unity
  !     .FALSE. => general case
  !     .TRUE. => unit diagonal
  !   use_transpose - indicates whether to use transpose
  !     .FALSE. => no
  !     .TRUE. => use transpose
  !   a - matrix in ELL format
  !
  !  In/Out:
  !   x - source vector on input, solution on output
  !
  !***************************************************************
  subroutine SolveLTriang_ELL (unit_diagonal, use_transpose, a, x)
    ! <arguments>
    logical(logical_type), intent(in) :: unit_diagonal, use_transpose
    type(ELL_matrix), intent(in) :: a
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, ji, n, maxnz

    n = SIZE(a%values,dim=1)
    maxnz = SIZE(a%values,dim=2)

    ! Forward elimination.
    if (use_transpose) then

       ! Use transpose.
       if (unit_diagonal) then

          ! Unit diagonal (UNUSED, UNTESTED).
          do i=2,n
             j_loop: do j=1,i-1

                ! Find ji such that a%map(j,ji) = i.  If one is found, we
                ! also know that a%values(j,ji) is non-zero.
                ji_loop: do ji=2,maxnz
                   if (a%map(j,ji) == 0) then

                      ! Out of elements, so a%values(j,ji) = 0, so bail out
                      ! and get another j.
                      cycle j_loop
                   else if (a%map(j,ji) == i) then

                      ! Found ji.
                      exit ji_loop
                   end if
                end do ji_loop
                x(i) = x(i) - a%values(j,ji)*x(j)
             end do j_loop
          end do
       else

          ! Non-unit diagonal (UNUSED, UNTESTED).
          x(1) = x(1) / a%values(1,1)
          do i=2,n
             j_loop2: do j=1,i-1

                ! Find ji such that a%map(j,ji) = i.  If one is found, we
                ! also know that a%values(j,ji) is non-zero.
                ji_loop2: do ji=2,maxnz
                   if (a%map(j,ji) == 0) then

                      ! Out of elements, so a%values(j,ji) = 0, so bail out
                      ! and get another j.
                      cycle j_loop2
                   else if (a%map(j,ji) == i) then

                      ! Found ji.
                      exit ji_loop2
                   end if
                end do ji_loop2
                x(i) = x(i) - a%values(j,ji)*x(j)
             end do j_loop2
             x(i) = x(i) / a%values(i,1)
          end do
       end if
    else

       ! Use coeff. as is.
       if (unit_diagonal) then

          ! Unit diagonal (used in ILU preconditioning).
          do i=2,n
             do j=2,maxnz
                if (a%map(i,j) == 0) then
                   exit
                else if (a%map(i,j) < i) then
                   x(i) = x(i) - a%values(i,j)*x(a%map(i,j))
                end if
             end do
          end do
       else

          ! Non-unit diagonal (used in IC preconditioning).
          x(1) = x(1) / a%values(1,1)
          do i=2,n
             do j=2,maxnz
                if (a%map(i,j) == 0) then
                   exit
                else if (a%map(i,j) < i) then
                   x(i) = x(i) - a%values(i,j)*x(a%map(i,j))
                end if
             end do
             x(i) = x(i) / a%values(i,1)
          end do
       end if
    end if

    return
  end subroutine SolveLTriang_ELL

  !*********************************************************************
  !
  ! <para>
  ! Solve upper triangular system stored in ELL format.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   unit_diagonal - indicates whether diagonal is unity
  !     .FALSE. => general case
  !     .TRUE. => unit diagonal
  !   use_transpose - indicates whether to use transpose
  !     .FALSE. => no
  !     .TRUE. => use transpose
  !   a - matrix in ELL format
  ! 
  !  In/Out:
  !   x - source vector on input, solution on output
  !
  !***************************************************************
  subroutine SolveUTriang_ELL (unit_diagonal, use_transpose, a, x)
    ! <arguments>
    logical(logical_type), intent(in) :: unit_diagonal, use_transpose
    type(ELL_matrix), intent(in) :: a
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: i, j, ji, n, maxnz

    n = SIZE(a%values,dim=1)
    maxnz = SIZE(a%values,dim=2)

    ! Back-substitution.
    if (use_transpose) then

       ! Use transpose.
       if (unit_diagonal) then

          ! Unit diagonal (UNUSED, UNTESTED).
          i_loop: do i=n-1,1,-1
             j_loop: do j=i+1,n

                ! Find ji such that a%map(j,ji) = i.  If one is found, we also
                ! know that a%values(j,ji) is non-zero.
                ji_loop: do ji=2,maxnz
                   if (a%map(j,ji) == 0) then

                      ! Out of elements, so a%values(j,ji) = 0, so bail out and get
                      ! another j.
                      cycle j_loop
                   else if (a%map(j,ji) == i) then

                      ! Found ji, modify i-th element of x.
                      x(i) = x(i) - a%values(j,ji)*x(j)
                      cycle j_loop
                   end if
                end do ji_loop
             end do j_loop
          end do i_loop
       else

          ! Non-unit diagonal (used in IC preconditioning).
          x(n) = x(n)/a%values(n,1)
          i_loop2: do i=n-1,1,-1
             j_loop2: do j=i+1,n

                ! Find ji such that a%map(j,ji) = i.  If one is found, we also
                ! know that a%values(j,ji) is non-zero.
                ji_loop2: do ji=2,maxnz
                   if (a%map(j,ji) == 0) then

                      ! Out of elements, so a%values(j,ji) = 0, so bail out and get
                      ! another j.
                      cycle j_loop2
                   else if (a%map(j,ji) == i) then

                      ! Found ji, modify i-th element of x.
                      x(i) = x(i) - a%values(j,ji)*x(j)
                      cycle j_loop2
                   end if
                end do ji_loop2
             end do j_loop2

             ! Finished with i-th element of x.
             x(i) = x(i) / a%values(i,1)
          end do i_loop2
       end if
    else

       ! Use coeff. as is.
       if (unit_diagonal) then

          ! Unit diagonal (UNUSED, UNTESTED).
          do i=n-1,1,-1
             do j=2,maxnz
                if (a%map(i,j) == 0) then
                   exit
                else if (a%map(i,j) > i) then
                   x(i) = x(i) - a%values(i,j)*x(a%map(i,j))
                end if
             end do
          end do
       else

          ! Non-unit diagonal (used in ILU preconditioning).
          x(n) = x(n)/a%values(n,1)
          do i=n-1,1,-1
             do j=2,maxnz
                if (a%map(i,j) == 0) then
                   exit
                else if (a%map(i,j) > i) then
                   x(i) = x(i) - a%values(i,j)*x(a%map(i,j))
                end if
             end do
             x(i) = x(i)/a%values(i,1)
          end do
       end if
    end if

    return
  end subroutine SolveUTriang_ELL
  
  !*******************************************************************
  !
  ! <para>
  ! Performs one iteration of SOR or SSOR for a system in ELL format.
  ! 
  ! Note that the main diagonal of the coefficient is assumed to be in
  ! the first column.
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
  subroutine SOR_ELL (imeth, omega, b, a, x)
    use ubik_constants_module
    use Ubik_Vector_module, only: Ubik_Dot_Product

    ! <arguments>
    integer(int_type), intent(in) :: imeth
    real(real_type), intent(in) :: omega
    real(real_type), intent(in), dimension(:) :: b
    type(ELL_matrix), intent(in) :: a
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    integer(int_type) :: i, n, maxnz, status
    real(real_type) :: dot
!    real(real_type), dimension(SIZE(a,dim=2)-1) :: row, col
    real(real_type), dimension(0:SIZE(x)) :: xtmp

    n = SIZE(x)
    maxnz = SIZE(a%values, dim=2)

    xtmp(1:) = x
    xtmp(0) = 0.0

    do i=1,n
!       where (a%map(i,2:maxnz) /= 0)
!          row = a%values(i,2:maxnz)
!          col = x(a%map(i,2:maxnz))
!       elsewhere
!          row = 0.0
!          col = 0.0
!       end where
!!       row = MERGE(a%values(i,2:maxnz), 0.0, a%map(i,2:maxnz) /= 0)
!!       col = MERGE(x(a%map(i,2:maxnz)), 0.0, a%map(i,2:maxnz) /= 0)
!       dot = Dot_Product( row, col )
       dot = Ubik_Dot_Product( a%values(i,2:maxnz), xtmp(a%map(i,2:maxnz)) )
       dot = ( b(i) - dot ) / a%values(i,1)
!       x(i) = x(i) + omega*( dot - x(i) )
       xtmp(i) = xtmp(i) + omega*( dot - xtmp(i) )
    end do

    if (imeth == method_SSOR) then
       do i=n,1,-1
!          where (a%map(i,2:maxnz) /= 0)
!             row = a%values(i,2:maxnz)
!             col = x(a%map(i,2:maxnz))
!          elsewhere
!             row = 0.0
!             col = 0.0
!          end where
!!          row = MERGE(a%values(i,2:maxnz), 0.0, a%map(i,2:maxnz) /= 0)
!!          col = MERGE(x(a%map(i,2:maxnz)), 0.0, a%map(i,2:maxnz) /= 0)
!          dot = Dot_Product( row, col )
          dot = Ubik_Dot_Product( a%values(i,2:maxnz), xtmp(a%map(i,2:maxnz)) )
          dot = ( b(i) - dot ) / a%values(i,1)
!          x(i) = x(i) + omega*( dot - x(i) )
          xtmp(i) = xtmp(i) + omega*( dot - xtmp(i) )
       end do
    end if

    x = xtmp(1:)

    return
  end subroutine SOR_ELL
  
  !*********************************************************************
  !
  ! <para>
  ! Write out ELL matrix.
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
  !  Spew
  !
  !*********************************************************************
  subroutine WriteMatrix_ELL (lu, a, title)
    use Ubik_WriteString_module
    use Ubik_Full_module, only: Spew

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    type(ELL_matrix), intent(in) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: status

    call Spew (lu, a%values, title)
    call Spew (lu, a%map, 'Column map for the above matrix:')

    return
  end subroutine WriteMatrix_ELL

  !*********************************************************************
  !
  ! <para>
  ! Write out entire linear system, Ax = b, for A stored in ELL format.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   lu - logical unit to which to write
  !   a - matrix
  !   x - unknown vector
  !   b - source vector
  !   title - character description
  !
  ! <SUBROUTINES REQUIRED>
  !
  !  Spew
  !
  !*********************************************************************
  subroutine WriteSystem_ELL (lu, a, x, b, title)
    use Ubik_Full_module, only: Spew

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:) :: b, x
    type(ELL_matrix), intent(in) :: a
    ! </arguments>

    ! Local:
    integer(int_type) :: status

    call Spew (lu, a%values, x, b, title)
    call Spew (lu, a%map, 'Column map for the above matrix:')

    return
  end subroutine WriteSystem_ELL

end module Ubik_ELL_module
