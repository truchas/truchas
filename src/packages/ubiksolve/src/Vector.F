! -*- Mode: f90 -*-

!*********************************************************************
! $Id: Vector.F,v 1.10 2002/12/03 04:55:05 turner Exp $
! 
! Vector routines.
!
!*********************************************************************
module Ubik_Vector_module
  use ubik_types_module
  implicit none

  ! Functions.

  interface Ubik_Dot_Product
     module procedure DotProduct
  end interface

  interface Ubik_MaxVal
     module procedure MaxVal_Vector
  end interface

  interface Norm
     module procedure Norm_Vector
  end interface

  interface Ubik_Sum
     module procedure Sum_Vector
  end interface

  ! Subroutines.

  interface Dump
     module procedure Dump_Vector
  end interface

  interface Load
     module procedure Load_Vector
  end interface

  interface Swap
     module procedure Swap_Vector
  end interface

  interface Spew
     module procedure Write1Vector
     module procedure Write2Vectors
     module procedure Write3Vectors
  end interface

  private
  public :: &
       Ubik_Dot_Product, &
       Ubik_MaxVal, &
       Norm, &
       Ubik_Sum, &
       Dump, &
       Load, &
       Swap, &
       Spew

contains
  
  !*********************************************************************
  !
  ! <para>
  ! Compute the dot product of two vectors.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   x, y - vectors to be dotted
  ! 
  !  Output:
  !   Dot_Product - dot product
  !
  !*********************************************************************
  function DotProduct(x, y)
#ifdef use_PGSLib
    use Ubik_ParallelInfo_module, only: PGSLib_scope_flag
    use PGSLib_module
#endif

    ! <arguments>
    real(real_type), dimension(:), intent(in) :: x, y
    real(real_type) :: DotProduct
    ! </arguments>

#ifdef use_PGSLib
    DotProduct = PGSLib_Global_Dot_Product(x, y, scope=PGSLib_scope_flag)
#else
    DotProduct = Dot_Product(x, y)
#endif

    return
  end function DotProduct

  !*********************************************************************
  !
  ! <para>
  ! Return the maximum value in a vector.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   x - array to be scanned
  !
  !  Output:
  !   MaxValVector
  !
  !*********************************************************************
  function MaxVal_Vector(x)
#ifdef use_PGSLib
    use Ubik_ParallelInfo_module, only: PGSLib_scope_flag
    use PGSLib_module
#endif

    ! <arguments>
    real(real_type), dimension(:), intent(in) :: x
    real(real_type) :: MaxVal_Vector
    ! </arguments>

#ifdef use_PGSLib
    MaxVal_Vector = PGSLib_Global_MAXVAL(x, scope=PGSLib_scope_flag)
#else
    MaxVal_Vector = MAXVAL(x)
#endif

    return
  end function MaxVal_Vector
  
  !*********************************************************************
  !
  ! <para>
  ! Compute the norm of a vector.
  !                                    
  !        infinity norm  ==>   max   | a |
  !                           1<=i<=n    i
  !                      n
  !                     ---
  !                     \  
  !        1-norm  ==>  /   | a |
  !                     ---    i
  !                     i=1
  !                           n
  !                          ---
  !                          \    2
  !        2-norm  ==>  sqrt /   a
  !                          ---  i
  !                          i=1
  ! </para>
  ! 
  ! <PARAMETER LIST>
  !  Input:
  !   inorm - determines which vector norm to use
  !      0  ==>  infinity norm (max absolute value)
  !      1  ==>  1-norm (sum of absolute values)
  !      2  ==>  2-norm, or Euclidean norm (square root of sum of
  !              squares)
  !   x - vector
  ! 
  !  Output:
  !   Norm - norm of x
  ! 
  ! <FUNCTIONS REQUIRED>
  ! 
  !  MaxVal
  !  Sum
  !
  !*********************************************************************
  function Norm_Vector(inorm, x)

    ! <arguments>
    integer(int_type), intent(in) :: inorm
    real(real_type), intent(in), dimension(:) :: x
    real(real_type) :: Norm_Vector
    ! </arguments>

    Norm_Vector = 0.0
    select case (inorm)
    case (0)
       Norm_Vector = Ubik_MaxVal( ABS(x) )
    case (1)
       Norm_Vector = Ubik_Sum( ABS(x) )
    case (2)
       Norm_Vector = SQRT( Ubik_Sum( x**2 ) )
    end select

    return
  end function Norm_Vector
  
  !*********************************************************************
  !
  ! <para>
  ! Compute the sum of the elements of a vector.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   x - array to be summed
  !
  !  Output:
  !   SumVector
  !
  !*********************************************************************
  function Sum_Vector(x)
#ifdef use_PGSLib
    use Ubik_ParallelInfo_module, only: PGSLib_scope_flag
    use PGSLib_module
#endif

    ! <arguments>
    real(real_type), dimension(:), intent(in) :: x
    real(real_type) :: Sum_Vector
    ! </arguments>

#ifdef use_PGSLib
    Sum_Vector = PGSLib_Global_SUM(x, scope=PGSLib_scope_flag)
#else
    Sum_Vector = SUM(x)
#endif

    return
  end function Sum_Vector

  !*********************************************************************
  !
  ! <para>
  ! Dump one, two, or three vectors.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit from which to read
  !
  !  Output:
  !   x, y, z - vectors
  !
  !*********************************************************************
  subroutine Dump_Vector (lu, x, y, z)
    use Ubik_WriteString_module

    ! <arguments>
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:) :: x
    real(real_type), intent(in), dimension(:), optional :: y, z
    ! </arguments>

    ! Local:
    character(len=200) :: out_string
    integer(int_type) :: i, n

    if (PRESENT(z)) then
       n = MIN( SIZE(x), SIZE(y), SIZE(z) )
       do i=1,n
          write(out_string,100) x(i), y(i), z(i)
          call Spew (lu, out_string)
       end do
    else if (PRESENT(y)) then
       n = MIN( SIZE(x), SIZE(y) )
       do i=1,n
          write(out_string,100) x(i), y(i)
          call Spew (lu, out_string)
       end do
    else
       n = SIZE(x)
       do i=1,n
          write(out_string,100) x(i)
          call Spew (lu, out_string)
       end do
    end if

9999 continue
    return
100 format (1p, 3(1x, g15.8))
  end subroutine Dump_Vector

  !*********************************************************************
  !
  ! <para>
  ! Load one, two, or three vectors.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit from which to read
  !
  !  Output:
  !   x, y, z - vectors
  !
  !*********************************************************************
  subroutine Load_Vector (lu, x, y, z)

    ! <arguments>
    integer(int_type), intent(in) :: lu
    real(real_type), intent(out), dimension(:) :: x
    real(real_type), intent(out), dimension(:), optional :: y, z
    ! </arguments>

    ! Local:
    integer(int_type) :: i, n

    if (PRESENT(z)) then
       n = MIN( SIZE(x), SIZE(y), SIZE(z) )
       do i=1,n
          read(lu,*,end=9999,err=9999) x(i), y(i), z(i)
       end do
    else if (PRESENT(y)) then
       n = MIN( SIZE(x), SIZE(y) )
       do i=1,n
          read(lu,*,end=9999,err=9999) x(i), y(i)
       end do
    else
       n = SIZE(x)
       do i=1,n
          read(lu,*,end=9999,err=9999) x(i)
       end do
    end if

9999 continue
    return
  end subroutine Load_Vector

  !*********************************************************************
  !
  ! <para>
  ! Swap the elements of two vectors.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In/Out:
  !   x - vector x
  !   y - vector y
  !
  !*********************************************************************
  subroutine Swap_Vector (x, y)

    ! <arguments>
    real(real_type), intent(inout), dimension(:) :: x, y
    ! </arguments>

    ! Local:
    real(real_type), dimension(SIZE(x)) :: save

    save = y
    y = x
    x = save

    return
  end subroutine Swap_Vector
  
  !*********************************************************************
  !
  ! <para>
  ! Write out vector.
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   lu - logical unit to which to write
  !   n - number of elements to write out
  !   x - vector
  !   title - character description
  !
  ! <SUBROUTINES REQUIRED>
  !
  !  WriteString
  !
  !*********************************************************************
  subroutine Write1Vector (lu, x, title)
    use Ubik_WriteString_module

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:) :: x
    ! </arguments>

    ! Local:
    character(len=132) :: out_string
    integer(int_type) :: i, status

    call Spew (lu, ' ')
    call Spew (lu, title)
    call Spew (lu, ' ')
    do i=1,SIZE(x)
       write(out_string,500) i, x(i)
       call Spew (lu, out_string)
    end do

500 format (t2, i7, 1p, 3x, e16.7)
    return
  end subroutine Write1Vector

  ! <para>
  ! Write out two vectors.
  ! </para>
  subroutine Write2Vectors (lu, x, y, title)
    use Ubik_WriteString_module

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:) :: x, y
    ! </arguments>

    ! Local:
    character(len=132) :: out_string
    integer(int_type) :: i, status

    call Spew (lu, ' ')
    call Spew (lu, title)
    call Spew (lu, ' ')
    do i=1,SIZE(x)
       write(out_string,500) i, x(i), y(i)
       call Spew (lu, out_string)
    end do

500 format (t2, i7, 1p, 2(3x, e16.7))
    return
  end subroutine Write2Vectors

  ! <para>
  ! Write out three vectors.
  ! </para>
  subroutine Write3Vectors (lu, x, y, z, title)
    use Ubik_WriteString_module

    ! <arguments>
    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:) :: x, y, z
    ! </arguments>

    ! Local:
    character(len=132) :: out_string
    integer(int_type) :: i, status

    call Spew (lu, ' ')
    call Spew (lu, title)
    call Spew (lu, ' ')
    do i=1,SIZE(x)
       write(out_string,500) i, x(i), y(i), z(i)
       call Spew (lu, out_string)
    end do

500 format (t2, i7, 1p, 3(3x, e16.7))
    return
  end subroutine Write3Vectors

end module Ubik_Vector_module
