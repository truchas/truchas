! -*- Mode: f90 -*-

!**********************************************************************
!    $Id: SetUpSaad.F,v 1.9 2003/10/19 22:28:18 turner Exp $
!
!    Subroutines to set up Saad test problem; Ax = b, where:
!
!                   B  -I                     4  a
!                  -I   B  -I                 c  4  a
!             A =      -I   B  -I         B =    c  4  a
!                          -I   B  -I               c  4  a
!                              -I   B                  c  4
!
!    where a = -1 + d and c = -1 - d.  The right-hand side is set to
!    b = Af, where f = [1,...,1]^T, so the solution is x = f.
!
!**********************************************************************
module SetUpSaad_module
  use UbikMatrix_module
  implicit none

  interface SetUpSaad
     module procedure SetUpSaad_Full
     module procedure SetUpSaad_ELL
  end interface

  private
  public :: SetUpSaad

contains
  
  !*******************************************************************
  !
  !    Subroutine to set up Saad test problem in full format.
  !
  !    <PARAMETER LIST>
  !
  !     Input:
  !      nsizeb - size of block matrix B above
  !      dvalue - parameter d, described above
  !
  !     Output:
  !      a - coeff. matrix
  !      b - source vector
  !      ref - reference (true) solution
  !
  !     Used: 
  !      x
  !
  !**********************************************************************
  subroutine SetUpSaad_Full (nsizeb, dvalue, a, b, ref, x)

    ! Args:
    integer(Ubik_int_type), intent(in) :: nsizeb
    real(Ubik_real_type), intent(in) :: dvalue
    real(Ubik_real_type), intent(inout), dimension(:,:) :: a
    real(Ubik_real_type), intent(out), dimension(:) :: x, b, ref

    ! Local:
    integer(Ubik_int_type) :: i, nsizea
    real(Ubik_real_type) :: aa, bb

    nsizea = SIZE(a, dim=1)

    aa = -1.0 + dvalue
    bb = -1.0 - dvalue
    a = 0.0d0
    do i=1,nsizea

       ! Main diagonal.
       a(i,i) = 4.0

       ! Inner diagonals.
       if (MOD(i,nsizeb) /= 0) then
          if (i < nsizea) then
             a(i,i+1) = aa
             a(i+1,i) = bb
          end if
       end if

       ! Outer upper diagonal.
       if (i <= nsizea-nsizeb) then
          a(i,i+nsizeb) = -1.0
       end if

       ! Outer lower diagonal.
       if (i > nsizeb) then
          a(i,i-nsizeb) = -1.0
       end if
    end do

!    call Ubik_Spew (6, a, 'blorp')

    ! Set up source vector.
    x = 1.0
    b = Ubik_MatVec(a, x)

    ! Set reference solution.
    ref = 1.0

    return
  end subroutine SetUpSaad_Full

  !**********************************************************************
  !
  !    Subroutine to set up Saad test problem in ELL format.
  !
  !    <PARAMETER LIST>
  !
  !     Input:
  !      nsizeb - size of block matrix B above
  !      dvalue - parameter d, described above
  !
  !     Output:
  !      a - coeff. matrix in ELL format
  !      b - source vector
  !      ref - reference (true) solution
  !
  !     Used: 
  !      x
  !
  !    <SUBROUTINES REQUIRED>
  !
  !     Ubik_y_eq_Ax_ELL
  !
  !**********************************************************************
  subroutine SetUpSaad_ELL (nsizeb, dvalue, a, b, ref, x)

    ! Args:
    integer(Ubik_int_type), intent(in) :: nsizeb
    real(Ubik_real_type), intent(in) :: dvalue
    type(Ubik_ELL_matrix), intent(inout) :: a
    real(Ubik_real_type), intent(out), dimension(:) :: x, b, ref

    ! Local:
    integer(Ubik_int_type) :: i, maxnz, nsizea
    real(Ubik_real_type) :: aa, bb

    nsizea = SIZE(a%values, dim=1)
    maxnz = SIZE(a%values, dim=2)

    aa = -1.0 + dvalue
    bb = -1.0 - dvalue
    a = 0.0d0
    a%values(1:nsizea,1) = 4.0
    do i=1,nsizea

       ! Main diagonal.
       a%map(i,1) = i

       ! Inner diagonals.
       if (MOD(i,nsizeb) /= 0) then
          if (i < nsizea) then
             a%values(i,2) = aa
             a%map(i,2) = i+1
             a%values(i+1,3) = bb
             a%map(i+1,3) = i
          end if
       end if

       ! Outer upper diagonal.
       if (maxnz >= 4) then
          if (i <= nsizea-nsizeb) then
             a%values(i,4) = -1.0
             a%map(i,4) = i + nsizeb
          end if

          ! Outer lower diagonal.
          if (maxnz >= 5) then
             if (i > nsizeb) then
                a%values(i,5) = -1.0
                a%map(i,5) = i - nsizeb
             end if
          end if
       end if
    end do

    ! Set up source vector.
    x = 1.0
    b = Ubik_MatVec(a, x)

    ! Set reference solution.
    ref = 1.0

    return
  end subroutine SetUpSaad_ELL

end module SetUpSaad_module
