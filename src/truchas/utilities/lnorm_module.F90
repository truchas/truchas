!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE LNORM_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define procedures for computing the l-norm of an vector.
  !
  ! Contains: L1NORM
  !           L2NORM
  !           LINORM
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  ! Public Subroutines
  public :: L1NORM
  public :: L2NORM
  public :: LINORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
CONTAINS

  FUNCTION L1NORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(1) norm of an (nx1) vector; xn=sum(abs(X)).
    !=======================================================================
    use parallel_communication, only: global_sum

    ! Argument List
    real(r8), dimension(:), intent(IN)  :: X
    real(r8) :: L1NORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    L1NORM = 0.0_r8

    ! Compute vector norm
    L1NORM = global_sum(ABS(X))

  END FUNCTION L1NORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION L2NORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(2) norm of an (nx1) vector; xn=sqrt(sum(X*X)).
    !=======================================================================
    use parallel_communication, only: global_dot_product

    ! Argument List
    real(r8), dimension(:), intent(IN)  :: X
    real(r8) :: L2NORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    L2NORM = 0.0_r8

    ! Compute vector norm
    L2NORM = SQRT(global_dot_product(X,X))

  END FUNCTION L2NORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION LINORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(infinity) norm of a (nx1) vector; xn=max(abs(X)).
    !=======================================================================
    use parallel_communication, only: global_maxval

    ! Argument List
    real(r8), dimension(:), intent(IN)  :: X
    real(r8) :: LINORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    LINORM = 0.0_r8

    ! Compute vector norm
    LINORM = global_maxval(ABS(X))

  END FUNCTION LINORM

END MODULE LNORM_MODULE
