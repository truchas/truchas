!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
    use pgslib_module, only: PGSLib_Global_SUM

    ! Argument List
    real(r8), dimension(:), intent(IN)  :: X
    real(r8) :: L1NORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    L1NORM = 0.0_r8

    ! Compute vector norm
    L1NORM = PGSLib_Global_SUM(ABS(X))

  END FUNCTION L1NORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION L2NORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(2) norm of an (nx1) vector; xn=sqrt(sum(X*X)).
    !=======================================================================
    use pgslib_module, only: PGSLib_Global_DOT_PRODUCT

    ! Argument List
    real(r8), dimension(:), intent(IN)  :: X
    real(r8) :: L2NORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    L2NORM = 0.0_r8

    ! Compute vector norm
    L2NORM = SQRT(PGSLib_Global_DOT_PRODUCT(X,X))

  END FUNCTION L2NORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION LINORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(infinity) norm of a (nx1) vector; xn=max(abs(X)).
    !=======================================================================
    use pgslib_module, only: PGSLib_Global_MAXVAL

    ! Argument List
    real(r8), dimension(:), intent(IN)  :: X
    real(r8) :: LINORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    LINORM = 0.0_r8

    ! Compute vector norm
    LINORM = PGSLib_Global_MAXVAL(ABS(X))

  END FUNCTION LINORM

END MODULE LNORM_MODULE
