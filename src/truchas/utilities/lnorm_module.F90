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
  use constants_module, only: zero
  use kind_module,      only: real_kind
  implicit none

  ! Private Module
  private

  ! Public Variables

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

    implicit none

    ! Argument List
    real(real_kind), dimension(:), intent(IN)  :: X
    real(real_kind) :: L1NORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    L1NORM = zero

    ! Compute vector norm
    L1NORM = PGSLib_Global_SUM(ABS(X))

    return
  END FUNCTION L1NORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION L2NORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(2) norm of an (nx1) vector; xn=sqrt(sum(X*X)).
    !=======================================================================
    use pgslib_module, only: PGSLib_Global_DOT_PRODUCT

    implicit none

    ! Argument List
    real(real_kind), dimension(:), intent(IN)  :: X
    real(real_kind) :: L2NORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    L2NORM = zero

    ! Compute vector norm
    L2NORM = SQRT(PGSLib_Global_DOT_PRODUCT(X,X))

    return
  END FUNCTION L2NORM

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION LINORM(X)
    !=======================================================================
    ! Purpose(s):
    !   Compute the l(infinity) norm of a (nx1) vector; xn=max(abs(X)).
    !=======================================================================
    use pgslib_module, only: PGSLib_Global_MAXVAL

    implicit none

    ! Argument List
    real(real_kind), dimension(:), intent(IN)  :: X
    real(real_kind) :: LINORM

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize vector norm
    LINORM = zero

    ! Compute vector norm
    LINORM = PGSLib_Global_MAXVAL(ABS(X))

    return
  END FUNCTION LINORM

END MODULE LNORM_MODULE
