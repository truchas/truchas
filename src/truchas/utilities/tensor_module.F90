!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TENSOR_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define tensor manipulation variables and routines.
  !
  ! Contains: TENSOR_MATRIX
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use parameter_module, only: ndim
  implicit none
  private

  ! Public Variables
  public :: Tensor

  ! Public Subroutines
  public :: TENSOR_MATRIX

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Tensor Matrix
  integer, dimension(ndim,ndim), save :: Tensor

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE TENSOR_MATRIX ()
    !=======================================================================
    ! Purpose(s):
    !   Evaluate the tensor matrix.
    !
    !   ndim = 1: Tensor = [1]
    !
    !                       _    _
    !   ndim = 2: Tensor = | 1, 2 |
    !                      |_2, 1_|
    !
    !                       _       _
    !                      | 1, 2, 3 |
    !   ndim = 3: Tensor = | 3, 1, 2 |
    !                      |_2, 3, 1_|
    !
    !=======================================================================

    ! Argument List

    ! Local Variables
    integer :: i, j, m, n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    Tensor = 0

    ! Tensor Matrix
    do i = 1, ndim
       n = 1
       do j = i, i+ndim-1
          m = MOD(j,ndim)
          if (m > 0) then
             Tensor(i,m) = n
          else
             Tensor(i,j) = n
          end if
          n = n + 1
       end do
    end do

  END SUBROUTINE TENSOR_MATRIX

END MODULE TENSOR_MODULE
