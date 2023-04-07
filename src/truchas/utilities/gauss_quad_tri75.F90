!!
!! GAUSS_QUAD_TRI75
!!
!! Barycentric coordinates and weights for 7-point, 5th order, Gaussian
!! quadrature formula on a triangle.
!!
!! Reference: ???
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Extracted from an old 2D MFE code
!!

module gauss_quad_tri75

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  real(r8), parameter :: a0 = 1.0_r8 / 3.0_r8
  real(r8), parameter :: a1 = 0.1012865073234563388_r8 ! (6 -   Sqrt[15]) / 21
  real(r8), parameter :: b1 = 0.7974269853530873224_r8 ! (9 + 2 Sqrt[15]) / 21
  real(r8), parameter :: a2 = 0.4701420641051150898_r8 ! (6 +   Sqrt[15]) / 21
  real(r8), parameter :: b2 = 0.0597158717897698204_r8 ! (9 - 2 Sqrt[15]) / 21
  
  real(r8), parameter :: w0 = 9.0_r8 / 40.0_r8
  real(r8), parameter :: w1 = 0.1259391805448271526_r8 ! (155 - Sqrt[15]) / 1200
  real(r8), parameter :: w2 = 0.1323941527885061807_r8 ! (155 + Sqrt[15]) / 1200

  ! Barycentric coordinates of the quadrature points  
  real(r8), parameter, public :: GQTRI75_phi(3,7) = reshape(shape=[3,7], &
      source=[a0, a0, a0, &
              b1, a1, a1, &
              a1, b1, a1, &
              a1, a1, b1, &
              b2, a2, a2, &
              a2, b2, a2, &
              a2, a2, b2] )

  ! Corresponding weights
  real(r8), parameter, public :: GQTRI75_wgt(7) = [w0, w1, w1, w1, w2, w2, w2]
  
end module gauss_quad_tri75
