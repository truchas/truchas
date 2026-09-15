!!
!! Unit tests for robust two-dimensional truncation-volume calculations.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program test_t2d_truncation_volume

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use t2d_truncation_volume_type
  implicit none

  call test_cartesian_triangle
  call test_quad_corner
  call test_quad_near_corner
  call test_quad_opposite_corner
  call test_nearly_degenerate_triangle
  call test_axisymmetric_degenerate_triangle

contains

  subroutine test_cartesian_triangle

    type(t2d_truncation_volume) :: trunc_vol
    real(r8) :: vol, x(2,3)

    x(:,1) = [0.0_r8, 0.0_r8]
    x(:,2) = [1.0_r8, 0.0_r8]
    x(:,3) = [0.0_r8, 1.0_r8]
    call trunc_vol%init(x, [1.0_r8, 0.0_r8], .false.)
    vol = trunc_vol%volume(0.5_r8)
    call require_volume(vol, 0.375_r8, 0.5_r8, 'Cartesian triangle truncation')

  end subroutine test_cartesian_triangle


  !! A line through the first corner cuts the unit-square swept polygon into
  !! two nondegenerate pieces. Two cell edges report the same corner.
  subroutine test_quad_corner

    type(t2d_truncation_volume) :: trunc_vol
    real(r8) :: normal(2), rho, vol, x(2,4)

    call unit_square(x)
    normal = [1.0_r8, -0.5_r8]
    normal = normal/norm2(normal)
    rho = dot_product(normal, x(:,1))
    call trunc_vol%init(x, normal, .false.)
    vol = trunc_vol%volume(rho)
    call require_volume(vol, 0.25_r8, 1.0_r8, 'quad corner intersection')

  end subroutine test_quad_corner


  !! Move the preceding line only a few floating-point spacings from the
  !! corner. It must remain a finite, bounded perturbation of the same cut.
  subroutine test_quad_near_corner

    type(t2d_truncation_volume) :: trunc_vol
    real(r8) :: normal(2), rho, vol, x(2,4)

    call unit_square(x)
    normal = [1.0_r8, -0.5_r8]
    normal = normal/norm2(normal)
    rho = dot_product(normal, x(:,1)) + 8.0_r8*epsilon(1.0_r8)
    call trunc_vol%init(x, normal, .false.)
    vol = trunc_vol%volume(rho)
    call require_volume(vol, 0.25_r8, 1.0_r8, 'near-corner quad intersection')

  end subroutine test_quad_near_corner


  !! The cutting line passes through the corner opposite the first vertex and
  !! leaves the complete swept polygon behind it.
  subroutine test_quad_opposite_corner

    type(t2d_truncation_volume) :: trunc_vol
    real(r8) :: normal(2), rho, vol, x(2,4)

    call unit_square(x)
    normal = [1.0_r8, 1.0_r8]
    normal = normal/norm2(normal)
    rho = dot_product(normal, x(:,3))
    call trunc_vol%init(x, normal, .false.)
    vol = trunc_vol%volume(rho)
    call require_volume(vol, 1.0_r8, 1.0_r8, 'opposite-corner quad intersection')

  end subroutine test_quad_opposite_corner


  !! This sliver reproduces the geometry that caused Kahan's side-length
  !! formula to signal an invalid operation in the Campbell flow example.
  subroutine test_nearly_degenerate_triangle

    type(t2d_truncation_volume) :: trunc_vol
    real(r8) :: vol, x(2,3)

    x(:,1) = [3.4752654155967111e-3_r8, 3.611974851854753e-1_r8]
    x(:,2) = [3.4752654155967115e-3_r8, 3.611974851854754e-1_r8]
    x(:,3) = [3.4756578947368428e-3_r8, 3.612500000000000e-1_r8]
    call trunc_vol%init(x, [1.0_r8, 0.0_r8], .false.)
    vol = trunc_vol%volume(1.0_r8)
    call require_volume(vol, 0.0_r8, 0.0_r8, 'nearly degenerate Cartesian triangle')

  end subroutine test_nearly_degenerate_triangle


  subroutine test_axisymmetric_degenerate_triangle

    type(t2d_truncation_volume) :: trunc_vol
    real(r8) :: vol, x(2,3)

    x(:,1) = [1.0_r8, 0.0_r8]
    x(:,2) = [1.0_r8, 0.0_r8]
    x(:,3) = [2.0_r8, 1.0_r8]
    call trunc_vol%init(x, [1.0_r8, 0.0_r8], .true.)
    vol = trunc_vol%volume(3.0_r8)
    call require_volume(vol, 0.0_r8, 0.0_r8, 'degenerate axisymmetric triangle')

  end subroutine test_axisymmetric_degenerate_triangle


  subroutine unit_square(x)

    real(r8), intent(out) :: x(2,4)

    x(:,1) = [0.0_r8, 0.0_r8]
    x(:,2) = [1.0_r8, 0.0_r8]
    x(:,3) = [1.0_r8, 1.0_r8]
    x(:,4) = [0.0_r8, 1.0_r8]

  end subroutine unit_square


  subroutine require_volume(actual, expected, total, message)

    real(r8), intent(in) :: actual, expected, total
    character(*), intent(in) :: message

    if (.not.ieee_is_finite(actual)) error stop 'FAIL: non-finite ' // message
    if (actual < -1.0e-14_r8 .or. actual > total+1.0e-14_r8) &
      error stop 'FAIL: unbounded ' // message
    if (abs(actual-expected) > 1.0e-13_r8) error stop 'FAIL: inaccurate ' // message

  end subroutine require_volume

end program test_t2d_truncation_volume
