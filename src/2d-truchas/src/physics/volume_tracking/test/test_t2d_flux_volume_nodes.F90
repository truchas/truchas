!!
!! Unit tests for swept-volume polygons that intersect cell corners.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program test_t2d_flux_volume_nodes

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use t2d_cell_geom_vof_type, only: t2d_cell_geom
  use t2d_flux_volume_nodes_function, only: flux_volume_nodes
  use t2d_truncation_volume_type, only: t2d_truncation_volume
  implicit none

  call test_first_opposite_corner
  call test_near_corner
  call test_second_opposite_corner

contains

  subroutine test_first_opposite_corner

    type(t2d_cell_geom) :: cell
    real(r8) :: x(2,4)

    x(:,1) = [0.0_r8, 0.0_r8]
    x(:,2) = [1.0_r8, 0.0_r8]
    x(:,3) = [1.0_r8, 0.25_r8]
    x(:,4) = [0.0_r8, 1.0_r8]
    call init_cell(cell, x)
    call check_corner_flux(cell, 0.25_r8, 'first opposite corner')

  end subroutine test_first_opposite_corner


  subroutine test_near_corner

    type(t2d_cell_geom) :: cell
    real(r8) :: delta, x(2,4)

    x(:,1) = [0.0_r8, 0.0_r8]
    x(:,2) = [1.0_r8, 0.0_r8]
    x(:,3) = [1.0_r8, 0.25_r8]
    x(:,4) = [0.0_r8, 1.0_r8]
    call init_cell(cell, x)
    delta = 8.0_r8*epsilon(1.0_r8)
    call check_corner_flux(cell, 0.25_r8+delta, 'near first opposite corner')

  end subroutine test_near_corner


  subroutine test_second_opposite_corner

    type(t2d_cell_geom) :: cell
    real(r8) :: x(2,4)

    x(:,1) = [0.0_r8, 0.0_r8]
    x(:,2) = [1.0_r8, 0.0_r8]
    x(:,3) = [1.0_r8, 1.0_r8]
    x(:,4) = [0.0_r8, 0.25_r8]
    call init_cell(cell, x)
    call check_corner_flux(cell, 0.25_r8, 'second opposite corner')

  end subroutine test_second_opposite_corner


  subroutine init_cell(cell, x)

    type(t2d_cell_geom), intent(out) :: cell
    real(r8), intent(in) :: x(2,4)

    real(r8) :: area, face_area(4), normal(2,4)

    area = 0.5_r8*abs(sum(x(1,:)*cshift(x(2,:),1)-cshift(x(1,:),1)*x(2,:)))
    face_area = 1.0_r8
    normal = 0.0_r8
    normal(:,1) = [0.0_r8, -1.0_r8]
    call cell%init(x, area, face_area, normal)

  end subroutine init_cell


  subroutine check_corner_flux(cell, target, message)

    type(t2d_cell_geom), intent(in) :: cell
    real(r8), intent(in) :: target
    character(*), intent(in) :: message

    type(t2d_truncation_volume) :: swept
    real(r8) :: nodes(2,4), vol
    integer :: nnode

    call flux_volume_nodes(1, cell, target, target, 1.0e-12_r8, nodes, nnode, .false.)
    if (nnode /= 4) error stop 'FAIL: duplicate node in ' // message
    call swept%init(nodes(:,:nnode), [1.0_r8, 0.0_r8], .false.)
    vol = swept%volume(2.0_r8)
    if (.not.ieee_is_finite(vol)) error stop 'FAIL: non-finite ' // message
    if (vol < 0.0_r8 .or. vol > cell%volume) error stop 'FAIL: unbounded ' // message
    if (abs(vol-target) > 1.0e-10_r8) error stop 'FAIL: inaccurate ' // message

  end subroutine check_corner_flux

end program test_t2d_flux_volume_nodes
