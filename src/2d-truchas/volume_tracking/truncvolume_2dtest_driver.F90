#include "f90_assert.fpp"

program truncvolume_2dtest_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use truncation_volume_2d_type
  use plane_2d_type
  implicit none

  real*8 :: xn(2,4), vp
  real*8 :: t_start, t_end
  type(plane) :: intplane
  type(truncation_volume) :: trunc_vol

  t_start = 0.0_r8
  t_end = 0.0_r8

  call cpu_time(t_start)

  !! cell vertices
  xn(:,1) = [0.0_r8, 0.0_r8]
  xn(:,2) = [1.0_r8, 0.0_r8]
  xn(:,3) = [1.0_r8, 1.0_r8]
  xn(:,4) = [0.0_r8, 1.0_r8]

  !! interface plane-1 definition
  intplane%normal = [1.0_r8, -0.5_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,1))

  write(*,*) "Plane-1 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-1 constant: ", intplane%rho

  !! truncation volume-1
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-1: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.25_r8) < 1e-13_r8)

  !! interface plane-2 definition
  intplane%normal = [0.5_r8, -1.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,1))

  write(*,*) "Plane-2 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-2 constant: ", intplane%rho

  !! truncation volume-2
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-2: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.75_r8) < 1e-13_r8)

  !! interface plane-3 definition
  intplane%normal = [1.0_r8, 0.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,1))

  write(*,*) "Plane-3 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-3 constant: ", intplane%rho

  !! truncation volume-3
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-3: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.0_r8) < 1e-13_r8)

  !! interface plane-4 definition
  intplane%normal = [-1.0_r8, 0.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,1))

  write(*,*) "Plane-4 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-4 constant: ", intplane%rho

  !! truncation volume-4
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-4: ", vp
  write(*,*) " "
  INSIST(abs(vp-1.0_r8) < 1e-13_r8)

  !! interface plane-5 definition
  intplane%normal = [1.0_r8, -1.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,1))

  write(*,*) "Plane-5 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-5 constant: ", intplane%rho

  !! truncation volume-5
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-5: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.5_r8) < 1e-13_r8)

  !! interface plane-6 definition
  intplane%normal = [1.0_r8, 0.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, [1.0_r8/3.0_r8, 0.0_r8])

  write(*,*) "Plane-6 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-6 constant: ", intplane%rho

  !! truncation volume-6
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-6: ", vp
  write(*,*) " "
  INSIST(abs(vp-1.0_r8/3.0_r8) < 1e-13_r8)

  !! interface plane-7 definition
  intplane%normal = [-1.0_r8, 0.5_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,1))

  write(*,*) "Plane-7 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-7 constant: ", intplane%rho

  !! truncation volume-7
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-7: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.75_r8) < 1e-13_r8)

  !! interface plane-8 definition
  intplane%normal = [1.0_r8, 1.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,2))

  write(*,*) "Plane-8 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-8 constant: ", intplane%rho

  !! truncation volume-8
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-8: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.5_r8) < 1e-13_r8)

  !! interface plane-9 definition
  intplane%normal = [1.0_r8, 1.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,3))

  write(*,*) "Plane-9 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-9 constant: ", intplane%rho

  !! truncation volume-9
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-9: ", vp
  write(*,*) " "
  INSIST(abs(vp-1.0_r8) < 1e-13_r8)

  !! interface plane-10 definition
  intplane%normal = [-1.0_r8, -1.0_r8]
  intplane%normal = intplane%normal / norm2(intplane%normal)
  intplane%rho = dot_product(intplane%normal, xn(:,3))

  write(*,*) "Plane-10 normal:   ", intplane%normal(1), intplane%normal(2)
  write(*,*) "Plane-10 constant: ", intplane%rho

  !! truncation volume-10
  call trunc_vol%init(xn, intplane%normal)
  vp = trunc_vol%volume(intplane%rho)

  write(*,*) "Truncation volume behind plane-10: ", vp
  write(*,*) " "
  INSIST(abs(vp-0.0_r8) < 1e-13_r8)

  call cpu_time(t_end)

  write(*,*) "Runtime: ", t_end-t_start


end program truncvolume_2dtest_driver
