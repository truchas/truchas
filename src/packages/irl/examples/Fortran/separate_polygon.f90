!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! This example demonstrates the ability
! to create an arbitrary polygon and then
! return individual triangles from its surface.
! The integration of surface moments for localized
! regions is also demonstrated.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)

  type(Poly_type) :: polygon
  type(Tri_type) :: triangle
  type(PlanarLoc_type) :: planar_localizer
  real(DP), dimension(3,5) :: polygon_pts

  integer :: n
  real(DP) :: surface_area, triangle_area
  real(DP), dimension(4) :: plane

  ! Initialize the different polygon objects
  call new(polygon)
  call new(triangle)


  ! Set the polygon to be a pentagon in the XZ Plane
  polygon_pts(:,1) = (/0.0_DP, 0.0_DP, 0.0_DP/)
  polygon_pts(:,2) = (/0.0_DP, 0.0_DP, 1.0_DP/)
  polygon_pts(:,3) = (/0.5_DP, 0.0_DP, 1.5_DP/)
  polygon_pts(:,4) = (/1.0_DP, 0.0_DP, 1.0_DP/)
  polygon_pts(:,5) = (/1.0_DP, 0.0_DP, 0.0_DP/)
  call construct(polygon, 5, polygon_pts)

  ! Calculate the plane the polygon exists on
  ! Following the right-hand-rule, this should be
  ! n = (0.0, 1.0, 0.0) and d = 0.0
  call calculateAndSetPlaneOfExistence(polygon)
  plane =  getPlaneOfExistence(polygon)
  write(*,'(A)')
  write(*,'(A)') '          Polygon usage demonstration           '
  write(*,'(A)') '================================================'
  write(*,'(A, 4F10.5)') ' Expected existence plane   : ', (/0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP/)
  write(*,'(A, 4F10.5)') ' Calculated existence plane : ', plane(1:4)

  ! Surface moments for this polygon can be computed
  ! It can also be used along with getVolumeMoments
  ! to localize the polygon first as well.
  surface_area = calculateVolume(polygon)
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Correct surface area    :  ', 1.0_DP + 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated surface area :  ', surface_area

  ! Localize the top triangle of the pentagon
  call new(planar_localizer)
  call setFromRectangularCuboid(planar_localizer, (/0.0_DP, -0.1_DP, 1.0_DP/), (/1.0_DP, 0.1_DP, 1.5_DP/))
  call getNormMoments(polygon, planar_localizer, surface_area)
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Correct localized surface area    :  ', 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated localized surface area :  ', surface_area


  ! Individual triangles can also be obtained from
  ! polygon as well.
  surface_area = 0.0_DP
  do n = 1, getNumberOfSimplicesInDecomposition(polygon)
    call getSimplexFromDecomposition(polygon, n-1, triangle)
    triangle_area = calculateVolume(triangle)
    surface_area = surface_area + triangle_area
  end do
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Correct surface area                   :  ', 1.0_DP + 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated surface area from triangles :  ', surface_area

  surface_area = 0.0_DP
  do n = 1, getNumberOfSimplicesInDecomposition(polygon)
    call getSimplexFromDecomposition(polygon, n-1, triangle)
    call getNormMoments(triangle, planar_localizer, triangle_area)
    surface_area = surface_area + triangle_area
  end do
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Correct surface area                             :  ', 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated localized surface area from triangles :  ', surface_area

end program main

