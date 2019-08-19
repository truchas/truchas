!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! In this example, the moment-of-fluid
! reconstruction (MoF) method implemented inside
! IRL will be demonstrated. First, a PlanarSeparator
! will be directly specified and used to obtain
! SeperatedVolumeMoments for a cell. This will
! then be used to perform a MoF reconstruction.
! The moments from the obtained reconstruction
! are compared to those given.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  real(DP), dimension(1:3,1:8) :: rectangular_cuboid_pts
  type(SepVM_type) :: phase_moments
  type(SepVM_type) :: found_phase_moments_2D, found_phase_moments_3D
  type(SepVM_type) :: found_phase_moments_2D_SetWeight, found_phase_moments_3D_SetWeight
  real(DP), dimension(1:3) :: normal
  type(RectCub_type) :: rectangular_cuboid
  type(PlanarSep_type) :: correct_planar_separator
  type(PlanarSep_type) :: found_planar_separator_2D, found_planar_separator_3D

  ! Define unit-cubic cell
  rectangular_cuboid_pts(1:3, 1) = (/ 0.5_DP, -0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 2) = (/ 0.5_DP,  0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 3) = (/ 0.5_DP,  0.5_DP,  0.5_DP/)
  rectangular_cuboid_pts(1:3, 4) = (/ 0.5_DP, -0.5_DP,  0.5_DP/)
  rectangular_cuboid_pts(1:3, 5) = (/-0.5_DP, -0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 6) = (/-0.5_DP,  0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 7) = (/-0.5_DP,  0.5_DP,  0.5_DP/)
  rectangular_cuboid_pts(1:3, 8) = (/-0.5_DP, -0.5_DP,  0.5_DP/)

  ! Allocate RectangularCuboid storage and construct to be unit cell
  call new(rectangular_cuboid)
  call construct(rectangular_cuboid, rectangular_cuboid_pts)


  ! Allocate the three PlanarSeparators we'll use
  call new(correct_planar_separator)
  call new(found_planar_separator_2D)
  call new(found_planar_separator_3D)

  ! Allocate SeparatedMoments<VolumeMoments> objects to store cutting results
  call new(phase_moments)
  call new(found_phase_moments_2D)
  call new(found_phase_moments_3D)
  call new(found_phase_moments_2D_SetWeight)
  call new(found_phase_moments_3D_SetWeight)


  ! Define interface reconstruction representing a slanted
  ! line across the cell
  normal = (/0.5_DP*sqrt(3.0_DP), 0.5_DP,0.0_DP/)
  call addPlane(correct_planar_separator,normal,-0.15_DP)

  ! Perform cutting through IRL library to obtain resulting SeparatedMoments<VolumeMoments>
  call getNormMoments(rectangular_cuboid, correct_planar_separator, phase_moments)

  ! Now lets use MoF to obtain an interface reconstruction given the cell
  ! (rectangular_cuboid) and the SeparatedPhaseMoments (phase_moments)

  ! First, treat as a 2D case
  call reconstructMOF2D(rectangular_cuboid, phase_moments, found_planar_separator_2D)

  ! Let's now perform 3D MoF
  call reconstructMOF3D(rectangular_cuboid, phase_moments, found_planar_separator_3D)

  ! Now obtain the SeparatedMoments<VolumeMoments> for each
  ! reconstruction in order to compare to the given SeparatedMoments<VolumeMoments>
  call getNormMoments(rectangular_cuboid, found_planar_separator_2D, found_phase_moments_2D)
  call getNormMoments(rectangular_cuboid, found_planar_separator_3D, found_phase_moments_3D)

  ! We can also supply our own weights for MoF to use, allowing us
  ! to bias matching one centroid over the other. Here, let's just 
  ! try to match the internal centroid.
  ! First, treat as a 2D case
  call reconstructMOF2D(rectangular_cuboid, &
      phase_moments, 1.0_DP, 0.0_DP,found_planar_separator_2D)

  ! Let's now perform 3D MoF
  call reconstructMOF3D(rectangular_cuboid, &
      phase_moments, 1.0_DP, 0.0_DP,found_planar_separator_3D)

  ! Once again obtain the SeparatedMoments<VolumeMoments> to compare against
  call getNormMoments(rectangular_cuboid, found_planar_separator_2D, found_phase_moments_2D_SetWeight)
  call getNormMoments(rectangular_cuboid, found_planar_separator_3D, found_phase_moments_3D_SetWeight)

  write(*,'(A)')
  write(*,'(A)') 'Comparison between given and computed results'
  write(*,'(A)') '================================================'
  write(*,'(A)') 'Internal volume centroids '
  write(*,'(A,3F10.5)') '   Given                         : ', getCentroid(phase_moments,0)
  write(*,'(A,3F10.5)') '   Computed (2D, Default Weight) : ',getCentroid(found_phase_moments_2D,0)
  write(*,'(A,3F10.5)') '   Computed (3D, Default Weight) : ',getCentroid(found_phase_moments_3D,0)
  write(*,'(A,3F10.5)') '   Computed (2D, Set Weight)     : ',getCentroid(found_phase_moments_2D_SetWeight,0)
  write(*,'(A,3F10.5)') '   Computed (3D, Set Weight)     : ',getCentroid(found_phase_moments_3D_SetWeight,0)
  write(*,'(A)') 'External volume centroids '
  write(*,'(A,3F10.5)') '   Given                         : ', getCentroid(phase_moments,1)
  write(*,'(A,3F10.5)') '   Computed (2D, Default Weight) : ',getCentroid(found_phase_moments_2D,1)
  write(*,'(A,3F10.5)') '   Computed (3D, Default Weight) : ',getCentroid(found_phase_moments_3D,1)
  write(*,'(A,3F10.5)') '   Computed (2D, Set Weight)     : ',getCentroid(found_phase_moments_2D_SetWeight,1)
  write(*,'(A,3F10.5)') '   Computed (3D, Set Weight)     : ',getCentroid(found_phase_moments_3D_SetWeight,1)
  write(*,'(A)')

end program main
