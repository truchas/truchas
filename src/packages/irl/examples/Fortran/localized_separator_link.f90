!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! This example covers how to create
! a LocalizedSeparatorLink network to
! cut a polyhedron by a discontinuous
! separator interface (localized inside
! PlanarLocalizers). This is done to directly obtain
! volumetric moments for the polyhedron
! internal/external to the separator for each
! localizer given.

! This will be shown using a Dodecahedron
! that lays across a collection of
! localizers representing RectangularCuboids and
! mimicing a rectilinear mesh.
! These planes will be set directly, and in general
! are capable of representing any mesh made up of
! convex polyhedron.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer :: i,j,k

  real(DP), dimension(1:3,1:8) :: dodecahedron_volume_pts
  type(TagAccVM_SepVM_type) :: phase_moments_LocalizedSeparatorLink

  type(Dod_type) :: dodecahedron_volume
  type(PlanarLoc_type), dimension(-1:1,-1:1,-1:1) :: planar_localizer
  type(PlanarSep_type), dimension(-1:1,-1:1,-1:1) :: planar_separator
  type(LocSepLink_type), dimension(-1:1,-1:1,-1:1) :: localized_separator_link
  integer :: unique_id

  real(DP), dimension(3) :: plane_normal
  real(DP) :: plane_distance

  ! Allocate the actual objects in IRL
  do k = -1,1
    do j = -1,1
      do i = -1,1
        ! These relate to implicit constructors,
        ! allocating the memory for the PlanarLocalizer
        ! and PlanarSeparator
        call new(planar_localizer(i,j,k))
        call new(planar_separator(i,j,k))

        ! During construction of LocalizedSeparatorLink
        !  we need to provide the localizer
        ! and separator it maps to.
        ! This is handled with pointers, so changes to the
        ! objects in the planar_localizer and planar_separator
        ! arrays will be reflected in the localized_separator and
        ! localized_separator_link.
        call new(localized_separator_link(i,j,k),&
                 planar_localizer(i,j,k),&
                 planar_separator(i,j,k))

      end do
    end do
  end do

  ! Let's mimic a [-1.5 x 1.5] uniform cubic grid
  ! with 3 cells in each direction using the
  ! planar_localizers.
  do k = -1,1
    do j = -1,1
      do i = -1,1
        ! This is simply a helper function for this example
        ! It is ``contained`` at the bottom of this program.
        call make_cubic_planar_localizer(i,j,k)
      end do
    end do
  end do

  ! Let's now setup the PlanarSeparators involved, which
  ! will separate the volume of the hexahedron 
  ! after it is localized in the region dictated by each PlanarLocalizer
  ! in the localized_separator and localized_separator_link.
  ! For the sake of simplicity,
  ! lets separate by a flat plane across the middle of
  ! the domain with a normal (0.0,1.0,0.0) and a distance of 0.0.
  ! In general, the PlanarSeparator can have any plane orientation,
  ! and require no continuity between the PlanarSeparators in
  ! the LocalizedSeparatorLink.
  do k = -1,1
    do j = -1,1
      do i = -1,1
        plane_normal = (/0.0_DP,1.0_DP,0.0_DP/)
        plane_distance = 0.0_DP
        call addPlane(planar_separator(i,j,k),plane_normal,plane_distance)
      end do
    end do
  end do


  ! Lastly we need to define the polyhedron to be
  ! separated by the localized_separator_link network.
  ! This will be a trapezoidal prism that intersects all cells
  ! in our mimic'ed mesh.
  dodecahedron_volume_pts(1:3, 1) = (/ 0.25_DP, -1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 2) = (/ 1.00_DP,  1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 3) = (/ 1.00_DP,  1.00_DP,  1.00_DP/)
  dodecahedron_volume_pts(1:3, 4) = (/ 0.25_DP, -1.00_DP,  1.00_DP/)
  dodecahedron_volume_pts(1:3, 5) = (/-0.25_DP, -1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 6) = (/-1.00_DP,  1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 7) = (/-1.00_DP,  1.00_DP,  1.00_DP/)
  dodecahedron_volume_pts(1:3, 8) = (/-0.25_DP, -1.00_DP,  1.00_DP/)

  ! Allocate Dodecahedron storage and construct to use vertices dodecahedron_volume_pts
  call new(dodecahedron_volume)
  call construct(dodecahedron_volume, dodecahedron_volume_pts)

    ! Allocate AccumulatedVolumeMoments_SeparatedVolumeMoments which
  ! will store the moments for each localizer in the LocalizedSeparatorLink
  ! network.
  call new(phase_moments_LocalizedSeparatorLink)

  ! In order to form a network from LocalizedSeparatorLink objects,
  ! they must be informed of their connectivity.
  ! This involves indicating which LocalizedSeparatorLink neighbor
  ! lays on the other side of each plane in the contained PlanarLocalizer, where
  ! boundary planes have neighbors of nullptr.
  ! This will be done using an auxiliary function, ``setupLinking`` at
  ! the end of this program.
  unique_id = 0
  do k = -1,1
    do j = -1,1
      do i = -1,1
        ! The unique id provides a tag number that will identify the
        ! LocalizedSeparatorLink region that each SeparatedMoments<VolumeMoments>
        ! in the Tagged_AccumVM_SepVM_type object belong to.
        unique_id = unique_id + 1
        call setId(localized_separator_link(i,j,k),unique_id)

        ! Setup the linking to neighbors
        call setupLinking(i,j,k)
      end do
    end do
  end do

  ! With linking setup, now perform all the calculation of
  ! localized SeparatedMoments<VolumeMoments> at once. To do this,
  ! only one link in the network needs to be passed. The network
  ! will then be traversed starting from that point.
  ! The moments for a specific PlanarLocalizer region can be obtained
  ! using it's unique Id with the functions getVolumeAtTag and getCentroidAtTag.
  call getNormMoments(dodecahedron_volume, localized_separator_link(-1,-1,-1), &
          phase_moments_LocalizedSeparatorLink)

  write(*,'(A)')
  write(*,'(A)') 'Comparison between Expected and Computed Moments'
  write(*,'(A)') '================================================================'
  write(*,'(A)') ' '
  write(*,'(A)')       'For Localizer(0,-1,-1) '
  write(*,'(A)')       'Internal volume          '
  write(*,'(A,F10.5)') '   Expected:   ', 11.0_DP/64.0_DP
  write(*,'(A,F10.5)') '   Computed:   ', &
    getVolumeAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,-1,-1)),0)
  write(*,'(A)')        'Internal volume centroid'
  write(*,'(A,3F10.5)') '   Expected:  ', (/0.0_DP, -96.0_DP/(11.0_DP*12.0_DP), -0.75_DP/)
  write(*,'(A,3F10.5)') '   Computed:  ', &
    getCentroidAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,-1,-1)),0)
  write(*,'(A)') ' '
  write(*,'(A)')        'For Localizer(0,0,0) '
  write(*,'(A)')        'Internal volume '
  write(*,'(A,F10.5)')  '   Expected:  ', 5.0_DP/32.0_DP + 1.0_DP/3.0_DP
  write(*,'(A,F10.5)')  '   Computed:  ', &
    getVolumeAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,0,0)),0)
  write(*,'(A)')        'Internal volume centroid'
  write(*,'(A,3F10.5)') '   Expected:  ', (/0.0_DP, -104.0_DP/423.0_DP, 0.0_DP/)
  write(*,'(A,3F10.5)') '   Computed:  ', &
    getCentroidAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,0,0)),0)

  write(*,'(A)')        'External volume         '
  write(*,'(A,F10.5)')  '   Expected:  ', 0.5_DP
  write(*,'(A,F10.5)')  '   Computed:  ', &
    getVolumeAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,0,0)),1)
  write(*,'(A)')        'External volume centroid'
  write(*,'(A,3F10.5)') '   Expected:  ', (/0.0_DP, 0.25_DP, 0.0_DP/)
  write(*,'(A,3F10.5)') '   Computed:  ', &
    getCentroidAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,0,0)),1)

contains

  ! This subroutine is used to build the 
  ! planar_localizer objects (pl here)
  ! to represent the hexahedron at the mesh-cell
  ! localizer i,j,k.
  subroutine make_cubic_planar_localizer(i,j,k)
    implicit none
    integer, intent(in) :: i,j,k
    ! use of planar_localizer variables from main program

    real(DP), dimension(1:3) :: lower_pt ! Lower point of hex we're representing
    real(DP), dimension(1:3) :: upper_pt ! Upper point of hex we're representing

    ! Cell centers would be at i,j,k.
    lower_pt = real((/i,j,k/),DP)
    upper_pt = real((/i,j,k/),DP)

    ! Shift by 0.5*dx = 0.5 to get cell vertex positions.
    lower_pt = lower_pt - (/0.5_DP,0.5_DP,0.5_DP/)
    upper_pt = upper_pt + (/0.5_DP,0.5_DP,0.5_DP/)

    ! Now add the planes to the planar localizer (pl)
    call addPlane(planar_localizer(i,j,k),(/-1.0_DP, 0.0_DP, 0.0_DP/),-lower_pt(1)) ! x- face
    call addPlane(planar_localizer(i,j,k),(/ 1.0_DP, 0.0_DP, 0.0_DP/), upper_pt(1)) ! x+ face

    call addPlane(planar_localizer(i,j,k),(/ 0.0_DP,-1.0_DP, 0.0_DP/),-lower_pt(2)) ! y- face
    call addPlane(planar_localizer(i,j,k),(/ 0.0_DP, 1.0_DP, 0.0_DP/), upper_pt(2)) ! y+ face

    call addPlane(planar_localizer(i,j,k),(/ 0.0_DP, 0.0_DP,-1.0_DP/),-lower_pt(3)) ! z- face
    call addPlane(planar_localizer(i,j,k),(/ 0.0_DP, 0.0_DP, 1.0_DP/), upper_pt(3)) ! z+ face

  end subroutine make_cubic_planar_localizer

  subroutine setupLinking(i,j,k)
    implicit none
    integer, intent(in) :: i,j,k
    ! use of localized_separator_link variables from main program

    ! Remember, we set the plane order as x-, x+, y-, y+, z-, z+

    ! First face, x-
    if(i-1 .lt. -1) then
      ! Boundary plane
      call setEdgeConnectivityNull(localized_separator_link(i,j,k),0)
    else
      ! Plane has neighbor
      call setEdgeConnectivity(localized_separator_link(i,j,k),0,localized_separator_link(i-1,j,k))
    end if

    ! Second face, x+
    if(i+1 .gt. 1) then
      ! Boundary plane
      call setEdgeConnectivityNull(localized_separator_link(i,j,k),1)
    else
      ! Plane has neighbor
      call setEdgeConnectivity(localized_separator_link(i,j,k),1,localized_separator_link(i+1,j,k))
    end if

    ! Third face, y-
    if(j-1 .lt. -1) then
      ! Boundary plane
      call setEdgeConnectivityNull(localized_separator_link(i,j,k),2)
    else
      ! Plane has neighbor
      call setEdgeConnectivity(localized_separator_link(i,j,k),2,localized_separator_link(i,j-1,k))
    end if

    ! Fourth face, y+
    if(j+1 .gt. 1) then
      ! Boundary plane
      call setEdgeConnectivityNull(localized_separator_link(i,j,k),3)
    else
      ! Plane has neighbor
      call setEdgeConnectivity(localized_separator_link(i,j,k),3,localized_separator_link(i,j+1,k))
    end if

    ! Fifth face, z-
    if(k-1 .lt. -1) then
      ! Boundary plane
      call setEdgeConnectivityNull(localized_separator_link(i,j,k),4)
    else
      ! Plane has neighbor
      call setEdgeConnectivity(localized_separator_link(i,j,k),4,localized_separator_link(i,j,k-1))
    end if

    ! Sixth face, z+
    if(k+1 .gt. 1) then
      ! Boundary plane
      call setEdgeConnectivityNull(localized_separator_link(i,j,k),5)
    else
      ! Plane has neighbor
      call setEdgeConnectivity(localized_separator_link(i,j,k),5,localized_separator_link(i,j,k+1))
    end if

  end subroutine setupLinking

end program main
