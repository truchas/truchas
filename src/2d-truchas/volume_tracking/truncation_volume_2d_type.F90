!!
!! This module defines the truncation volume type required for Volume-Of-Fluid
!! calculations in 2D Cartesian and axisymmetric co-ordinates.
!!
!! Aditya K. Pandare <apandare@lanl.gov>
!! Jan 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truncation_volume_2d_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use geom_axisymmetric
  implicit none

  type, public :: truncation_volume
    private
    integer :: nfc
    real(r8), allocatable :: nodex(:,:)
    ! set of nodes for triangular elements resulting from element splitting:
    ! first index corresponds to how many triangles resulted from the splitting;
    ! second index is the (x,y) coordinate index;
    ! third index is the local node number in the component triangle corresponding to the
    ! first index
    real(r8), allocatable :: node_set(:,:,:)
    real(r8) :: plane_normal(2)
    logical :: is_axisymmetric
  contains
    procedure :: init => init_truncation_volume
    procedure :: volume
    procedure, private :: split_quad4
    procedure, private :: trunc_tri_volume
  end type truncation_volume

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_truncation_volume(this, nodex, plane_normal, axisym)

    use cell_topology_2d

    class(truncation_volume), intent(out) :: this
    real(r8), intent(in) :: nodex(:,:), plane_normal(:)
    logical, intent(in) :: axisym

    allocate(this%nodex(2, size(nodex, dim=2)))
    this%nodex = nodex
    this%plane_normal = plane_normal

    this%is_axisymmetric = axisym

    ! get the number of faces for this cell type and split quadrilaterals into triangles
    select case (size(nodex, dim=2))
    case (3) ! triangle
      this%nfc = 3
      allocate(this%node_set(1,2,3))
      this%node_set(1,:,:) = this%nodex
    case (4) ! quadrilateral
      this%nfc = 4
      allocate(this%node_set(2,2,3))
      call this%split_quad4()
    case default
      call TLS_panic('unaccounted topology in truncation_volume_2d_type')
    end select

  end subroutine init_truncation_volume

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! calculates the truncation volume
  real(r8) function volume(this, plane_rho)

    class(truncation_volume), intent(in) :: this
    real(r8), intent(in) :: plane_rho

    integer :: i
    real(r8), allocatable :: vol_sub(:)

    volume = 0.0_r8

    select case (this%nfc)
    case (3) ! triangle
      allocate(vol_sub(1))
    case (4) ! quadrilateral
      allocate(vol_sub(2))
    case default
      call TLS_panic('unaccounted topology in truncation_volume_2d_type')
    end select

    ! determine intersection of plane with component triangles and get truncated volume
    do i = 1, size(this%node_set, dim=1)
      vol_sub(i) = this%trunc_tri_volume(this%node_set(i,:,:), plane_rho)
    end do !i

    ! total volume behind interface plane (truncation volume for this element) is
    ! the sum of truncated volumes of component triangles
    volume = sum(vol_sub)

  end function volume

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! splits a quadratic element (quad4) into two triangular elements (tri3)
  subroutine split_quad4(this)
    class(truncation_volume), intent(inout) :: this
    this%node_set(1,:,1:3) = this%nodex(:,1:3)
    this%node_set(2,:,1) = this%nodex(:,1)
    this%node_set(2,:,2) = this%nodex(:,3)
    this%node_set(2,:,3) = this%nodex(:,4)
  end subroutine split_quad4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! checks intersection of plane with component triangle, and then calculates the
  ! truncated volume of the component triangle
  real(r8) function trunc_tri_volume(this, node_set, plane_rho)

    use cell_geometry
    use plane_2d_type

    class(truncation_volume), intent(in) :: this
    real(r8), intent(in) :: node_set(:,:), plane_rho

    integer :: f, icount, icut, on_point(3), fid(3)
    real(r8) :: xf(2,2), xint(2,3), xt(2,3), vol_full_tri, vol_sub_tri
    type(plane) :: int_plane
    logical :: entire_element, cut_plane(3), f_intersection(3)

    int_plane%rho = plane_rho
    int_plane%normal = this%plane_normal

    if (this%is_axisymmetric) then
      vol_full_tri = polygon_volume_axisym(node_set)
    else
      vol_full_tri = tri_area(node_set)
    end if

    ! 1. check intersection of plane with each face of triangle.
    ! this loop also detects if the face was intersected at its end-points (nodes) or
    ! if it was cut into two (i.e. intersected between its end-points). The cut_plane
    ! and icut data keeps track of this.
    xt = 0.0_r8
    on_point = 0
    icut = 0
    cut_plane(:) = .false.
    f_intersection(1:3) = .false.
    do f = 1, 3
      xf(:,1) = node_set(:,f)
      xf(:,2) = node_set(:,mod(f,3)+1)
      f_intersection(f) = int_plane%intersects(xf)
      if (f_intersection(f)) then
        call int_plane%intersection_point(xint(:,f), on_point(f), xf)
        ! if this face has been cut, mark it as a cut_plane
        if (on_point(f) == 0) then
          cut_plane(f) = .true.
          icut = icut + 1
        end if
      end if
    end do

    INSIST(.not.all(cut_plane))
    INSIST(icut < 3)

    ! 2. if no cut_planes have been found, then, the entire cell is behind the plane
    if (any(cut_plane)) then
      entire_element = .false.
    else
      INSIST(icut == 0)
      entire_element = .true.
    end if

    ! 3. determine the vertices of the truncated sub-triangle, and its area
    if (entire_element) then
    ! if the entire element is behind/in-front-of the plane, use the entire triangle
      xt = node_set
      vol_sub_tri = vol_full_tri

    else
    ! if cut_planes have been found, then, determine the truncated triangle
      fid(:) = 0
      do f = 1, 3
        if (cut_plane(f)) fid(f) = f
      end do !f

      if (icut == 2) then
      ! two faces have been cut
        icount = 0
        do f = 1, 3
          if (fid(f) /= 0) then
            icount = icount + 1
            xt(:,icount) = xint(:,fid(f))

            ! if two intersections have already been accounted for, then assign the
            ! third vertex of the truncated triangle
            if (icount == 2 .and. f == 3 .and. fid(2) == 0) then
              xt(:,3) = node_set(:,1)
            else if (icount == 2 .and. f == 3 .and. fid(2) == 2) then
              xt(:,3) = node_set(:,3)
            else if (icount == 2) then
              xt(:,3) = node_set(:,fid(f))
            end if
          end if
        end do !f

      else
      ! one face has been cut
        INSIST(icut == 1)
        do f = 1, 3
          if (fid(f) /= 0) then
            xt(:,1) = xint(:,fid(f))
            xt(:,2) = node_set(:,mod(fid(f),3)+1)
            xt(:,3) = node_set(:,mod(fid(f)+1,3)+1)
            exit
          end if
        end do !f

      end if

      ! truncated sub-triangle area based on above vertices
      if (this%is_axisymmetric) then
        vol_sub_tri = polygon_volume_axisym(xt)
      else
        vol_sub_tri = tri_area(xt)
      end if

    end if

    ! 4. determine volume 'behind' the plane based on the signed_distance function
    trunc_tri_volume = vol_sub_tri
    if (int_plane%signed_distance(sum(xt(:,:), dim=2)/3.0_r8) > 0.0_r8) then
      if (entire_element) then
        trunc_tri_volume = 0.0_r8
      else
        trunc_tri_volume = vol_full_tri-vol_sub_tri
      end if
    end if

  end function trunc_tri_volume

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module truncation_volume_2d_type
