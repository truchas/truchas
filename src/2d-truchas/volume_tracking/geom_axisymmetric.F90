!!
!! This module provides various functions necessary to transform 2D cartesian geometric
!! entites (volumes or face-areas) to 2D axisymmetric ones. All axisymmetric entities are
!! assumed to be calculated for a sweep of 2*pi rad.
!!
!! Aditya K. Pandare <apandare@lanl.gov>
!! Jan 2020
!!

#include "f90_assert.fpp"

module geom_axisymmetric

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_2d_mesh_type
  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mesh_axisymmetry_mod(mesh)

    type(unstr_2d_mesh), target :: mesh

    integer :: i

    ! transform cell-volumes to axisymmetrically swept ones
    do i = 1, mesh%ncell

      associate (cn => mesh%cnode(mesh%cstart(i):mesh%cstart(i+1)-1))
        mesh%volume(i) = polygon_volume_axisym(mesh%x(:,cn))
      end associate

    end do !i

    ! transform face-areas and area-weighted normals to axisymmetrically swept ones
    do i = 1, mesh%nface
      mesh%normal(:,i) = mesh%normal(:,i)/mesh%area(i)
      mesh%area(i) = face_area_axisym(mesh%x(:,mesh%fnode(:,i)))
      mesh%normal(:,i) = mesh%normal(:,i)*mesh%area(i)
    end do !i

  end subroutine mesh_axisymmetry_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Computes the swept volume obtained by rotating a polygon axisymmetrically. The
  ! relation from Rider and Kothe's paper is used to get the axisymmetric volume. However,
  ! note that this function computes the axisymmetric volume obtained after sweeping by a
  ! full 2*pi rad; as opposed to sweeping only by pi rad. Also, this function can handle
  ! clockwise or anti-clockwise ordering of the polygon nodes.
  ! See W. J. Rider, & D. B. Kothe (1998). Reconstructing volume tracking. Journal of
  ! computational physics, 141(2), 112-152.
  real(r8) function polygon_volume_axisym(x)

    use cell_geometry

    real(r8), intent(in) :: x(:,:)

    integer :: nnode, i
    real(r8), allocatable :: ccx(:,:)
    real(r8) :: pi, ax, ay, bx, by, cross

    nnode = size(x, dim=2)

    INSIST(nnode > 2)

    ! allocate an array to hold reordered nodes
    allocate(ccx(2,nnode+1))

    pi = 4.0_r8*atan(1.0_r8)

    ! check if the node-set is in counter-clockwise order. If not, reorder to be counter-clockwise.
    ax = x(1,2)-x(1,1)
    ay = x(2,2)-x(2,1)
    bx = x(1,3)-x(1,1)
    by = x(2,3)-x(2,1)

    cross = ax*by - ay*bx

    ! reorder if necessary
    if (cross < 0.0_r8) then
      ccx(:,1) = x(:,1)
      do i = 1, nnode-1
        ccx(:,i+1) = x(:,nnode-(i-1))
      end do !i
      ! the first node is repeated for the sake of simplicity in the volume calculation loop
      ccx(:,nnode+1) = x(:,1)
    else
      ccx(:,1:nnode) = x(:,:)
      ccx(:,nnode+1) = x(:,1)
    end if

    ! Rider and Kothe's formula assumes counter-clockwise ordering. The above reordering
    ! should have made the node-set order counter-clockwise at this point. If not, the
    ! following relation will give an incorrect axisymmetric volume. Also, this is the
    ! volume obtained by sweeping 2*pi rad, not pi rad.
    polygon_volume_axisym = 0.0_r8

    do i = 1, nnode
      polygon_volume_axisym = polygon_volume_axisym + &
        (ccx(1,i)+ccx(1,i+1)) * (ccx(1,i)*ccx(2,i+1)-ccx(1,i+1)*ccx(2,i))
    end do !i

    polygon_volume_axisym = pi/3.0_r8 * polygon_volume_axisym

  end function polygon_volume_axisym

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Computes the swept area obtained by rotating a face (line) axisymmetrically by 2*pi rad.
  real(r8) function face_area_axisym(x)

    real(r8), intent(in) :: x(:,:)

    real(r8) :: pi, dr, dz, rmin, k, l, incline

    pi = 4.0_r8*atan(1.0_r8)

    INSIST(size(x, dim=2) == 2)

    dr = x(1,1)-x(1,2)
    dz = x(2,1)-x(2,2)

    rmin = min(abs(x(1,1)), abs(x(1,2)))

    ! length of face
    l = sqrt(dr*dr + dz*dz)

    if (abs(dr) < 1e-14_r8) then
      incline = 0.0_r8
    else
      ! slope of face
      k = dz/dr
      incline = l/(2.0_r8*sqrt(1.0_r8+k*k))
    end if

    ! swept area
    face_area_axisym = 2.0_r8*pi*l*(rmin + incline)

  end function face_area_axisym

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module geom_axisymmetric
