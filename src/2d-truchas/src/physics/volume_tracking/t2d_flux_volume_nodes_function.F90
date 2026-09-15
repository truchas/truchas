!!
!! Constructs the polygon swept through a cell face during volume tracking.
!!
!! Aditya K. Pandare <apandare@lanl.gov>, January 2020
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flux_volume_nodes_function

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use t2d_cell_geom_vof_type, only: t2d_cell_geom
  use t2d_locate_plane_os_function, only: locate_plane_os
  use t2d_plane_type, only: t2d_plane, alpha
  use near_zero_function, only: near_zero
  implicit none
  private

  public :: flux_volume_nodes

contains

  !! Given the volume moving through a cell face, returns the ordered unique
  !! vertices of the swept polygon. DIST supplies the initial estimate of the
  !! distance traveled by the fluxing plane.
  subroutine flux_volume_nodes(face, cell, dist, flux_vol, cutoff, nodes, nnode, axisym)

    integer, intent(in) :: face
    type(t2d_cell_geom), intent(in) :: cell
    real(r8), intent(in) :: dist, flux_vol, cutoff
    real(r8), intent(out) :: nodes(:,:)
    integer, intent(out) :: nnode
    logical, intent(in) :: axisym

    integer, parameter :: max_iter = 10

    integer :: f, fc, i, on_point
    real(r8) :: xfc(2), rho, edge(2,2), xint(2)
    type(t2d_plane) :: plane

    xfc = 0.5_r8*(cell%node(:,face)+cell%node(:,mod(face,cell%nfc)+1))
    rho = dot_product(dist*cell%face_normal(:,face)-xfc, cell%face_normal(:,face))

    ! locate_plane_os finds the region behind the plane; negating the outward
    ! face normal selects the swept region immediately inside the cell.
    plane%normal = -cell%face_normal(:,face)
    plane%rho = rho
    call locate_plane_os(plane%normal, flux_vol/cell%volume, cell%volume, cell%node, &
      cutoff, max_iter, axisym, plane, guess=.true.)

    nnode = 2
    nodes(:,1) = cell%node(:,face)
    nodes(:,2) = cell%node(:,mod(face,cell%nfc)+1)
    do f = 1, cell%nfc-1
      fc = modulo(face+f-1, cell%nfc)+1
      edge(:,1) = cell%node(:,fc)
      edge(:,2) = cell%node(:,mod(fc,cell%nfc)+1)
      if (.not.plane%intersects(edge)) cycle

      call plane%intersection_point(xint, on_point, edge)
      if (any([(all(near_zero(xint-nodes(:,i), alpha)), i=1,nnode)])) cycle
      nnode = nnode+1
      INSIST(nnode <= size(nodes, dim=2))
      nodes(:,nnode) = xint
    end do

  end subroutine flux_volume_nodes

end module t2d_flux_volume_nodes_function
