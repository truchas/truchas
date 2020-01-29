!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module gradient_2d_cc_function

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_2d_mesh_type
  use index_partitioning
  implicit none
  private

  public :: gradient_2d_cc, gradient_rz_cc

contains

  ! 2d gradient of cell-centered scalar `x` evaluated at cell centers
  subroutine gradient_2d_cc(mesh, x, gx, gy, w_node0, w_node1, w_face)

    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: x(:)
    real(r8), intent(out) :: gx(:), gy(:), w_node0(:), w_node1(:), w_face(:)

    integer :: i, j, k

    w_face = 0.0_r8

    ! node average data stored in w_node0 workspace array
    call node_avg(mesh, x, w_node0, w_node1)
    call gather_boundary(mesh%node_ip, w_node0)

    do i = 1, mesh%nface
      ASSERT(size(mesh%fnode(1:2,i)) == 2)
      ! linear interpolation of vertex values to face centroid is arithmetic mean
      ! for triangle and quadrilateral faces
      w_face(i) = sum(w_node0(mesh%fnode(1:2,i)))/real(size(mesh%fnode(1:2,i)),r8)
    end do

    do i = 1, mesh%ncell_onP
      gx(i) = 0.0_r8
      gy(i) = 0.0_r8
      associate (cface => mesh%cface(mesh%cstart(i):mesh%cstart(i+1)-1))
        do j = 1, size(cface)
          k = cface(j)
          if (btest(mesh%cfpar(i),pos=j)) then ! true if normal points inward
            ! note the thae normal associate with the mesh object has already been scaled
            ! by the face area, so we do not do it again
            gx(i) = gx(i) - mesh%normal(1,k)*w_face(k)
            gy(i) = gy(i) - mesh%normal(2,k)*w_face(k)
          else
            gx(i) = gx(i) + mesh%normal(1,k)*w_face(k)
            gy(i) = gy(i) + mesh%normal(2,k)*w_face(k)
          end if
        end do
      end associate
      gx(i) = gx(i)/mesh%volume(i)
      gy(i) = gy(i)/mesh%volume(i)
    end do

  end subroutine gradient_2d_cc

  ! 2d-axisymmetric gradient of cell-centered scalar `x` evaluated at cell centers
  subroutine gradient_rz_cc(mesh, x, gx, gy, w_node0, w_node1, w_face)

    use cell_geometry, only: cell_volume

    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: x(:)
    real(r8), intent(out) :: gx(:), gy(:), w_node0(:), w_node1(:), w_face(:)

    integer :: i, j, k
    real(r8) :: xf(2,2), jacobian, vol

    w_face = 0.0_r8

    ! node average data stored in w_node0 workspace array
    call node_avg(mesh, x, w_node0, w_node1)
    call gather_boundary(mesh%node_ip, w_node0)

    do i = 1, mesh%nface
      ASSERT(size(mesh%fnode(1:2,i)) == 2)
      ! linear interpolation of vertex values to face centroid is arithmetic mean
      ! for triangle and quadrilateral faces
      w_face(i) = sum(w_node0(mesh%fnode(1:2,i)))/real(size(mesh%fnode(1:2,i)),r8)
    end do

    do i = 1, mesh%ncell_onP
      gx(i) = 0.0_r8
      gy(i) = 0.0_r8
      associate (cface => mesh%cface(mesh%cstart(i):mesh%cstart(i+1)-1))
        do j = 1, size(cface)
          k = cface(j)
          xf(:,1) = mesh%x(:,mesh%fnode(1,k))
          xf(:,2) = mesh%x(:,mesh%fnode(2,k))
          ! Jacobian of transformation to use the Gauss-divergence theorem
          jacobian = face_area(xf) / mesh%area(k)
          if (btest(mesh%cfpar(i),pos=j)) then ! true if normal points inward
            ! note the thae normal associate with the mesh object has already been scaled
            ! by the face area, so we do not do it again
            gx(i) = gx(i) - mesh%normal(1,k)*w_face(k)*jacobian
            gy(i) = gy(i) - mesh%normal(2,k)*w_face(k)*jacobian
          else
            gx(i) = gx(i) + mesh%normal(1,k)*w_face(k)*jacobian
            gy(i) = gy(i) + mesh%normal(2,k)*w_face(k)*jacobian
          end if
        end do
      end associate

      associate (cell_nodes => mesh%cnode(mesh%cstart(i):mesh%cstart(i+1)-1))
        vol = cell_volume(mesh%x(:,cell_nodes))
      end associate
      gx(i) = gx(i)/vol
      gy(i) = gy(i)/vol
    end do

  end subroutine gradient_rz_cc


  subroutine node_avg(mesh, x_cell, x_node, w_node)

    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: x_cell(:)
    real(r8), intent(out) :: x_node(:), w_node(:)

    integer :: i, j

    x_node = 0.0_r8
    w_node = 0.0_r8

    do i = 1, mesh%ncell
      associate (cnode => mesh%cnode(mesh%cstart(i):mesh%cstart(i+1)-1))
        do j = 1, size(cnode)
          x_node(cnode(j)) = x_node(cnode(j)) + mesh%volume(i)*x_cell(i)
          w_node(cnode(j)) = w_node(cnode(j)) + mesh%volume(i)
        end do
      end associate
    end do

    do i = 1, mesh%nnode_onP
      x_node(i) = x_node(i)/w_node(i)
    end do

  end subroutine node_avg

  real(r8) function face_area(x)
    use cell_geometry, only: vector_length
    real(r8), intent(in) :: x(2,2)
    real(r8) :: normal(2)
    normal(1) = x(2,2) - x(2,1)
    normal(2) = x(1,1) - x(1,2)
    face_area = vector_length(normal)
  end function face_area

end module gradient_2d_cc_function
