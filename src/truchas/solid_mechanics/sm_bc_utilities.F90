!!
!! SM_BC_UTILITIES
!!
!! This module collects several functions used for sm_*_bc_types.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_utilities

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  implicit none
  private

  public :: compute_index_connectivity
  public :: compute_ip_normals
  public :: rotation_matrix

contains

  !! Compute the connectivity between face-index and node-index, so that from a
  !! given face index we can directly reach the node indices attached to that
  !! face. Note these indices should be distinguished from faces and nodes
  !! themselves. What we are concerned about here is the connection between the
  !! fi and ni that go into bff%index(fi) and index(ni). This will let us
  !! directly access related elements in these separate compressed data
  !! structures, one compressed against faces (the user-supplied condition on
  !! face-sets) and the other on nodes (the discretization center) with
  !! integration points connecting the two. These structures will discard
  !! off-process nodes, as we do not need to reach them here.
  !!
  !! The structure produced here is the pair of rank-1 arrays, fini and xfini,
  !! analogous to cnode and xcnode. These will take a face index and return a
  !! node index.
  subroutine compute_index_connectivity(mesh, face_index, fini, xfini, node_index)

    type(unstr_mesh), intent(in) :: mesh
    integer, intent(in) :: face_index(:)
    integer, intent(out), allocatable :: fini(:), xfini(:), node_index(:)

    integer :: f, n, fi, xn, xfi, nfi, count
    integer :: ni_(mesh%nnode_onP)

    if (allocated(fini)) deallocate(fini)
    if (allocated(xfini)) deallocate(xfini)
    if (allocated(node_index)) deallocate(node_index)

    nfi = size(face_index) ! number of face indices

    ! compute the offsets by counting the valid nodes for each face
    allocate(xfini(nfi+1))
    xfini(1) = 1
    do fi = 1, nfi
      count = 0
      f = face_index(fi)
      do xn = mesh%xfnode(f), mesh%xfnode(f+1)-1
        n = mesh%fnode(xn)
        if (n <= mesh%nnode_onP) count = count + 1
      end do
      xfini(fi+1) = xfini(f) + count
    end do

    ! count the unique nodes
    count = 0
    ni_ = 0
    do fi = 1, nfi
      f = face_index(fi)
      do xn = mesh%xfnode(f), mesh%xfnode(f+1)-1
        n = mesh%fnode(xn)
        if (n > mesh%nnode_onP) cycle
        if (ni_(n) /= 0) cycle
        count = count + 1
        ni_(n) = count
      end do
    end do
    
    ! compute the connectivity
    allocate(fini(xfini(nfi+1)-1), node_index(count))
    count = 0
    ni_ = 0
    do fi = 1, nfi
      f = face_index(fi)
      xfi = 0
      do xn = mesh%xfnode(f), mesh%xfnode(f+1)-1
        n = mesh%fnode(xn)
        if (n > mesh%nnode_onP) cycle
        if (ni_(n) == 0) then
          count = count + 1
          ni_(n) = count
          node_index(count) = n
        end if
        xfi = xfi + 1
        fini(xfi) = ni_(n)
      end do
    end do

  end subroutine compute_index_connectivity


  !! Compute the oriented face areas for the integration points
  subroutine compute_ip_normals(face_index, xfini, mesh, ig, normal)

    use integration_geometry_type

    integer, intent(in) :: face_index(:), xfini(:)
    type(unstr_mesh), intent(in) :: mesh
    type(integration_geometry), intent(in) :: ig
    real(r8), intent(out) :: normal(:,:)

    integer :: fi, f, xn, xni
    real(r8) :: normal_ip(3,4)

    do fi = 1, size(face_index)
      f = face_index(fi)
      call compute_face_ip_normals(mesh, ig, f, normal_ip)
      do xni = xfini(fi), xfini(fi+1)-1
        xn = xni - xfini(fi) + 1
        normal(:,xni) = normal_ip(:,xn)
      end do
    end do

  end subroutine compute_ip_normals


  !! Compute the oriented face areas for each integration surface on a face.
  !! Each IP surface is associated with a node.
  subroutine compute_face_ip_normals(mesh, ig, f, normal)

    use integration_geometry_type
    use integration_cell_type

    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    integer, intent(in) :: f
    real(r8), intent(out) :: normal(:,:)

    integer :: i, j, n, xf
    integer, pointer :: cn(:) => null()
    type(integration_cell) :: ic

    xf = ig%fface(1,f) ! local face index
    j = mesh%fcell(1,f) ! cell index
    cn => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1)
    call ic%init(mesh%x(:,cn))

    ASSERT(size(normal,dim=2) >= ic%fsize(xf))
    i = 0
    do j = 1, ic%fsize(xf)
      n = mesh%fnode(mesh%xfnode(f)+j-1)
      if (n > mesh%ncell_onP) cycle
      i = i + 1
      normal(:,i) = ic%normal_boundary(xf,j)
    end do

  end subroutine compute_face_ip_normals


  !! Generate matrix which rotates such that normal is in the "z" direction.
  !! This routine does not assume that the input vector is normalized.
  function rotation_matrix(normal)

    use cell_geometry, only: cross_product, normalized

    real(r8), intent(in) :: normal(:)
    real(r8) :: rotation_matrix(3,3)

    real(r8), parameter :: zdir(3) = [0.0_r8, 0.0_r8, 1.0_r8]
    real(r8) :: a(3), n(3)

    n = normalized(normal)
    a = normalized(cross_product(zdir, n))

    if (norm2(a) == 0) then
      rotation_matrix = 0
      rotation_matrix(1,1) = 1
      rotation_matrix(2,2) = 1
      rotation_matrix(3,3) = 1
    else
      rotation_matrix(1,:) = cross_product(a, n)
      rotation_matrix(2,:) = a
      rotation_matrix(3,:) = n
    end if

  end function rotation_matrix

end module sm_bc_utilities
