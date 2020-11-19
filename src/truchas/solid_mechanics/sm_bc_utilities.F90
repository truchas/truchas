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
  public :: compute_link_index_connectivity
  public :: compute_ip_normals
  public :: compute_node_normals
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
      xfini(fi+1) = xfini(fi) + count
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
    xfi = 0
    do fi = 1, nfi
      f = face_index(fi)
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

    ASSERT(all(fini > 0))
    ASSERT(all(fini <= mesh%nnode_onP))

  end subroutine compute_index_connectivity


  !! This routine is almost identical to above, except instead of operating on
  !! literal faces and nodes, it will operate on a given connectivity structure
  !! xfnode/fnode. This provides the generality needed to work on connections
  !! between link faces and link nodes. As a consequence, it does not check if
  !! nodes are on-process before adding them to the new structure.
  !!
  !! Note that it's important here that the node index be listed in the correct
  !! order: we are given a list of faces in a certain order, and each node
  !! is added to the list as we iterate through those faces, and only if that
  !! node hasn't been added before. In the sm_gap_contact_bc type, we expect
  !! this behavior so that we can make a semblance of node groups from face
  !! groups.
  subroutine compute_link_index_connectivity(fnode, xfnode, face_index, fini, xfini, node_index)

    integer, intent(in) :: fnode(:), xfnode(:), face_index(:)
    integer, intent(out), allocatable :: fini(:), xfini(:), node_index(:)

    integer :: f, n, fi, xn, xfi, nfi, count, nnode, stat
    integer, allocatable :: ni_(:)

    nnode = maxval(fnode)
    allocate(ni_(nnode), stat=stat)
    INSIST(stat == 0)

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
      do xn = xfnode(f), xfnode(f+1)-1
        n = fnode(xn)
        count = count + 1
      end do
      xfini(fi+1) = xfini(fi) + count
    end do

    ! count the unique nodes
    count = 0
    ni_ = 0
    do fi = 1, nfi
      f = face_index(fi)
      do xn = xfnode(f), xfnode(f+1)-1
        n = fnode(xn)
        if (ni_(n) /= 0) cycle
        count = count + 1
        ni_(n) = count
      end do
    end do

    ! compute the connectivity
    allocate(fini(xfini(nfi+1)-1), node_index(count))
    count = 0
    ni_ = 0
    xfi = 0
    do fi = 1, nfi
      f = face_index(fi)
      do xn = xfnode(f), xfnode(f+1)-1
        n = fnode(xn)
        if (ni_(n) == 0) then
          count = count + 1
          ni_(n) = count
          node_index(count) = n
        end if
        xfi = xfi + 1
        fini(xfi) = ni_(n)
      end do
    end do

    ASSERT(all(fini > 0))

  end subroutine compute_link_index_connectivity


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
      if (n > mesh%nnode_onP) cycle
      i = i + 1
      normal(:,i) = ic%normal_boundary(xf,j)
    end do

  end subroutine compute_face_ip_normals


  !! Compute the "node normals", evaluated from an area-weighted average of the
  !! surrounding integration point face normals. The routine takes a mapping
  !! between faces and nodes, with each face-node pair indicating an integration
  !! point.
  subroutine compute_node_normals(fini, xfini, normal_ip, normal_node)

    integer, intent(in) :: fini(:), xfini(:)
    real(r8), intent(in) :: normal_ip(:,:)
    real(r8), intent(out) :: normal_node(:,:)

    integer :: nface, fi, xni, ni

    nface = size(xfini)-1

    normal_node = 0
    do fi = 1, nface
      do xni = xfini(fi), xfini(fi+1)-1
        ni = fini(xni)
        normal_node(:,ni) = normal_node(:,ni) + normal_ip(:,xni)
      end do
    end do

  end subroutine compute_node_normals


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
