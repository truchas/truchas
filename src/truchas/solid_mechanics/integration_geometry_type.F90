!!
!! This module provides a type which contains integration geometry, held
!! similarly to the unstr_mesh type. The integration geometry consists of a
!! collection of integration points (IPs) surrounding each mesh node. These
!! points are analagous to a face-center in the mesh, having an associated
!! location, area, and normal vector.
!!
!! This collection does not include boundary integration points, since those are
!! included as part of a special RHS calculation in the solid mechanics solver.
!!
!! Programming interface:
!!
!!   x - the rank-2 real array of IP coordinates; x(:,j) is the position in R^3
!!       of integration point j. Its shape is [3,npoint]
!!
!!   n - the rank-2 real array of oriented IP face areas; n(:,j) is the oriented
!!       area of IP j. Its shape is [3,npoint]
!!
!!   npoint, xnpoint - pair of rank-1 integer arrays storing the node-IP data:
!!       npoint(xnpoint(j):xnpoint(j+1)-1) is the unordered list of IP indices
!!       surrounding mesh node j. The shape of xnpoint is [nnode+1] and the
!!       shape of npoint is [xnpoint(nnode+1)-1].
!!
!!   nppar - an integer bit mask array storing the relative node-IP
!!        orientations: btest(nppar(j),k) is true when IP k of node j is inward
!!        oriented with respect to node j, and false when it is outward
!!        oriented. The shape of nnpar is [nnode].
!!
!! References:
!!
!! - Truchas Physics & Algorithms handbook.
!!
!! - C. Bailey and M. Cross. A finite volume procedure to solve elastic solid
!! mechanics problems in three dimensions on an unstructured mesh. International
!! Journal for Numerical Methods in Engineering, 38:1757–1776, 1995.
!!
!! - O. C. Zienkiewicz. The Finite Element Method. McGraw-Hill, New York, NY,
!! 1977.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! July 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module integration_geometry_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use integration_geometry_type
  implicit none
  private

  type, public :: integration_geometry
    type(unstr_mesh), pointer, private :: mesh => null() ! unowned reference

    integer :: npt
    real(r8), allocatable :: x(:,:), n(:,:), volume(:), subvolume(:)
    real(r8), allocatable, private :: xi(:,:), jacobian_inverse(:,:,:)
    integer, allocatable :: nppar(:)
    integer, allocatable :: npoint(:), xnpoint(:) ! node to IP connectivity
    integer, allocatable :: cpoint(:), xcpoint(:) ! cell to IP connectivity
    integer, allocatable :: pcell(:) ! IP to cell connectivity

    type(matrix_box), allocatable :: grad_shape(:)

    real(r8), private :: tet4_xi_ip(3,4), pyr5_xi_ip(3,5), wed6_xi_ip(3,6), hex8_xi_ip(3,8)
    real(r8), private :: tet4_grad_shape_l(3,4,4), pyr5_grad_shape_l(3,5,5), wed6_grad_shape_l(3,6,5), hex8_grad_shape_l(3,8,6)
  contains
    procedure :: init
    procedure :: compute_linear_grad
  end type integration_geometry


  type matrix_box
    real(r8), allocatable :: p(:,:)
  end type matrix_box


  real(r8), parameter :: tet4_xi_node(3,4) = reshape([], [3,4])
  real(r8), parameter :: pyr5_xi_node(3,5) = reshape([], [3,5])
  real(r8), parameter :: wed6_xi_node(3,6) = reshape([], [3,6])
  real(r8), parameter :: hex8_xi_node(3,8) = reshape([&
      & -1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8,  1.0_r8, -1.0_r8, &
      & -1.0_r8,  1.0_r8, -1.0_r8, &
      & -1.0_r8, -1.0_r8,  1.0_r8, &
      &  1.0_r8, -1.0_r8,  1.0_r8, &
      &  1.0_r8,  1.0_r8,  1.0_r8, &
      & -1.0_r8,  1.0_r8,  1.0_r8], [3,8])

  real(r8), parameter :: tet4_shape_coeff(4) = []
  real(r8), parameter :: pyr5_shape_coeff(5) = []
  real(r8), parameter :: wed6_shape_coeff(6) = []
  real(r8), parameter :: hex8_shape_coeff(8) = 1.0_r8 / 8.0_r8

contains

  subroutine init(this, mesh)

    class(integration_geometry), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh

    integer :: f, n, i, j
    integer, allocatable :: nface(:), xnface(:), k(:)

    this%npt = count(mesh%fcell(:,:this%nface_onP) > 0)
    allocate(this%x(3,this%npt), this%n(3,this%npt), this%nppar(this%npt), k(this%nnode_onP), &
        this%xnpoint(mesh%nnode_onP+1), this%volume(this%nnode_onP), this%subvolume(this%npt))

    ! Get the xnpoint offsets.
    ! Each node is associated with j internal integration points. There
    ! is one integration point for each pair (face,cell). Domain
    ! boundaries are not accociated with a cell center, and are skipped.
    call compute_nface(mesh, nface, xnface) ! TODO might be able to get rid of this using k(n) and a loop like the next one
    j = 1
    this%xnpoint(1) = 1
    do n = 1, mesh%nnode_onP
      j = 0
      do xf = xnface(n), xnface(n+1)-1
        f = nface(xf)
        j = j + count(mesh%fcell(:,f) > 0)
      end do
      this%xnpoint(n+1) = this%xnpoint(n) + j
    end do

    ! For each integration point, compute the position, area vector and
    ! node connectivity
    j = 1
    k = 0
    do f = 1, this%nface_onP
      do k = 1, 2
        i = mesh%fcell(k,f)
        if (i < 1) cycle
        this%x(:,j) = (this%face_centroid(:,f) + this%cell_centroid(:,i)) / 2
        this%n(:,j) =
        do xn = mesh%xfnode(f), mesh%xfnode(f+1)-1
          n = mesh%fnode(xn)
          if (n <= mesh%nnode_onP) then
            this%npoint(this%xnpoint(n)+k(n)) = j
            k(n) = k(n) + 1
          end if
        end do
        j = j + 1
      end do
    end do

    do n = 1, this%nnode_onP
      this%volume(n) = sum(this%subvolume(this%npoint(this%xnpoint(n):this%xnpoint(n+1)-1)))
    end do

  end subroutine init


  !! Compute the node-face connectivity from the mesh data.
  subroutine compute_nface(mesh, nface, xnface)

    type(unstr_mesh), intent(in) :: mesh
    integer, intent(out), allocatable :: nface(:), xnface(:)

    integer :: f, n, j, k(mesh%nnode_onP)

    n = 0
    do f = 1, mesh%nface
      n = n + count(mesh%fnode(xfnode(f):xfnode(f+1)-1) <= mesh%nnode_onP)
    end do
    allocate(xnface(mesh%nnode_onP+1), nface(n))

    ! Set up offsets (xnface)
    k = 0
    do f = 1, this%nface
      do xn = mesh%xfnode(f), mesh%xfnode(f+1)-1
        n = mesh%fnode(xn)
        if (n <= mesh%nnode_onP) k(n) = k(n) + 1
      end do
    end do

    xnface(1) = 1
    do n = 1, this%nnode_onP
      xnface(n+1) = xnface(n) + k(n)
    end do

    ! Write nface data
    j = 1
    k = 0
    do f = 1, mesh%nface
      do xn = mesh%xfnode(f), mesh%xfnode(f+1)-1
        n = mesh%fnode(xn)
        if (n <= mesh%nnode_onP) then
          nface(xnface(n) + k(n)) = f
          k(n) = k(n) + 1
        end if
      end do
    end do

  end subroutine compute_nface


  !! Compute the gradient of the shape function in global coordinates for each
  !! node of each cell at each of the integration points associated with that
  !! cell. This involves first computing the Jacobian as described on pg 1762 of
  !! Bailey & Cross 1995, and inverting it. The inverse Jacobian is multiplied
  !! against the shape function gradients in the reference coordinate system
  !! for the cell type.
  subroutine compute_grad_shape(this)

    class(integration_geometry), intent(inout) :: this

    integer :: j, d, p, xp
    real(r8) :: jacobian(3,3), jacobian_inverse(3,3)

    call this%compute_reference_quantities

    allocate(this%grad_shape(this%npt))

    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
          cp => this%cpoint(this%xcpoint(j):this%xcpoint(j+1)-1))
        do xp = 1, this%xcpoint(j+1)-this%xcpoint(j)
          p = this%xcpoint(j) + xp - 1
          do d = 1, 3
            !call compute_linear_grad(this%mesh%x(d,cn), p, jacobian(:,d), global=.false.)
            jacobian(:,d) = matmul(grad_shape_l(:,:,xp), this%mesh%x(d,cn))
          end do

          call invert_3x3(jacobian, jacobian_inverse)

          this%grad_shape(p)%p = matmul(jacobian_inverse, grad_shape_l(:,:,p))
        end do
      end associate
    end do

  end subroutine compute_grad_shape


  subroutine compute_reference_quantities(this)

    use cell_topology

    class(integration_geometry), intent(inout) :: this

    integer :: t

    ! compute the integration point coordinates for each cell type in reference coordinate-space
    call compute_xi_ip(tet4_xi_node, tet4_faces, tet4_xface, tet4_fsize, this%tet4_xi_ip)
    call compute_xi_ip(pyr5_xi_node, pyr5_faces, pyr5_xface, pyr5_fsize, this%pyr5_xi_ip)
    call compute_xi_ip(wed6_xi_node, wed6_faces, wed6_xface, wed6_fsize, this%wed6_xi_ip)
    call compute_xi_ip(hex8_xi_node, hex8_faces, hex8_xface, hex8_fsize, this%hex8_xi_ip)

    call compute_grad_shape_l(tet4_shape_coeff, tet4_xi_node, tet4_xi_ip, this%tet4_grad_shape_l)
    call compute_grad_shape_l(pyr5_shape_coeff, pyr5_xi_node, pyr5_xi_ip, this%pyr5_grad_shape_l)
    call compute_grad_shape_l(wed6_shape_coeff, wed6_xi_node, wed6_xi_ip, this%wed6_grad_shape_l)
    call compute_grad_shape_l(hex8_shape_coeff, hex8_xi_node, hex8_xi_ip, this%hex8_grad_shape_l)

  end subroutine compute_reference_quantities


  !! Returns the local coordinate xi of the integration point p of a cell with
  !! nnode nodes. p is expected to correspond to a local face ID. Since this
  !! occurs on a reference element in a local coordinate space, the
  !! physical-space coordinates of the actual cell are not needed.
  subroutine compute_xi_ip(xi_node, faces, xface, fsize, xi_ip)

    real(r8), intent(in) :: xi_node(:,:)
    integer, intent(in) :: faces(:), xface(:), fsize(:)
    real(r8), intent(out) :: xi_ip(:,:)

    integer :: p
    real(r8) :: xi_cc(3), xi_fc(3)

    xi_cc = sum(xi_node, dim=2) / size(xi_node, dim=2)
    do p = 1, size(fsize)
      xi_fc = sum(xi_node(:,faces(xface(p):xface(p+1)-1)), dim=2) / fsize(p)
      xi_ip(:,p) = (xi_cc + xi_fc) / 2
    end do

  end subroutine compute_xi_ip


  subroutine compute_grad_shape_l(coeff, xi_node, xi_ip, grad_shape_l)

    real(r8), intent(in) :: coeff(:), xi_node(:,:), xi_ip(:,:)
    real(r8), intent(out) :: grad_shape_l(:,:,:)

    integer :: p, n

    do p = 1, size(xi_ip, dim=2)
      do n = 1, size(xi_node, dim=2)
        associate (xi_n => xi_node(:,n), xi_p => xi_ip(:,p))
          grad_shape_l(1,n,p) = coeff(n) * xi_n(1) * (1 + xi_n(2) * xi_p(2)) * (1 + xi_n(3) * xi_p(3))
          grad_shape_l(2,n,p) = coeff(n) * xi_n(2) * (1 + xi_n(1) * xi_p(1)) * (1 + xi_n(3) * xi_p(3))
          grad_shape_l(3,n,p) = coeff(n) * xi_n(3) * (1 + xi_n(1) * xi_p(1)) * (1 + xi_n(2) * xi_p(2))
        end associate
      end do
    end do

  end subroutine compute_grad_shape_l

end module integration_geometry_type
