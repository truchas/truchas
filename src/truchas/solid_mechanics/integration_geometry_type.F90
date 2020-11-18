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
!!   fface - a rank-2 integer array storing cell-local face IDs for a given
!!       global face ID. For global face ID f corresponding to local face ID k
!!       of cell j, such that f = cface(xcface(j) - 1 + k), k is given by the
!!       fface array. fcell(1,f) gives the local face ID for the cell opposite
!!       the face normal direction, and fcell(2,f) gives the local face ID for
!!       the cell in the face normal direction.
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
!!   xpxn - a rank-1 integer array storing cell-local node IDs, ordered by
!!        node-local IP IDs. That is,
!!        mesh%cnode(mesh%xcnode(c)-1 + xpxn(xnpoint(n)-1+k)) = n
!!        For global cell ID c, global node ID n for a node associated with c,
!!        and k is a node-local IP ID.
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
!! August 2020
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
  implicit none
  private

  type matrix_box
    real(r8), allocatable :: p(:,:)
  end type matrix_box

  type, public :: integration_geometry
    integer :: npt
    real(r8), allocatable :: n(:,:), volume(:)
    type(matrix_box), allocatable :: grad_shape(:)
    integer, allocatable :: nppar(:)
    integer, allocatable :: npoint(:), xnpoint(:) ! node to IP connectivity
    integer, allocatable :: xcpoint(:), pcell(:) ! cell to IP connectivity
    integer, allocatable :: xpxn(:) ! cell-local node ID from node-local IP ID (see above)
    integer, allocatable :: fface(:,:) ! cell-local face IDs for a given global face

    integer :: nlnode
    integer, allocatable :: lnode(:,:)
    integer, allocatable :: xlflnode(:), lflnode(:) ! link-face to link-node connectivity

    type(unstr_mesh), pointer, private :: mesh => null() ! unowned reference
  contains
    procedure :: init
  end type integration_geometry

  real(r8), parameter :: tet4_xi_node(3,4) = reshape([&
      & -1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8, -1.0_r8, -1.0_r8, &
      &  0.0_r8,  1.0_r8, -1.0_r8, & ! coalesced hex nodes 3 & 4
      &  0.0_r8,  0.0_r8,  1.0_r8], [3,4]) ! coalesced hex nodes 5, 6, 7, & 8
  real(r8), parameter :: pyr5_xi_node(3,5) = reshape([&
      & -1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8,  1.0_r8, -1.0_r8, &
      & -1.0_r8,  1.0_r8, -1.0_r8, &
      &  0.0_r8,  0.0_r8,  1.0_r8], [3,5]) ! coalesced hex nodes 5, 6, 7, & 8
  real(r8), parameter :: wed6_xi_node(3,6) = reshape([&
      & -1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8, -1.0_r8, -1.0_r8, &
      &  0.0_r8,  1.0_r8, -1.0_r8, & ! coalesced hex nodes 3 & 4
      & -1.0_r8, -1.0_r8,  1.0_r8, &
      &  1.0_r8, -1.0_r8,  1.0_r8, &
      &  0.0_r8,  1.0_r8,  1.0_r8], [3,6]) ! coalesced hex nodes 7 & 8
  real(r8), parameter :: hex8_xi_node(3,8) = reshape([&
      & -1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8, -1.0_r8, -1.0_r8, &
      &  1.0_r8,  1.0_r8, -1.0_r8, &
      & -1.0_r8,  1.0_r8, -1.0_r8, &
      & -1.0_r8, -1.0_r8,  1.0_r8, &
      &  1.0_r8, -1.0_r8,  1.0_r8, &
      &  1.0_r8,  1.0_r8,  1.0_r8, &
      & -1.0_r8,  1.0_r8,  1.0_r8], [3,8])

  real(r8), parameter :: tet4_shape_coeff(4) = [0.125_r8, 0.125_r8, 0.25_r8, 0.5_r8]
  real(r8), parameter :: pyr5_shape_coeff(5) = [0.125_r8, 0.125_r8, 0.125_r8, 0.125_r8, 0.5_r8]
  real(r8), parameter :: wed6_shape_coeff(6) = [0.125_r8, 0.125_r8, 0.25_r8, 0.125_r8, 0.125_r8, 0.25_r8]
  real(r8), parameter :: hex8_shape_coeff(8) = 0.125_r8

  interface select_topology
    module procedure select_topology_r8xxx
  end interface select_topology

contains

  subroutine init(this, mesh)

    use integration_cell_type

    class(integration_geometry), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh

    integer :: j, xp, p, xn, n
    type(integration_cell) :: ic

    this%mesh => mesh

    call compute_connectivity(this)
    call compute_link_nodes(this)
    call compute_grad_shape(this)

    ! Accumulate the node volume and compute the control volume face areas
    ! associated with each integration point.
    allocate(this%n(3,this%npt), this%volume(this%mesh%nnode_onP))
    this%volume = 0
    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        call ic%init(mesh%x(:,cn))

        do xn = 1, size(cn)
          n = cn(xn)
          if (n > this%mesh%nnode_onP) cycle
          this%volume(n) = this%volume(n) + ic%subvolume(xn)
        end do

        do xp = 1, this%xcpoint(j+1)-this%xcpoint(j)
          p = xp + this%xcpoint(j) - 1
          this%n(:,p) = ic%normal(xp)
        end do
      end associate
    end do

  end subroutine init


  !! Each node is associated with j integration points. There is one integration
  !! point for each pair (edge,cell). Domain boundaries are not accociated with
  !! a cell center, and are skipped.
  !!
  !! On output, this routine will have computed npt, xcpoint, xnpoint, npoint,
  !! and nppar.
  subroutine compute_connectivity(this)

    use cell_topology, only: cell_edges

    class(integration_geometry), intent(inout) :: this

    integer :: j, n, p, xf, f, xp, node1, node2, k(this%mesh%nnode_onP)
    integer, pointer :: edges(:,:) => null()

    ! cell-IP connectivity
    allocate(this%xcpoint(this%mesh%ncell+1))
    this%xcpoint(1) = 1
    do j = 1, this%mesh%ncell
      edges => cell_edges(this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
      this%xcpoint(j+1) = this%xcpoint(j) + size(edges, dim=2)
    end do
    this%npt = this%xcpoint(this%mesh%ncell+1)-1
    allocate(this%pcell(this%npt))
    do j = 1, this%mesh%ncell
      this%pcell(this%xcpoint(j):this%xcpoint(j+1)-1) = j
    end do

    ! IP-node connectivity
    k = 0
    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        edges => cell_edges(cn)
        do p = this%xcpoint(j), this%xcpoint(j+1)-1
          xp = p - this%xcpoint(j) + 1
          node1 = cn(edges(1,xp))
          node2 = cn(edges(2,xp))
          if (node1 <= this%mesh%nnode_onP) k(node1) = k(node1) + 1
          if (node2 <= this%mesh%nnode_onP) k(node2) = k(node2) + 1
        end do
      end associate
    end do

    allocate(this%xnpoint(this%mesh%nnode_onP+1))
    this%xnpoint(1) = 1
    do n = 1, this%mesh%nnode_onP
      this%xnpoint(n+1) = this%xnpoint(n) + k(n)
    end do

    allocate(this%npoint(this%xnpoint(this%mesh%nnode_onP+1)-1), &
        this%xpxn(this%xnpoint(this%mesh%nnode_onP+1)-1), this%nppar(this%mesh%nnode_onP))
    this%npoint = -1
    this%nppar = 0
    k = 0
    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        edges => cell_edges(cn)
        do p = this%xcpoint(j), this%xcpoint(j+1)-1
          xp = p - this%xcpoint(j) + 1
          node1 = cn(edges(1,xp))
          node2 = cn(edges(2,xp))

          if (node1 <= this%mesh%nnode_onP) then
            this%npoint(this%xnpoint(node1)+k(node1)) = p
            this%xpxn(this%xnpoint(node1)+k(node1)) = edges(1,xp)
            k(node1) = k(node1) + 1
          end if
          if (node2 <= this%mesh%nnode_onP) then
            ! By convention, the normal is oriented
            ! toward the 2nd node of the associated edge.
            this%nppar(node2) = ibset(this%nppar(node2),k(node2)+1)
            this%npoint(this%xnpoint(node2)+k(node2)) = p
            this%xpxn(this%xnpoint(node2)+k(node2)) = edges(2,xp)
            k(node2) = k(node2) + 1
          end if
        end do
      end associate
    end do

    ! Face-cell local face connectivity
    allocate(this%fface(2,this%mesh%nface))
    this%fface = 0
    do j = 1, this%mesh%ncell
      associate(cf => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        do xf = 1, size(cf)
          f = cf(xf)
          n = merge(2, 1, btest(this%mesh%cfpar(j),xf))
          this%fface(n,f) = xf
        end do
      end associate
    end do

    ASSERT(all(this%npoint > 0))
    ASSERT(all(this%npoint <= this%npt))

  end subroutine compute_connectivity


  !! TODO-WARN: need xlflnode/lflnode
  !! Finds linked node pairs from linked face pairs. Places linked nodes into
  !! appropriate link set IDs, and ensures linked nodes aren't duplicated.
  subroutine compute_link_nodes(this)

    class(integration_geometry), intent(inout) :: this

    logical :: touched(this%mesh%nnode)
    integer :: count, l, n
    integer, allocatable :: lnode(:,:)

    count = 0
    touched = .false.
    do l = 1, this%mesh%nlink
      lnode = face_lnode(this%mesh%lface(:,l), this%mesh)
      do n = 1, size(lnode,dim=2)
        if (.not.touched(lnode(1,n))) then
          touched(lnode(1,n)) = .true.
          count = count + 1
        end if
      end do
    end do

    this%nlnode = count
    allocate(this%lnode(2,count))
    count = 0
    touched = .false.
    do l = 1, this%mesh%nlink
      lnode = face_lnode(this%mesh%lface(:,l), this%mesh)
      do n = 1, size(lnode,dim=2)
        if (.not.touched(lnode(1,n))) then
          touched(lnode(1,n)) = .true.
          count = count + 1
          this%lnode(:,count) = lnode(:,n)
        end if
      end do
    end do

  contains

    !! Compute the global IDs for linked nodes on either side of a given face
    !! link, for the entire face. This algorithm searches the 2 face for a node
    !! matching the position of the first node of the 1 face. Then it matches up
    !! the following nodes relying on the assumption that these faces are
    !! oriented in opposite directions.
    !!
    !! Example: If node 1 of face 1 is linked against node 1 of face 2, then
    !! node 2 of face 1 would be linked against the last node of face 2.
    function face_lnode(lface, mesh) result(lnode)

      integer, intent(in) :: lface(:)
      type(unstr_mesh), intent(in) :: mesh
      integer, allocatable :: lnode(:,:)

      integer :: i, j, n2, fsize
      real(r8) :: x1(3), x2(3)

      fsize = mesh%xfnode(lface(1)+1) - mesh%xfnode(lface(1))
      allocate(lnode(2,fsize))

      ! Find the match in the second face to the first node of the first face.
      ! Put the result in j.
      x1 = mesh%x(:,mesh%fnode(mesh%xfnode(lface(1))))
      do j = 1, fsize
        n2 = mesh%fnode(mesh%xfnode(lface(2))+j-1)
        ASSERT(n2 <= mesh%nnode)
        x2 = mesh%x(:,n2)
        if (all(x1 == x2)) exit
      end do
      ASSERT(j <= fsize)

      ! Increment j backwards as i goes forwards. This relies on the faces being
      ! oriented in opposite direction.
      do i = 1, fsize
        !j = modulo(i + offset - 1, fsize) + 1
        lnode(1,i) = mesh%fnode(mesh%xfnode(lface(1))+i-1)
        lnode(2,i) = mesh%fnode(mesh%xfnode(lface(2))+j-1)
        j = merge(j-1, fsize, j>1) 
        ASSERT(all(mesh%x(:,lnode(1,i)) == mesh%x(:,lnode(2,i))))
      end do

    end function face_lnode
    

  end subroutine compute_link_nodes


  !! Compute the gradient of the shape function in global coordinates for each
  !! node of each cell at each of the integration points associated with that
  !! cell. This involves first computing the Jacobian as described on pg 1762 of
  !! Bailey & Cross 1995, and inverting it. The inverse Jacobian is multiplied
  !! against the shape function gradients in the reference coordinate system
  !! for the cell type.
  subroutine compute_grad_shape(this)

    external dgesv ! LAPACK general system solve

    class(integration_geometry), intent(inout) :: this

    integer :: j, d, p, xp, ipiv(3), stat, nnode
    real(r8) :: jacobian(3,3)
    real(r8), target :: tet4_grad_shape_l(3,4,6), pyr5_grad_shape_l(3,5,8), &
        wed6_grad_shape_l(3,6,9), hex8_grad_shape_l(3,8,12)
    real(r8), pointer :: grad_shape_l(:,:,:) => null()

    call compute_reference_quantities

    allocate(this%grad_shape(this%npt))

    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        nnode = size(cn)
        grad_shape_l => select_topology(nnode, tet4_grad_shape_l, pyr5_grad_shape_l, &
            wed6_grad_shape_l, hex8_grad_shape_l)

        do xp = 1, this%xcpoint(j+1)-this%xcpoint(j)
          do d = 1, 3
            jacobian(:,d) = matmul(grad_shape_l(:,:,xp), this%mesh%x(d,cn))
          end do

          ! Compute G = J^-1 * L, where G are the shape function gradients
          ! in global coordinates and L are the shape function gradients in
          ! reference coordinates, and J is the Jacobian.
          p = xp + this%xcpoint(j) - 1
          this%grad_shape(p)%p = grad_shape_l(:,:,xp)
          call dgesv(3, nnode, jacobian, 3, ipiv, this%grad_shape(p)%p, 3, stat)
          ASSERT(stat == 0)
        end do
      end associate
    end do

#ifndef NDEBUG
    ! Ensure each integration point's matrix is allocated.
    do p = 1, this%npt
      if (.not.allocated(this%grad_shape(p)%p)) then
        print *, p
        ASSERT(.false.)
      end if
    end do
#endif

  contains

    !! For each cell type, compute the integration point coordinates and shape
    !! function gradients in reference coordinate-space. Outputs are *_xi_ip and
    !! *_grad_shape_l arrays
    subroutine compute_reference_quantities()

      use cell_topology

      real(r8) :: tet4_xi_ip(3,6), pyr5_xi_ip(3,8), wed6_xi_ip(3,9), hex8_xi_ip(3,12)

      call compute_xi_ip(tet4_xi_node, tet4_edges, tet4_xi_ip)
      call compute_xi_ip(pyr5_xi_node, pyr5_edges, pyr5_xi_ip)
      call compute_xi_ip(wed6_xi_node, wed6_edges, wed6_xi_ip)
      call compute_xi_ip(hex8_xi_node, hex8_edges, hex8_xi_ip)

      call compute_grad_shape_l(tet4_shape_coeff, tet4_xi_node, tet4_xi_ip, tet4_grad_shape_l)
      call compute_grad_shape_l(pyr5_shape_coeff, pyr5_xi_node, pyr5_xi_ip, pyr5_grad_shape_l)
      call compute_grad_shape_l(wed6_shape_coeff, wed6_xi_node, wed6_xi_ip, wed6_grad_shape_l)
      call compute_grad_shape_l(hex8_shape_coeff, hex8_xi_node, hex8_xi_ip, hex8_grad_shape_l)

    end subroutine compute_reference_quantities

  end subroutine compute_grad_shape


  !! Returns the local coordinates xi of the integration points of a cell. Since
  !! this occurs on a reference element in a local coordinate space, the
  !! physical-space coordinates of an actual cell are not needed.
  subroutine compute_xi_ip(xi_node, edges, xi_ip)

    real(r8), intent(in) :: xi_node(:,:)
    integer, intent(in) :: edges(:,:)
    real(r8), intent(out) :: xi_ip(:,:)

    integer :: p
    real(r8) :: xi_cc(3), xi_ec(3)

    xi_cc = sum(xi_node, dim=2) / size(xi_node, dim=2)
    do p = 1, size(edges, dim=2)
      !xi_fc = sum(xi_node(:,faces(xface(p):xface(p+1)-1)), dim=2) / fsize(p)
      xi_ec = sum(xi_node(:,edges(:,p)), dim=2) / 2
      xi_ip(:,p) = (xi_cc + xi_ec) / 2
    end do

  end subroutine compute_xi_ip


  !! Returns the gradients of the shape functions in reference coordinates. This
  !! is done for a list of integration points, for the shape functions
  !! associated with the provided list of node coordinates, and provided
  !! coefficients associated with those nodes.
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


  function select_topology_r8xxx(nnode, tet4, pyr5, wed6, hex8) result(s)
    integer, intent(in) :: nnode
    real(r8), intent(in), target :: tet4(:,:,:), pyr5(:,:,:), wed6(:,:,:), hex8(:,:,:)
    real(r8), pointer :: s(:,:,:)
    select case (nnode)
    case (4)
      s => tet4
    case (5)
      s => pyr5
    case (6)
      s => wed6
    case (8)
      s => hex8
    case default
      s => null()
      ASSERT(.false.)
    end select
  end function

end module integration_geometry_type
