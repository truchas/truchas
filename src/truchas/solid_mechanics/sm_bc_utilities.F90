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
  public :: compute_inverted_connectivity
  public :: compute_ip_normals
  public :: compute_node_normals
  public :: rotation_matrix
  public :: contact_factor, derivative_contact_factor
  public :: compute_gradient_node_to_cell
  public :: check_if_matching_node
  public :: compute_stress, von_mises_stress, compute_deviatoric_stress

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
    integer :: ni_(mesh%nnode)

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
        !if (n <= mesh%nnode_onP) count = count + 1
        count = count + 1
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
        !if (n > mesh%nnode_onP) cycle
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
        !if (n > mesh%nnode_onP) cycle
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
    !ASSERT(all(fini <= mesh%nnode_onP))

  end subroutine compute_index_connectivity


  ! subroutine compute_inverted_connectivity(nnode, fini, xfini, nifi, xnifi)

  !   integer, intent(in) :: nnode, fini(:), xfini(:)
  !   integer, intent(out), allocatable :: nifi(:), xnifi(:)



  !   nface = 0
  !   do fi = 1, size(xfini)-1
  !     do xni = xfini(fi), xfini(fi+1)-1
  !       ni = fini(xni)
  !       nface(ni) = nface(ni) + 1
  !     end do
  !   end do

  !   allocate(xnifi(nnode+1))
  !   xnifi(1) = 1
  !   do ni = 1, nnode
  !     xnifi(ni+1) = xnifi(ni) + nface(ni)
  !   end do

  !   allocate(nifi(xnifi(nnode+1)-1))
  !   nface = 0
  !   do fi = 1, size(xfini)-1
  !     do xni = xfini(fi), xfini(fi+1)-1
  !       ni = fini(xni)
  !       nifi(xnifi(ni)+nface(ni)) = fi
  !       nface(ni) = nface(ni) + 1
  !     end do
  !   end do

  ! end subroutine compute_inverted_connectivity

  !! Given a mapping from a to b (ab and xab), compute a mapping
  !! from b to a (ba and xba). It is assumed that both a and b start
  !! at 1. The maximum value of b (nb) must be provided.
  !!
  !! Example usage:
  !!   call compute_inverted_connectivity(mesh%fnode, mesh%xfnode, mesh%nnode, nface, xnface)
  subroutine compute_inverted_connectivity(ab, xab, nb, ba, xba)

    integer, intent(in) :: ab(:), xab(:), nb
    integer, intent(out), allocatable :: ba(:), xba(:)

    integer :: a, b, x, na, j(nb)

    na = size(xab)-1

    ! compute offsets
    allocate(xba(nb+1))
    xba(1) = 1
    j = 0
    do a = 1, na
      do x = xab(a), xab(a+1)
        b = ab(x)
        j(b) = j(b) + 1
      end do
    end do
    do b = 1, nb
      xba(b+1) = xba(b) + j(b)
    end do

    ! compute inverse map
    allocate(ba(xba(size(xba))))
    j = 0
    do a = 1, na
      do x = xab(a), xab(a+1)
        b = ab(x)
        ba(xba(b)+j(b)) = a
        j(b) = j(b) + 1
      end do
    end do

  end subroutine compute_inverted_connectivity


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
      !if (n > mesh%nnode_onP) cycle
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


  ! !! For each node, compute a collection of normals. Each node is surrounded by
  ! !! a number of faces, each with a potentially different normal. These faces
  ! !! each belong to some face set. For each node, we collect the surrounding
  ! !! faces into one group for each face set. Faces belonging
  ! subroutine compute_node_multi_normals(mesh, face, fini, xfini, nnode, normal_ip, normal_node)

  !   type(unstr_mesh), intent(in) :: mesh
  !   integer, intent(in) :: face(:), fini(:), xfini(:), nnode
  !   real(r8), intent(in) :: normal_ip(:,:)
  !   real(r8), intent(out), allocatable :: normal_node(:,:)

  !   integer :: nface, nset, f, fi, xni, ni, j(nnode)
  !   logical :: touched(nnode)

  !   nface = size(xfini)-1
  !   nset = size(mesh%face_set_id)
  !   j = 0
  !   touched = .false.
  !   do k = 1, nset
  !     do fi = 1, nface
  !       f = face(fi)
  !       if (btest(mesh%face_set_mask(f),k)) then
  !         do xni = xfini(fi), xfini(fi+1)-1
  !           ni = fini(xni)
  !           if (.not.touched(ni)) then
  !             j(fi) = j(fi) + 1
  !             touched(ni) = .true.
  !           end if
  !         end do
  !       end if
  !     end do
  !   end do
  !   j = -1
  !   touched = .false.
  !   do k = 1, nset
  !     do fi = 1, nface
  !       f = face(fi)
  !       if (btest(mesh%face_set_mask(f),k)) then
  !         do xni = xfini(fi), xfini(fi+1)-1
  !           ni = fini(xni)
  !           if (.not.touched(ni)) then
  !             j(fi) = j(fi) + 1
  !             touched(ni) = .true.
  !           end if
  !           normal_node(:,ni+j(ni)) = normal_node(:,ni+j(ni)) + normal_ip(:,xni)
  !         end do
  !       end if
  !     end do
  !   end do


  ! !   nface = size(xfini)-1
  ! !   do fi = 1, nface
  ! !     f = face(fi)
  ! !     face_set = face_set_id(mesh, f)

  ! !     do xni = xfini(fi), xfini(fi+1)-1
  ! !       ni = fini(xni)
  ! !       normal_node(:,ni) = normal_node(:,ni) + normal_ip(:,xni)
  ! !     end do
  ! !   end do

  ! ! contains

  ! !   pure function face_set_id(mesh, f) result(k)
  ! !     type(unstr_mesh), intent(in) :: mesh
  ! !     integer, intent(in) :: f
  ! !     integer :: k
  ! !     integer :: nsets
  ! !     nsets = size(mesh%face_set_id)
  ! !     do k = 1, nsets
  ! !       if (btest(mesh%face_set_mask(f),k)) exit
  ! !     end do
  ! !     INSIST(k <= nsets)
  ! !   end function face_set_id

  ! end subroutine compute_node_multi_normals


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


  !! Contact helper functions
  ! !! Compute the contact residual term
  ! pure subroutine compute_contact(penalty, distance, traction, &
  !     displ, ftot, stress_factor, residual)

  !   real(r8), intent(in) :: penalty, distance, traction
  !   real(r8), intent(in) :: displ(:,:), ftot(:,:), stress_factor(:), normal(:)
  !   real(r8), intent(out) :: residual(:,:)

  !   real(r8) :: stress1, stress2, x1, x2, s, tn, l, v(2)

  !   stress1 = dot_product(normal, ftot(:,1))
  !   stress2 = dot_product(normal, ftot(:,2))
  !   x1 = dot_product(normal, displ(:,1))
  !   x2 = dot_product(normal, displ(:,2))

  !   s = x2 - x1
  !   tn = - stress1 / this%area(i)
  !   l = contact_factor(s, tn, distance, traction)

  !   v(1) = stress2 + penalty*(x2 - x1) * stress_factor(1)
  !   v(2) = stress1 + penalty*(x1 - x2) * stress_factor(2)

  !   residual(:,1) = normal * l * v(1)
  !   residual(:,2) = normal * l * v(2)

  ! end subroutine compute_contact


  !! Given the difference in normal displacements s, and the tensile force normal
  !! to the surface tn, compute the contact factor.
  !!
  !! If s <= 0, then the nodes are in contact or inside one another.
  !! If tn <= 0, the normal traction is compressive
  real(r8) function contact_factor(s, tn, distance, traction)

    real(r8), intent(in) :: s, tn, distance, traction

    real(r8) :: ls, lt, x

    if (s <= 0) then
      ls = 1
    else if (s >= distance) then
      ls = 0
    else
      x = s / distance - 1
      ls = x**2 * (2*x + 3)
    end if

    if (tn <= 0) then
      lt = 1
    else if (tn >= traction) then
      lt = 0
    else
      x = tn / traction - 1
      lt = x**2 * (2*x + 3)
    end if

    contact_factor = ls * lt
    ASSERT(contact_factor >= 0 .and. contact_factor <= 1)

  end function contact_factor


  pure function derivative_contact_factor(s, tn, distance, traction) result(dl)

    real(r8), intent(in) :: s, tn, distance, traction
    real(r8) :: dl(2)

    real(r8) :: ls, lt, x

    dl = 0

    if (s <= 0) then
      ls = 1
    else if (s >= distance) then
      ls = 0
    else
      x = s / distance - 1
      ls = x**2 * (2*x + 3)
      dl(1) = 6 * x * (x + 1) / distance
    end if

    if (tn <= 0) then
      lt = 1
    else if (tn >= traction) then
      lt = 0
    else
      x = tn / traction - 1
      lt = x**2 * (2*x + 3)
      dl(2) = 6 * x * (x + 1) / traction
    end if

    dl(1) = dl(1) * lt
    dl(2) = dl(2) * ls

  end function derivative_contact_factor


  !! Calculate the cell-centered gradient {dQ_dx, dQ_dy, dQ_dz} of the
  !! vertex-centered displacement component Q by averaging the gradient,
  !! defined by the solution of a 3x3 system of equations, over the
  !! volume. This averaging procedure is accomplished through an
  !! integral over the cell volume. The gradient at any logical point
  !! within the cell is a ratio of Jacobian-like determinants.
  subroutine compute_gradient_node_to_cell(x, volume, f, g)

    real(r8), intent(in) :: x(:,:), volume, f(:) ! node-centered
    real(r8), intent(out) :: g(:) ! cell-centered

    integer :: d
    real(r8) :: phi(3,size(f))

    ASSERT(size(f) == size(x,dim=2))
    do d = 1, 3
      phi = x
      phi(d,:) = f
      g(d) = determinant_vol_avg(phi) / volume
    end do

  contains

    !! Compute the volume-averaged value of a Jacobian determinant |J|:
    !! Avg = Integral(|J|dV), where J is given by column vectors X, Y,
    !! and Z:             -                       -
    !!                    |  X_xi   Y_xi   Z_xi   |
    !!              J =   |  X_eta  Y_eta  Z_eta  |
    !!                    |  X_zeta Y_zeta Z_zeta |
    !!                    -                       -
    !! where X_xi == dX/dxi, X_eta == dX/deta, X_zeta == dX/dzeta, etc.
    !! The volume integral is converted to a surface integral and
    !! evaluated following J. Dukowicz, JCP 74: 493-496 (1988).
    !!
    !! Note 1: There is a minus sign introduced here because the formula
    !!         presented by Dukowicz orients faces in the opposite
    !!         direction.
    function determinant_vol_avg(phi) result(avg)

      use cell_geometry, only: triple_product

      real(r8), intent(in) :: phi(:,:)
      real(r8) :: avg

      integer, pointer :: faces(:) => null(), xface(:) => null()
      integer :: f
      real(r8) :: phif(3,4), x1(3), x2(3), x3(3)

      call get_cell_faces(size(phi,dim=2), faces, xface)

      avg = 0
      do f = 1, size(xface)-1
        associate (node => faces(xface(f):xface(f+1)-1))
          phif(:,:size(node)) = phi(:,node)
          if (size(node) == 3) phif(:,4) = phif(:,3)
          x1 = phif(:,2) + phif(:,3)
          x2 = phif(:,1) + phif(:,2)
          x3 = phif(:,3) + phif(:,4)
          avg = avg - triple_product(x1, x2, x3) ! See note 1.
        end associate
      end do
      avg = avg / 12

    end function determinant_vol_avg

  end subroutine compute_gradient_node_to_cell


  subroutine get_cell_faces(nnode, faces, xface)
    use cell_topology
    integer, intent(in) :: nnode
    integer, pointer :: faces(:), xface(:)
    select case (nnode)
    case (4)  ! tet
      faces => TET4_FACES
      xface => TET4_XFACE
    case (5)  ! pyramid
      faces => PYR5_FACES
      xface => PYR5_XFACE
    case (6)  ! wface
      faces => WED6_FACES
      xface => WED6_XFACE
    case (8)  ! hex
      faces => HEX8_FACES
      xface => HEX8_XFACE
    case default
      faces => null()
      xface => null()
    end select
  end subroutine get_cell_faces


  !! This routine is used by the various sm_bc_cXdY types, to identify
  !! nodes for each individual combination used by those types. This
  !! routine expects icontact and idispl arrays with lengths
  !! corresponding to the requested type. E.g., icontact of size 0 and
  !! idispl of size 2 for BC-combos with no contact, and two
  !! displacements. On output, the routine will provide a logical value
  !! indicating whether the node matches the requested type, and indexes
  !! for those BCs into sm_bc_node_list::bcid.
  subroutine check_if_matching_node(ni, n, bcid, xbcid, nnode_onP, xcontact, &
      icontact, idispl, matching_node)

    integer, intent(in) :: n, ni, bcid(:), xbcid(:), nnode_onP, xcontact
    integer, intent(out) :: icontact(:), idispl(:)
    logical, intent(out) :: matching_node

    integer :: b, nbc, ibc, xibc

    matching_node = .false.
    icontact = 0
    idispl = 0

    if (n > nnode_onP) return ! only consider owned nodes

    ! Count the BCs applied to this node. If the node doesn't have
    ! exactly the requested number of BCs, it is disqualified.
    nbc = xbcid(ni+1) - xbcid(ni)
    if (nbc /= size(icontact)+size(idispl)) return

    do xibc = xbcid(ni), xbcid(ni+1)-1
      ibc = bcid(xibc)
      if (ibc >= xcontact) then
        ! append to the list, or exit if we've found too many
        do b = 1, size(icontact)
          if (icontact(b) == 0) then
            icontact(b) = xibc
            exit
          end if
        end do
        if (b > size(icontact)) return

      else
        ! append to the list, or exit if we've found too many
        do b = 1, size(idispl)
          if (idispl(b) == 0) then
            idispl(b) = xibc
            exit
          end if
        end do
        if (b > size(idispl)) return
      end if
    end do

    ! If all requested conditions have been set, and no other conditions
    ! have triggered an early exit, we have a matching node.
    if (all(icontact /= 0) .and. all(idispl /= 0)) matching_node = .true.

  end subroutine check_if_matching_node


  !! Computes the stress from the elastic strain and Lame parameters.
  pure subroutine compute_stress(lame1, lame2, strain, stress)
    real(r8), intent(in) :: lame1, lame2, strain(:)
    real(r8), intent(out) :: stress(:)
    ! ASSERT(size(strain) == 6)
    ! ASSERT(size(stress) == 6)
    stress = 2 * lame2 * strain
    stress(:3) = stress(:3) + lame1 * sum(strain(:3))
  end subroutine compute_stress


  real(r8) pure function von_mises_stress(stress)
    real(r8), intent(in) :: stress(:)
    !ASSERT(size(stress) == 6)
    von_mises_stress = sqrt(((stress(1) - stress(2))**2 &
        &                  + (stress(2) - stress(3))**2 &
        &                  + (stress(3) - stress(1))**2 &
        &                  + 6 * sum(stress(4:6)**2)) / 2)
  end function von_mises_stress


  pure subroutine compute_deviatoric_stress(stress, deviatoric_stress)
    real(r8), intent(in) :: stress(:)
    real(r8), intent(out) :: deviatoric_stress(:)
    real(r8) :: mean_stress
    ! ASSERT(size(stress) == 6)
    ! ASSERT(size(deviatoric_stress) == 6)
    mean_stress = sum(stress(1:3)) / 3
    deviatoric_stress = stress
    deviatoric_stress(1:3) = deviatoric_stress(1:3) - mean_stress
  end subroutine compute_deviatoric_stress

end module sm_bc_utilities
