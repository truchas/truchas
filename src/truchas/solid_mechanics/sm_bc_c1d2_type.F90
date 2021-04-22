!!
!! SM_BC_NODE_CONTACT_TYPES
!!
!! TODO
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_c1d2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication, only: global_sum
  use truchas_logging_services
  use scalar_func_containers, only: scalar_func_ptr
  use unstr_mesh_type
  use sm_bc_list_type
  use sm_bc_node_list_type
  use sm_bc_utilities, only: contact_factor, derivative_contact_factor
  implicit none
  private

  !! One contact and two displacements
  type, public :: sm_bc_c1d2
    private
    integer, allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:,:), dvalue(:,:), tangent(:,:)
    logical, public :: enabled

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal_d(:,:,:), normal_gap(:,:), alpha(:)
    type(scalar_func_ptr), allocatable :: displf(:,:)
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_deriv
  end type sm_bc_c1d2

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized

    class(sm_bc_c1d2), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, li, lni, ni, icontact, idispl1, idispl2
    logical :: matching_node

    this%mesh => mesh
    this%penalty = penalty
    this%distance = distance
    this%normal_traction = traction

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      do lni = 1, 2
        call check_if_c1d2_node(nodebc%link(lni,li), icontact, idispl1, idispl2, matching_node)
        if (matching_node) nnode = nnode + 1
      end do
    end do
    allocate(this%index(nnode), this%value(3,nnode), this%dvalue(3,nnode))
    allocate(this%normal_d(3,2,nnode), this%normal_gap(3,nnode), this%tangent(3,nnode))
    allocate(this%linked_node(nnode), this%area(nnode), this%displf(2,nnode), this%alpha(nnode))

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      do lni = 1, 2
        ni = nodebc%link(lni,li)
        call check_if_c1d2_node(ni, icontact, idispl1, idispl2, matching_node)
        if (.not.matching_node) cycle
        nnode = nnode + 1
        this%index(nnode) = nodebc%node(ni)
        this%linked_node(nnode) = nodebc%node(nodebc%link(merge(1,2,lni==2),li))

        this%displf(1,nnode)%f => bc%displacement(nodebc%bcid(idispl1))%f
        this%displf(2,nnode)%f => bc%displacement(nodebc%bcid(idispl2))%f
        this%normal_d(:,1,nnode)  = nodebc%normal(:,idispl1) / norm2(nodebc%normal(:,idispl1))
        this%normal_d(:,2,nnode)  = nodebc%normal(:,idispl2) / norm2(nodebc%normal(:,idispl2))
        this%normal_gap(:,nnode) = nodebc%normal(:,icontact) / norm2(nodebc%normal(:,icontact))
        this%area(nnode) = norm2(nodebc%normal(:,icontact))

        this%tangent(:,nnode) = normalized(cross_product(this%normal_d(:,1,nnode), this%normal_d(:,2,nnode)))
        this%alpha(nnode) = dot_product(this%tangent(:,nnode), this%normal_gap(:,nnode))**2
      end do
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C1D2 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  contains

    !! Check if the given link index refers to a one contact + one displacement
    !! link pair. Here both sides of the link will have some
    subroutine check_if_c1d2_node(ni, icontact, idispl1, idispl2, matching_node)

      integer, intent(in) :: ni
      integer, intent(out) :: icontact, idispl1, idispl2
      logical, intent(out) :: matching_node

      integer :: nbc, bcid, xbcid

      matching_node = .false.
      icontact = 0
      idispl1 = 0
      idispl2 = 0

      ! Count the BCs applied to the two nodes on either side of the link.
      ! One must have 2 BCs, the other must have
      nbc = nodebc%xbcid(ni+1) - nodebc%xbcid(ni)
      if (nbc /= 3) return ! if either node has >2 or <1 BCs, the link is disqualified
      do xbcid = nodebc%xbcid(ni), nodebc%xbcid(ni+1)-1
        bcid = nodebc%bcid(xbcid)
        if (bcid > bc%dsz) then
          ! Assign exactly one contact condition.
          ! if more are found, disqualify this node.
          if (icontact == 0) then
            icontact = xbcid
          else
            return
          end if
        else
          if (idispl1 == 0) then
            idispl1 = xbcid
          else if (idispl2 == 0) then
            if (bcid /= nodebc%bcid(idispl1)) idispl2 = xbcid
          else
            return
          end if
        end if
      end do

      ! If all three conditions have been set, and no other conditions
      ! have triggered an early exit, we have a matching node.
      if (icontact /= 0 .and. idispl1 /= 0 .and. idispl2 /= 0) matching_node = .true.

    end subroutine check_if_c1d2_node

  end subroutine init


  subroutine compute(this, time, displ, ftot, stress_factor)

    class(sm_bc_c1d2), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)

    integer :: i, n1, n2, k
    real(r8) :: stress1, stress2, x1, x2, s, tn, l, v
    real(r8) :: args(0:3), displbc(2), x(3)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(i)

      ! displacement part
      args(1:) = this%mesh%x(:,n1)
      displbc(1) = this%displf(1,i)%eval(args) ! associated with this%normal(:,1,i)
      displbc(2) = this%displf(2,i)%eval(args) ! associated with this%normal(:,2,i)
      x = displ(:,n1) - displacement_vector(this%normal_d(:,:,i), displbc)
      x = x - dot_product(x, this%tangent(:,i)) * this%tangent(:,i)
      this%value(:,i) = this%penalty * stress_factor(n1) * x

      ! contact part
      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      s = x2 - x1
      tn = - stress1 / this%area(i) ! TODO-WARN: is the sign right?
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      stress1 = dot_product(this%tangent(:,i), ftot(:,n1))
      stress2 = dot_product(this%tangent(:,i), ftot(:,n2))
      x1 = dot_product(this%tangent(:,i), displ(:,n1))
      x2 = dot_product(this%tangent(:,i), displ(:,n2))
      v = stress2 + this%penalty * this%alpha(i) * (x2 - x1) * stress_factor(n1)

      this%value(:,i) = this%value(:,i) + this%tangent(:,i) * l * v
    end do

  contains

    !! in (d1,d2) displacements in the given normal directions,
    !! out a displacement vector a = b1*n1 + b2*n2 such that
    !!   dot(n1, a) = d1
    !!   dot(n2, a) = d2
    function displacement_vector(normal, displ) result(a)
      real(r8), intent(in) :: normal(:,:), displ(:)
      real(r8) :: a(3)
      real(r8) :: matrix(2,2), displacements(2)
      matrix(1,1) = 1
      matrix(2,1) = dot_product(normal(:,1),normal(:,2))
      matrix(2,2) = 1
      matrix(1,2) = matrix(2,1)
      displacements = displ
      call solve2x2(matrix, displacements)
      a = matmul(normal, displacements)
    end function displacement_vector

    !! Solve 2x2 a*x=b linear system.
    subroutine solve2x2(a, x)

      real(r8), intent(inout) :: a(:,:) ! in a out inv(a)
      real(r8), intent(inout) :: x(:) ! in b out x

      real(r8) :: det, s

      ASSERT(size(x) == 2)
      ASSERT(size(a,dim=1) == 2 .and. size(a,dim=2) == 2)

      det = 1 - a(1,2)*a(2,1)
      ASSERT(det /= 0)

      s = a(1,1)
      a(1,1) = a(2,2)
      a(2,2) = s

      a(1,2) = -a(1,2)
      a(2,1) = -a(2,1)

      a = a/det

      x = matmul(a, x)

    end subroutine solve2x2

  end subroutine compute


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine compute_deriv(this, time, displ, ftot, stress_factor, diag, F)

    class(sm_bc_c1d2), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), diag(:,:), F(:,:,:)

    integer :: i, nl, n, d
    real(r8) :: x(3)

    do i = 1, size(this%index)
      n = this%index(i)
      if (n > this%mesh%nnode_onP) cycle
      do d = 1,3
        x(d) = dot_product(this%tangent(:,i), F(:,d,n))
      end do
      this%dvalue(:,i) = this%tangent(:,i) * x &
          - this%penalty * stress_factor(n) * (1 - this%tangent(:,i)**2)
      ! diag(:,n) = diag(:,n) - this%normal(:,i) * x
      ! !diag(:,n) = diag(:,n) - diag(:,n) * normal(:,i)**2
      ! diag(:,n) = diag(:,n) - this%penalty * stress_factor(n) * this%normal(:,i)**2
    end do

  end subroutine compute_deriv

end module sm_bc_c1d2_type
