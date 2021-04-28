!!
!! SM_BC_NODE_CONTACT_TYPES
!!
!! TODO
!!
!! Zach Jibben <zjibben@lanl.gov>
!! March 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_node_contact_types

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

  !! One contact
  type, public :: sm_bc_c1
    private
    integer, allocatable, public :: index(:,:)
    real(r8), allocatable, public :: value(:,:,:), dvalue(:,:,:)
    logical, public :: enabled

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal(:,:)
  contains
    procedure :: init => c1_init
    procedure :: compute => c1_compute
    procedure :: compute_deriv => c1_compute_deriv
  end type sm_bc_c1

  !! One contact and one displacement
  type, public :: sm_bc_c1d1
    private
    integer, allocatable, public :: index(:,:)
    real(r8), allocatable, public :: value(:,:,:), dvalue(:,:,:)
    logical, public :: enabled

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal(:,:), normal_gap(:,:), align(:,:), alpha(:)
    type(scalar_func_ptr), allocatable :: displf(:)
  contains
    procedure :: init => c1d1_init
    procedure :: compute => c1d1_compute
    procedure :: compute_deriv => c1d1_compute_deriv
  end type sm_bc_c1d1

contains

  !! Initialize the single-normal displacement BC nodes.
  !!
  !! The main difficulty here, and in the other below routines for initializing
  !! node-BCs is identifying the appropriate nodes. Here we take a
  !! face_displacement object which holds BCs on faces, along with some helpful
  !! mapping arrays for looking at nodes connected to those faces, and quickly
  !! searching BC faces neighboring each given node. Thus one private procedure
  !! uses faces surrounding a given node to decide whether that node meets the
  !! conditions for this type.
  !!
  !! Here we're looking for nodes touching only faces with a displacement-BC in
  !! a single normal direction. That normal direction may change smoothly over
  !! a curved surface, so there is a small tolerance for deciding what qualifies
  !! as 'distinct' normal vectors between faces. Currently the tolerance is set
  !! to about 25 degrees.
  subroutine c1_init(this, mesh, nodebc, bc, penalty, distance, traction)

    !use integration_geometry_type

    class(sm_bc_c1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    !type(integration_geometry), intent(in) :: ig
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, li, ni, xbcid

    this%mesh => mesh
    this%penalty = penalty
    this%distance = distance
    this%normal_traction = traction

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      if (is_1contact_node(li)) nnode = nnode + 1
    end do
    allocate(this%index(2,nnode), this%value(3,2,nnode), this%dvalue(3,2,nnode), &
        this%normal(3,nnode), this%area(nnode))

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      if (.not.is_1contact_node(li)) cycle
      nnode = nnode + 1
      this%index(:,nnode) = nodebc%node(nodebc%link(:,li))
      ni = nodebc%link(1,li)
      xbcid = nodebc%xbcid(ni)
      this%area(nnode) = norm2(nodebc%normal(:,xbcid))
      this%normal(:,nnode) = nodebc%normal(:,xbcid) / norm2(nodebc%normal(:,xbcid))
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C1 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  contains

    !! To qualify as a 1-contact node, the nodes on both sides of the link may
    !! only have a single BC: a contact. If either has another contact or a
    !! displacement BC, it's out.
    !!
    !! Note that there only can be one contact BC for a given link pair. All
    !! contact BCs are identical, so there's no reason to support multiple
    !! contacts for a pair.
    pure logical function is_1contact_node(li)
      integer, intent(in) :: li
      integer :: i, ni, xfi, bc1

      is_1contact_node = .true.
      do i = 1, 2
        ni = nodebc%link(i,li)
        xfi = nodebc%xbcid(ni)
        bc1 = nodebc%bcid(xfi)

        ! If a non-contact bc is here, skip.
        ! If there is any other BC here at all, skip.
        is_1contact_node = is_1contact_node .and. &
            bc1 > bc%dsz &
            .and. nodebc%xbcid(ni+1)-1 == nodebc%xbcid(ni)
      end do

    end function is_1contact_node

  end subroutine c1_init


  subroutine c1_compute(this, displ, ftot, stress_factor)

    class(sm_bc_c1), intent(inout) :: this
    real(r8), intent(in) :: displ(:,:), ftot(:,:), stress_factor(:)

    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, s, tn, l, v(2)

    do i = 1, size(this%index, dim=2)
      n1 = this%index(1,i)
      n2 = this%index(2,i)
      stress1 = dot_product(this%normal(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal(:,i), ftot(:,n2))
      x1 = dot_product(this%normal(:,i), displ(:,n1))
      x2 = dot_product(this%normal(:,i), displ(:,n2))

      s = x2 - x1
      tn = - stress1 / this%area(i) ! TODO-WARN: is the sign right?
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      v(1) = stress2 + this%penalty*(x2 - x1) * stress_factor(n1)
      v(2) = stress1 + this%penalty*(x1 - x2) * stress_factor(n2)

      this%value(:,1,i) = this%normal(:,i) * l * v(1)
      this%value(:,2,i) = this%normal(:,i) * l * v(2)
    end do

  end subroutine c1_compute


  subroutine c1_compute_deriv(this, displ, ftot, stress_factor, diag)

    class(sm_bc_c1), intent(inout) :: this
    real(r8), intent(in) :: displ(:,:), ftot(:,:), stress_factor(:), diag(:,:)

    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, dldu1(3), dldu2(3), diag1(3) !, diag2(3)
    real(r8) :: s, tn, l, dl(2), v(2)

    do i = 1, size(this%index, dim=2)
      ! n1 = this%index(1,i)
      ! n2 = this%index(2,i)
      ! stress1 = dot_product(this%normal(:,i), ftot(:,n1))
      ! stress2 = dot_product(this%normal(:,i), ftot(:,n2))
      ! x1 = dot_product(this%normal(:,i), displ(:,n1))
      ! x2 = dot_product(this%normal(:,i), displ(:,n2))
      ! diag1 = diag(:,n1) * this%normal(:,i) !* stress_factor(n1)
      ! !diag2 = diag(:,n2) * this%normal(:,i) !* stress_factor(n2)

      ! s = x2 - x1
      ! tn = - stress1 / this%area(i)
      ! l = contact_factor(s, tn, this%distance, this%normal_traction)
      ! dl = derivative_contact_factor(s, tn, this%distance, this%normal_traction)
      ! dldu1 = -dl(1)*this%normal(:,i) - dl(2)*diag1 / this%area(i)
      ! dldu2 =  dl(1)*this%normal(:,i)

      ! v(1) = stress2 + this%penalty*(x2 - x1) * stress_factor(n1)
      ! v(2) = stress1 + this%penalty*(x1 - x2) * stress_factor(n2)

      ! this%dvalue(:,1,i) = this%normal(:,i) * (-l*this%penalty*this%normal(:,i) * stress_factor(n1) + dldu1*v(1))
      ! this%dvalue(:,2,i) = this%normal(:,i) * (-l*this%penalty*this%normal(:,i) * stress_factor(n2) + dldu2*v(2))

      ! ! this%dvalue(:,1,i) = - this%normal(:,i)**2 * this%penalty * stress_factor(n1)
      ! ! this%dvalue(:,2,i) = - this%normal(:,i)**2 * this%penalty * stress_factor(n2)

      this%dvalue(:,:,i) = 0
    end do

  end subroutine c1_compute_deriv


  subroutine c1d1_init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized

    class(sm_bc_c1d1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, li, idispl, icontact
    logical :: matching_link

    this%mesh => mesh
    this%penalty = penalty
    this%distance = distance
    this%normal_traction = traction

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      call check_if_c1d1_link(li, idispl, icontact, matching_link)
      if (matching_link) nnode = nnode + 1
    end do
    allocate(this%index(2,nnode), this%value(3,2,nnode), this%dvalue(3,2,nnode))
    allocate(this%normal(3,nnode), this%normal_gap(3,nnode), this%align(3,nnode))
    allocate(this%area(nnode), this%displf(nnode), this%alpha(nnode))

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      call check_if_c1d1_link(li, idispl, icontact, matching_link)
      if (.not.matching_link) cycle
      nnode = nnode + 1
      this%index(:,nnode) = nodebc%node(nodebc%link(:,li))

      this%displf(nnode)%f => bc%displacement(nodebc%bcid(idispl))%f
      this%normal(:,nnode) = nodebc%normal(:,idispl) / norm2(nodebc%normal(:,idispl))
      this%normal_gap(:,nnode) = nodebc%normal(:,icontact) / norm2(nodebc%normal(:,icontact))
      this%area(nnode) = norm2(nodebc%normal(:,icontact))

      this%align(:,nnode) = normalized(cross_product(this%normal(:,nnode), this%normal_gap(:,nnode)))
      this%align(:,nnode) = cross_product(this%align(:,nnode), this%normal(:,nnode))
      this%alpha(nnode) = dot_product(this%align(:,nnode), this%normal_gap(:,nnode))**2
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C1D1 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  contains

    !! Check if the given link index refers to a one contact + one displacement
    !! link pair. Here both sides of the link will have some
    subroutine check_if_c1d1_link(li, idispl, icontact, matching_link)

      integer, intent(in) :: li
      integer, intent(out) :: idispl, icontact
      logical, intent(out) :: matching_link

      integer :: nbc, bcid, xbcid, ni, bc1, bc2

      matching_link = .false.
      idispl = 0
      icontact = 0

      ! Count the BCs applied to the two nodes on either side of the link.
      ! One must have 2 BCs, the other must have
      ni = nodebc%link(1,li)
      nbc = nodebc%xbcid(ni+1) - nodebc%xbcid(ni)
      if (nbc < 1 .or. nbc > 2) return ! if either node has >2 or <1 BCs, the link is disqualified
      do xbcid = nodebc%xbcid(ni), nodebc%xbcid(ni+1)-1
        bcid = nodebc%bcid(xbcid)
        if (bcid > bc%dsz) then
          icontact = xbcid
        else
          idispl = xbcid
        end if
      end do

      ! If one side has two BCs, and they aren't exactly one contact
      ! and one displacement, then this is not a qualifying link.
      if (nbc == 2 .and. idispl == 0) return

      ! If a contact BC wasn't identified on the first node, this mesh
      ! gap doesn't have any contact BC.
      if (icontact == 0) return

      ! Next check the second node in the link. If the displacement condition
      ! hasn't been found yet, it must be found here. We also check that the
      ! contact condition on this node matches the contact condition on the
      ! other. Same for any displacement condition found.
      ni = nodebc%link(2,li)
      nbc = nodebc%xbcid(ni+1) - nodebc%xbcid(ni)
      if (nbc < 1 .or. nbc > 2) return ! if either node has >2 or <1 BCs, the link is disqualified
      do xbcid = nodebc%xbcid(ni), nodebc%xbcid(ni+1)-1
        bcid = nodebc%bcid(xbcid)
        if (bcid > bc%dsz) then
          bc1 = nodebc%bcid(xbcid)
          bc2 = nodebc%bcid(icontact)
          if (bc1 /= bc2) return
        else
          if (idispl == 0) then
            idispl = xbcid
          else
            bc1 = nodebc%bcid(xbcid)
            bc2 = nodebc%bcid(idispl)
            if (bc1 /= bc2) return
          end if
        end if
      end do
      if (idispl == 0) return ! Didn't find a displacement condition

      ! If the link hasn't been disqualified, mark it as valid
      matching_link = .true.

    end subroutine check_if_c1d1_link

  end subroutine c1d1_init

  subroutine c1d1_compute(this, time, displ, ftot, stress_factor)

    class(sm_bc_c1d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)

    integer :: i, n1, n2, k
    real(r8) :: stress1, stress2, x1, x2, s, tn, l, v(2), a(3), b(3), args(0:3)

    args(0) = time
    do i = 1, size(this%index, dim=2)
      n1 = this%index(1,i)
      n2 = this%index(2,i)

      ! displacement part
      k = n1
      this%value(:,1,i) = -dot_product(ftot(:,k), this%normal(:,i)) * this%normal(:,i)
      this%value(:,2,i) = -dot_product(ftot(:,n2), this%normal(:,i)) * this%normal(:,i)

      args(1:) = this%mesh%x(:,k)
      s = this%penalty * stress_factor(k)
      a = this%normal(:,i) * this%displf(i)%eval(args)
      b = dot_product(displ(:,k), this%normal(:,i)) * this%normal(:,i)
      this%value(:,1,i) = this%value(:,1,i) - s * (b - a)
      s = this%penalty * stress_factor(n2)
      b = dot_product(displ(:,n2), this%normal(:,i)) * this%normal(:,i)
      this%value(:,2,i) = this%value(:,2,i) - s * (b - a)

      ! contact part
      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      s = x2 - x1
      tn = - stress1 / this%area(i) ! TODO-WARN: is the sign right?
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      stress1 = dot_product(this%align(:,i), ftot(:,n1))
      stress2 = dot_product(this%align(:,i), ftot(:,n2))
      x1 = dot_product(this%align(:,i), displ(:,n1))
      x2 = dot_product(this%align(:,i), displ(:,n2))
      v(1) = stress2 + this%penalty*this%alpha(i)*(x2 - x1) * stress_factor(n1)
      v(2) = stress1 + this%penalty*this%alpha(i)*(x1 - x2) * stress_factor(n2)

      this%value(:,1,i) = this%value(:,1,i) + this%align(:,i) * l * v(1)
      this%value(:,2,i) = this%value(:,2,i) + this%align(:,i) * l * v(2)
    end do

  end subroutine c1d1_compute

  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine c1d1_compute_deriv(this, time, displ, ftot, stress_factor, diag, F)

    class(sm_bc_c1d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), diag(:,:), F(:,:,:)

    integer :: i, nl, n, d
    real(r8) :: x(3)

    do i = 1, size(this%index, dim=2)
      do nl = 1, 2
        n = this%index(nl,i)
        if (n > this%mesh%nnode_onP) cycle

        this%dvalue(:,nl,i) = 0
        !this%dvalue(:,nl,i) = - this%align(:,i)**2 * this%penalty * stress_factor(n)

        do d = 1,3
          x(d) = dot_product(this%normal(:,i), F(:,d,n))
        end do
        this%dvalue(:,nl,i) = this%dvalue(:,nl,i) &
            - this%normal(:,i) * x &
            - this%penalty * stress_factor(n) * this%normal(:,i)**2
      end do
    end do

  end subroutine c1d1_compute_deriv

end module sm_bc_node_contact_types
