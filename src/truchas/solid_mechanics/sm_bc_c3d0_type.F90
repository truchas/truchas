!!
!! SM_BC_C3D0_TYPE
!!
!! This module implements a type for applying boundary conditions with three
!! gap-contact and zero Dirichlet displacements.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! June 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_c3d0_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication, only: global_sum
  use truchas_logging_services
  use unstr_mesh_type
  use sm_bc_list_type
  use sm_bc_node_list_type
  use sm_bc_utilities, only: contact_factor
  use sm_bc_class
  implicit none
  private

  !! Three contacts and zero displacements
  type, extends(sm_bc), public :: sm_bc_c3d0
    private
    real(r8), allocatable, public :: value(:,:), dvalue(:,:)

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:,:), normal(:,:,:), tangent(:,:,:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: compute_deriv_diag
    procedure :: compute_deriv_full
  end type sm_bc_c3d0

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized
    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c3d0), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, k, icontact(3), idispl(0), linked_node(3)
    logical :: matching_node

    this%mesh => mesh
    this%penalty = penalty
    this%distance = distance
    this%normal_traction = traction

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (matching_node) nnode = nnode + 1
    end do
    allocate(this%index(nnode), this%value(3,nnode), this%dvalue(3,nnode))
    allocate(this%normal(3,3,nnode), this%tangent(3,3,nnode))
    allocate(this%linked_node(nnode), this%area(3,nnode))
    allocate(this%displacement(nnode), this%traction(nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      this%linked_node(nnode) = nodebc%linked_node(icontact(1))
      do k = 1, 3
        linked_node(k) = nodebc%linked_node(icontact(k))
        this%area(k,nnode) = norm2(nodebc%normal(:,icontact(k)))
        this%normal(:,k,nnode) = nodebc%normal(:,icontact(k)) / this%area(k,nnode)
      end do

      ! This will zero out the tangent if the normals are identical.
      ! We might want to change this to zero out the tangent if the
      ! normals are within some angle of each other.
      this%tangent(:,1,nnode) = normalized(cross_product(this%normal(:,1,nnode), this%normal(:,2,nnode)))
      this%tangent(:,2,nnode) = normalized(cross_product(this%normal(:,2,nnode), this%normal(:,3,nnode)))
      this%tangent(:,3,nnode) = normalized(cross_product(this%normal(:,3,nnode), this%normal(:,1,nnode)))

      ! Restricted to the case where there is one node across all gap surfaces.
      INSIST(linked_node(1) == linked_node(2) .and. linked_node(1) == linked_node(3))
      ASSERT(this%linked_node(nnode) /= 0)
      ASSERT(all(this%area(:,nnode) /= 0))
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C3D0 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c3d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1, n2, li
    real(r8) :: stress1, stress2, delta, lambda(3), v, stress_penalty
    real(r8) :: x(3), xn(3,3), xt(3,3)

    do i = 1, size(this%index)
      ! Restricted to the case where there is one node across all gap surfaces.
      n1 = this%index(i)
      n2 = this%linked_node(i)
      stress_penalty = this%penalty * stress_factor(n1)

      ! See the Physics & Algorithms manual, SM BC appendix.
      do li = 1, 3
        call compute_contact_variables(li, i, lambda(li), delta, stress1, stress2)
        v = stress2 + stress_penalty * delta
        r(:,n1) = r(:,n1) + lambda(li) * this%normal(:,li,i) * v
      end do

      x = ftot(:,n2) + stress_penalty * (displ(:,n2) - displ(:,n1))
      do li = 1, 3
        xn(:,li) = this%normal(:,li,i) * dot_product(this%normal(:,li,i), x)
        xt(:,li) = this%tangent(:,li,i) * dot_product(this%tangent(:,li,i), x)
      end do
      r(:,n1) = r(:,n1) &
          + lambda(1)*lambda(2)*(x - xt(:,1) - xn(:,1) - xn(:,2)) &
          + lambda(2)*lambda(3)*(x - xt(:,2) - xn(:,2) - xn(:,3)) &
          + lambda(3)*lambda(1)*(x - xt(:,3) - xn(:,3) - xn(:,1)) &
          + lambda(1)*lambda(2)*lambda(3)*(2*x - sum(xt,dim=2) - sum(xn,dim=2))
    end do

  contains

    subroutine compute_contact_variables(li, i, lambda, delta, stress1, stress2)
      integer, intent(in) :: li, i
      real(r8), intent(out) :: lambda, delta, stress1, stress2
      integer :: n1, n2
      real(r8) :: x1, x2, tn
      n1 = this%index(i)
      n2 = this%linked_node(i)
      stress1 = dot_product(this%normal(:,li,i), ftot(:,n1))
      stress2 = dot_product(this%normal(:,li,i), ftot(:,n2))
      x1 = dot_product(this%normal(:,li,i), displ(:,n1))
      x2 = dot_product(this%normal(:,li,i), displ(:,n2))
      delta = x2 - x1
      tn = - stress1 / this%area(li,i)
      lambda = contact_factor(delta, tn, this%distance, this%normal_traction)

      ! visualization storage
      ! Last one written is what gets stored to output. See comments in
      ! sm_bc_type::compute_viz_fields.
      this%displacement(i) = delta
      this%traction(i) = tn
    end subroutine compute_contact_variables

  end subroutine apply


  !! Contact preconditioner contribution currently not implemented.
  subroutine compute_deriv_diag(this, time, displ, ftot, stress_factor, F, diag)
    class(sm_bc_c3d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)
    ! do nothing
  end subroutine compute_deriv_diag


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine compute_deriv_full(this, time, displ, ftot, stress_factor, Aforce, A)
    use pcsr_matrix_type
    class(sm_bc_c3d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    type(pcsr_matrix), intent(in) :: Aforce
    type(pcsr_matrix), intent(inout) :: A
    ! no-op
  end subroutine compute_deriv_full

end module sm_bc_c3d0_type
