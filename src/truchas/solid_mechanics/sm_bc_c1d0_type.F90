!!
!! SM_BC_C1D0_TYPE
!!
!! This module implements a type for applying boundary conditions with one
!! gap-contact surface and no Dirichlet displacements.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_c1d0_type

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

  !! One contact and zero displacements
  type, extends(sm_bc), public :: sm_bc_c1d0
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal_gap(:,:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: compute_deriv_diag
    procedure :: compute_deriv_full
  end type sm_bc_c1d0

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c1d0), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, icontact(1), idispl(0)
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
    allocate(this%index(nnode), this%normal_gap(3,nnode), &
        this%linked_node(nnode), this%area(nnode))
    allocate(this%displacement(nnode), this%traction(nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      this%linked_node(nnode) = nodebc%linked_node(icontact(1))
      this%area(nnode) = norm2(nodebc%normal(:,icontact(1)))
      this%normal_gap(:,nnode) = nodebc%normal(:,icontact(1)) / this%area(nnode)
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C1D0 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c1d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1, n2
    real(r8) :: stress_penalty, v, l, tn, stress1, stress2, s, x1, x2

    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(i)
      stress_penalty = this%penalty * stress_factor(n1)

      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      s = x2 - x1
      tn = - stress1 / this%area(i)
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      v = stress2 + stress_penalty * (x2 - x1)

      r(:,n1) = r(:,n1) + this%normal_gap(:,i) * l * v

      ! visualization storage
      this%displacement(i) = s
      this%traction(i) = tn
    end do

  end subroutine apply


  !! This includes the contact contribution to the preconditioner. It is
  !! the only contact type which does so, but it is also likely the most
  !! common gap situation in a simulation.
  subroutine compute_deriv_diag(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c1d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

!NNC, 7/27/2022. Per ZJIBBEN this experiment to improve the PC actually made it
!worse, so its been commented out so as to be a no-op.
!
!    integer :: i, n1, n2
!    real(r8) :: stress_penalty, v, l, tn, stress1, stress2, s, x1, x2, dldu(3), dl(2)
!
!    do i = 1, size(this%index)
!      n1 = this%index(i)
!      n2 = this%linked_node(i)
!      stress_penalty = this%penalty * stress_factor(n1)
!
!      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
!      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
!      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
!      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
!      s = x2 - x1
!      tn = - stress1 / this%area(i)
!      l = contact_factor(s, tn, this%distance, this%normal_traction)
!
!      v = stress2 + stress_penalty * (x2 - x1)
!
!      dl = derivative_contact_factor(s, tn, this%distance, this%normal_traction)
!      dldu = -dl(1)*this%normal_gap(:,i) - dl(2)*this%normal_gap(:,i)*diag(:,n1) / this%area(i)
!
!      ! diag(:,n1) = diag(:,n1) - this%normal_gap(:,i)**2 * l * stress_penalty
!      ! diag(:,n1) = diag(:,n1) + this%normal_gap(:,i) * v * dldu
!    end do

  end subroutine compute_deriv_diag


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine compute_deriv_full(this, time, stress_factor, A)
    use pcsr_matrix_type
    class(sm_bc_c1d0), intent(inout) :: this
    real(r8), intent(in) :: time, stress_factor(:)
    type(pcsr_matrix), intent(inout) :: A
    ! no-op
  end subroutine compute_deriv_full

end module sm_bc_c1d0_type
