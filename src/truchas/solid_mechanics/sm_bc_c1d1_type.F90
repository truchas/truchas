!!
!! SM_BC_C1D1_TYPE
!!
!! This module implements a type for applying boundary conditions with one
!! gap-contact surface and one Dirichlet displacement.
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

module sm_bc_c1d1_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication, only: global_sum
  use truchas_logging_services
  use scalar_func_containers, only: scalar_func_ptr
  use unstr_mesh_type
  use sm_bc_list_type
  use sm_bc_node_list_type
  use sm_bc_utilities, only: contact_factor, derivative_contact_factor
  use sm_bc_class
  implicit none
  private

  !! One contact and one displacement
  type, extends(sm_bc), public :: sm_bc_c1d1
    private
    real(r8), allocatable, public :: value(:,:), dvalue(:,:), normal_d(:,:)

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal_gap(:,:), align(:,:), alpha(:)
    type(scalar_func_ptr), allocatable :: displf(:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_deriv
  end type sm_bc_c1d1

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized
    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c1d1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, k, icontact(1), idispl(1)
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
    allocate(this%index(nnode), this%value(3,nnode), this%dvalue(3,nnode), &
        this%normal_d(3,nnode), this%displf(nnode), &
        this%normal_gap(3,nnode), this%align(3,nnode), &
        this%linked_node(nnode), this%area(nnode), this%alpha(nnode))
    allocate(this%displacement(nnode), this%traction(nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      this%displf(nnode)%f => bc%displacement(nodebc%bcid(idispl(1)))%f
      this%normal_d(:,nnode)  = nodebc%normal(:,idispl(1)) / norm2(nodebc%normal(:,idispl(1)))

      this%linked_node(nnode) = nodebc%linked_node(icontact(1))
      this%area(nnode) = norm2(nodebc%normal(:,icontact(1)))
      this%normal_gap(:,nnode) = nodebc%normal(:,icontact(1)) / this%area(nnode)

      this%align(:,nnode) = normalized(cross_product(this%normal_d(:,nnode), this%normal_gap(:,nnode)))
      this%align(:,nnode) = normalized(cross_product(this%align(:,nnode), this%normal_d(:,nnode)))
      this%alpha(nnode) = dot_product(this%align(:,nnode), this%normal_gap(:,nnode))**2
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C1D1 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c1d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1, n2
    real(r8) :: args(0:3), stress_penalty, v, l, tn, stress1, stress2, s, x1, x2, x(3)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(i)
      stress_penalty = this%penalty * stress_factor(n1)
      r(:,n1) = r(:,n1) - dot_product(r(:,n1), this%normal_d(:,i)) * this%normal_d(:,i)

      ! displacement part
      args(1:) = this%mesh%x(:,n1)
      v = dot_product(this%normal_d(:,i), displ(:,n1)) - this%displf(i)%eval(args)
      x = this%normal_d(:,i) * v
      r(:,n1) = r(:,n1) - stress_penalty * x

      ! contact part
      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      s = x2 - x1
      tn = - stress1 / this%area(i)
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      stress1 = dot_product(this%align(:,i), ftot(:,n1))
      stress2 = dot_product(this%align(:,i), ftot(:,n2))
      x1 = dot_product(this%align(:,i), displ(:,n1))
      x2 = dot_product(this%align(:,i), displ(:,n2))
      v = stress2 + stress_penalty * this%alpha(i) * (x2 - x1)

      r(:,n1) = r(:,n1) + this%align(:,i) * l * v

      ! visualization storage
      this%displacement(i) = s
      this%traction(i) = tn
    end do

  end subroutine apply


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine apply_deriv(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c1d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, nl, n, d
    real(r8) :: x(3)

    do i = 1, size(this%index)
      n = this%index(i)

      !diag(:,n) = 0
      !diag(:,n) = - this%tangent(:,i)**2 * this%penalty * stress_factor(n)

      do d = 1,3
        x(d) = dot_product(this%normal_d(:,i), F(:,d,n))
      end do
      diag(:,n) = diag(:,n) - this%normal_d(:,i) * x &
          &                 - this%penalty * stress_factor(n) * this%normal_d(:,i)**2
    end do

  end subroutine apply_deriv

end module sm_bc_c1d1_type
