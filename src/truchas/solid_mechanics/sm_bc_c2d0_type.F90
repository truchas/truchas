!!
!! SM_BC_C2D0_TYPE
!!
!! This module implements a type for applying boundary conditions with two
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

module sm_bc_c2d0_type

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

  !! Two contacts and zero displacements
  type, extends(sm_bc), public :: sm_bc_c2d0
    private
    real(r8), allocatable, public :: value(:,:), dvalue(:,:), tangent(:,:)

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:,:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:,:), normal(:,:,:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_deriv
  end type sm_bc_c2d0

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized
    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c2d0), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, d, icontact(2), idispl(0)
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
    allocate(this%normal(3,2,nnode), this%tangent(3,nnode))
    allocate(this%linked_node(2,nnode), this%area(2,nnode))
    allocate(this%displacement(nnode), this%traction(nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      do d = 1, 2
        this%linked_node(d,nnode) = nodebc%linked_node(icontact(d))
        this%area(d,nnode) = norm2(nodebc%normal(:,icontact(d)))
        this%normal(:,d,nnode) = nodebc%normal(:,icontact(d)) / this%area(d,nnode)
      end do

      ! This will zero out the tangent if the normals are identical.
      ! We might want to change this to zero out the tangent if the
      ! normals are within some angle of each other.
      this%tangent(:,nnode) = normalized(cross_product(this%normal(:,1,nnode), this%normal(:,2,nnode)))

      ASSERT(all(this%linked_node(:,nnode) /= 0))
      ASSERT(all(this%area(:,nnode) /= 0))
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C2D0 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c2d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1, n2, n3, li
    real(r8) :: stress2(2), delta(2), lambda(2), v, stress_penalty

    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(1,i)
      n3 = this%linked_node(2,i)
      stress_penalty = this%penalty * stress_factor(n1)

      call compute_contact_variables(1, i, lambda(1), delta(1), stress2(1))
      call compute_contact_variables(2, i, lambda(2), delta(2), stress2(2))

      ! Here are all the different cases that can show up... See the
      ! reference manual for details.
      if (all(this%tangent(:,i) == 0)) then
        ! Two surfaces, two nodes, but only one normal.
        do li = 1, 2
          v = stress2(li) + stress_penalty * delta(li)
          r(:,n1) = r(:,n1) + lambda(li) * this%normal(:,li,i) * v
        end do

      else if (any(lambda == 0)) then
        ! Only in contact with one surface, or no surfaces.
        li = maxloc(lambda, dim=1)
        v = stress2(li) + stress_penalty * delta(li)
        r(:,n1) = r(:,n1) + lambda(li) * this%normal(:,li,i) * v

      else if (n2 == n3) then
        ! Corner contact, two surfaces but one node -- see reference manual.
        do li = 1, 2
          v = stress2(li) + stress_penalty * delta(li)
          r(:,n1) = r(:,n1) + lambda(li) * this%normal(:,li,i) * v
        end do

        block
          real(r8) :: x(3), xn1(3), xn2(3), xt(3)
          x = ftot(:,n2) + stress_penalty * (displ(:,n2) - displ(:,n1))
          xn1 = this%normal(:,1,i) * dot_product(this%normal(:,1,i), x)
          xn2 = this%normal(:,2,i) * dot_product(this%normal(:,2,i), x)
          xt = this%tangent(:,i) * dot_product(this%tangent(:,i), x)
          x = x - xn1 - xn2 - xt
          r(:,n1) = r(:,n1) + lambda(1)*lambda(2)*x
        end block

      else
        ! Two surfaces, two nodes, two normals
        do li = 1, 2
          v = stress2(li) + stress_penalty * delta(li)
          r(:,n1) = r(:,n1) + lambda(li) * this%normal(:,li,i) * v
        end do

        block
          real(r8) :: x(3), c21(3), c31(3)
          c21 = ftot(:,n2) + stress_penalty * (displ(:,n2) - displ(:,n1))
          c31 = ftot(:,n3) + stress_penalty * (displ(:,n3) - displ(:,n1))
          x = c21 + c31
          x = x - this%tangent(:,i) * dot_product(this%tangent(:,i), x)
          x = x - this%normal(:,1,i) * dot_product(this%normal(:,1,i), c21)
          x = x - this%normal(:,2,i) * dot_product(this%normal(:,2,i), c31)
          r(:,n1) = r(:,n1) + lambda(1)*lambda(2)*x
        end block

      end if
    end do

  contains

    subroutine compute_contact_variables(li, i, lambda, delta, stress2)
      integer, intent(in) :: li, i
      real(r8), intent(out) :: lambda, delta, stress2
      integer :: n1, n2
      real(r8) :: x1, x2, tn, stress1
      n1 = this%index(i)
      n2 = this%linked_node(li,i)
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
  subroutine apply_deriv(this, time, displ, ftot, stress_factor, F, diag)
    class(sm_bc_c2d0), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)
    ! do nothing
  end subroutine apply_deriv

end module sm_bc_c2d0_type
