!!
!! SM_BC_C2D1_TYPE
!!
!! This module implements a type for applying boundary conditions with two
!! gap-contact surfaces and one Dirichlet displacement.
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

module sm_bc_c2d1_type

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

  !! Two contacts and one displacement
  type, extends(sm_bc), public :: sm_bc_c2d1
    private
    real(r8), allocatable, public :: value(:,:), dvalue(:,:), normal_d(:,:)

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:,:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:,:), normal_gap(:,:,:), ortho(:,:,:), alpha(:,:)
    type(scalar_func_ptr), allocatable :: displf(:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_deriv
  end type sm_bc_c2d1

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized
    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c2d1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, k, icontact(2), idispl(1)
    logical :: matching_node
    real(r8) :: tangent(3)

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
        this%normal_gap(3,2,nnode), this%ortho(3,2,nnode), &
        this%linked_node(2,nnode), this%area(2,nnode), this%alpha(2,nnode))
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
      do k = 1, 2
        this%linked_node(k,nnode) = nodebc%linked_node(icontact(k))
        this%area(k,nnode) = norm2(nodebc%normal(:,icontact(k)))
        this%normal_gap(:,k,nnode) = nodebc%normal(:,icontact(k)) / this%area(k,nnode)

        tangent = normalized(cross_product(this%normal_d(:,nnode), this%normal_gap(:,k,nnode)))
        this%ortho(:,k,nnode) = normalized(cross_product(tangent, this%normal_d(:,nnode)))
        this%alpha(k,nnode) = dot_product(this%ortho(:,k,nnode), this%normal_gap(:,k,nnode))**2
      end do
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C2D1 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c2d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, li, n1, n2, n3
    real(r8) :: stress_penalty, v
    real(r8), dimension(2) :: lambda, delta, stress1, stress2
    real(r8) :: args(0:3), x(3), xn(3), c(3,2)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(1,i)
      n3 = this%linked_node(2,i)
      stress_penalty = this%penalty * stress_factor(n1)
      r(:,n1) = r(:,n1) - this%normal_d(:,i) * dot_product(this%normal_d(:,i), r(:,n1))

      ! displacement part
      args(1:) = this%mesh%x(:,n1)
      v = dot_product(this%normal_d(:,i), displ(:,n1)) - this%displf(i)%eval(args)
      x = this%normal_d(:,i) * v
      r(:,n1) = r(:,n1) - stress_penalty * x

      ! contact part
      call compute_contact_variables(1, i, lambda(1), delta(1), stress1(1), stress2(1))
      call compute_contact_variables(2, i, lambda(2), delta(2), stress1(2), stress2(2))

      if (any(lambda == 0)) then
        ! contact with one surface, or no surfaces
        li = maxloc(lambda, dim=1)
        x = ftot(:,n2) + stress_penalty * this%alpha(li,i) * (displ(:,n2) - displ(:,n1))
        x = this%ortho(:,1,i) * dot_product(this%ortho(:,1,i), x)
        r(:,n1) = r(:,n1) + lambda(li) * x

      else if (n2 == n3) then
        ! contact with two surfaces, one node
        do li = 1, 2
          x = ftot(:,n2) + this%alpha(li,i) * stress_penalty * (displ(:,n2) - displ(:,n1))
          xn = this%ortho(:,li,i) * dot_product(this%ortho(:,li,i), x)
          r(:,n1) = r(:,n1) + lambda(li) * xn
        end do

        ! cross terms
        x = ftot(:,n2)
        r(:,n1) = r(:,n1) + lambda(1)*lambda(2) * &
            (x &
            - this%normal_d(:,i)*dot_product(this%normal_d(:,i),x) &
            - this%ortho(:,1,i)*dot_product(this%ortho(:,1,i),x)  &
            - this%ortho(:,2,i)*dot_product(this%ortho(:,2,i),x))

        x = stress_penalty * (displ(:,n2) - displ(:,n1))
        r(:,n1) = r(:,n1) + lambda(1)*lambda(2) * &
            (sum(this%alpha(:,i))*(x - this%normal_d(:,i)*dot_product(this%normal_d(:,i),x)) &
            + this%alpha(1,i)*this%ortho(:,1,i)*dot_product(this%ortho(:,1,i),x) &
            + this%alpha(2,i)*this%ortho(:,2,i)*dot_product(this%ortho(:,2,i),x))

      else
        ! contact with two surfaces, two nodes
        c(:,1) = ftot(:,n2) + this%alpha(1,i) * stress_penalty * (displ(:,n2) - displ(:,n1))
        c(:,2) = ftot(:,n3) + this%alpha(2,i) * stress_penalty * (displ(:,n3) - displ(:,n1))
        do li = 1, 2
          xn = this%ortho(:,li,i) * dot_product(this%ortho(:,li,i), c(:,li))
          r(:,n1) = r(:,n1) + lambda(li) * xn
        end do

        if (dot_product(this%normal_gap(:,1,i), this%normal_gap(:,2,i)) /= 1) then
          ! two normals -- include cross terms to disallow sliding
          xn = this%normal_d(:,i) * dot_product(this%normal_d(:,i), c(:,1) + c(:,2))
          r(:,n1) = r(:,n1) + lambda(1)*lambda(2) * (c(:,1) + c(:,2) - xn)
          do li = 1, 2
            xn = this%ortho(:,li,i) * dot_product(this%ortho(:,li,i), c(:,li))
            r(:,n1) = r(:,n1) - lambda(1)*lambda(2) * xn
          end do
        end if
      end if
    end do

  contains

    subroutine compute_contact_variables(li, i, lambda, delta, stress1, stress2)
      integer, intent(in) :: li, i
      real(r8), intent(out) :: lambda, delta, stress1, stress2
      integer :: n1, n2
      real(r8) :: x1, x2, tn
      n1 = this%index(i)
      n2 = this%linked_node(li,i)
      stress1 = dot_product(this%normal_gap(:,li,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,li,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,li,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,li,i), displ(:,n2))
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


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine apply_deriv(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c2d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, nl, n, d
    real(r8) :: x(3)

    do i = 1, size(this%index)
      n = this%index(i)

      diag(:,n) = diag(:,n)
      !diag(:,n) = - this%tangent(:,i)**2 * this%penalty * stress_factor(n)

      do d = 1,3
        x(d) = dot_product(this%normal_d(:,i), F(:,d,n))
      end do
      diag(:,n) = diag(:,n) - this%normal_d(:,i) * x &
          &                 - this%penalty * stress_factor(n) * this%normal_d(:,i)**2
    end do

  end subroutine apply_deriv

end module sm_bc_c2d1_type
