!!
!! SM_BC_C0D2_TYPE
!!
!! This module implements a type for applying boundary conditions with two
!! Dirichlet displacements.
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

module sm_bc_c0d2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication, only: global_sum
  use truchas_logging_services
  use scalar_func_containers, only: scalar_func_ptr
  use unstr_mesh_type
  use sm_bc_list_type
  use sm_bc_node_list_type
  use sm_bc_class
  implicit none
  private

  !! Zero contact and two displacements
  type, extends(sm_bc), public :: sm_bc_c0d2
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: penalty
    real(r8), allocatable :: normal_d(:,:,:), tangent(:,:)
    type(scalar_func_ptr), allocatable :: displf(:,:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_deriv
  end type sm_bc_c0d2

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized
    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c0d2), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, icontact(0), idispl(2), d
    logical :: matching_node

    this%mesh => mesh
    this%penalty = penalty

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (matching_node) nnode = nnode + 1
    end do
    allocate(this%index(nnode), this%normal_d(3,2,nnode), this%tangent(3,nnode), &
        this%displf(2,nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      do d = 1, 2
        this%displf(d,nnode)%f => bc%displacement(nodebc%bcid(idispl(d)))%f
        this%normal_d(:,d,nnode) = nodebc%normal(:,idispl(d)) / norm2(nodebc%normal(:,idispl(d)))
      end do
      this%tangent(:,nnode) = normalized(cross_product(this%normal_d(:,1,nnode), this%normal_d(:,2,nnode)))
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C0D2 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  ! 2-normal displacement BC
  ! Here the residual is projected into the space tangent to the two given
  ! displacement directions, and a new term in the "displacement space" is
  ! added to the residual, which ensures the node's displacement matches
  ! the BCs.
  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c0d2), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1
    real(r8) :: args(0:3), stress_penalty, y(2), x(3)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      stress_penalty = this%penalty * stress_factor(n1)
      r(:,n1) = dot_product(r(:,n1), this%tangent(:,i)) * this%tangent(:,i)

      args(1:) = this%mesh%x(:,n1)

      y(1) = this%displf(1,i)%eval(args) ! associated with this%normal(:,1,i)
      y(2) = this%displf(2,i)%eval(args) ! associated with this%normal(:,2,i)

      x = displ(:,n1) - displacement_vector(this%normal_d(:,:,i), y)
      x = x - dot_product(x, this%tangent(:,i)) * this%tangent(:,i)
      r(:,n1) = r(:,n1) - stress_penalty * x
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

  end subroutine apply


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine apply_deriv(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c0d2), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, n, d
    real(r8) :: x(3)

    do i = 1, size(this%index)
      n = this%index(i)
      do d = 1,3
        x(d) = dot_product(this%tangent(:,i), F(:,d,n))
      end do
      diag(:,n) = this%tangent(:,i) * x &
          &       - this%penalty * stress_factor(n) * (1 - this%tangent(:,i)**2)
    end do

  end subroutine apply_deriv

end module sm_bc_c0d2_type
