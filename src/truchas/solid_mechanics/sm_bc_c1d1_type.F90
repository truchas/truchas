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
  use sm_bc_utilities, only: contact_factor
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
    real(r8), allocatable :: area(:), normal_gap(:,:), align(:,:), alpha(:), dot(:)
    type(scalar_func_ptr), allocatable :: displf(:)
  contains
    procedure :: init
    procedure :: add_graph_links
    procedure :: apply
    procedure :: compute_deriv_diag
    procedure :: compute_deriv_full
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
    integer :: nnode, ni, icontact(1), idispl(1)
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

    nnode = size(this%index)
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
      ! the stress_factor gets divided out in compute_residual
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


  !! Only the displacement part is currently implemented in the Jacobi preconditioner.
  subroutine compute_deriv_diag(this, time, displ, ftot, stress_factor, F, diag)

    use sm_bc_utilities, only: derivative_contact_factor

    class(sm_bc_c1d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, d, n1, n2
    real(r8) :: x(3), v, l, tn, stress1, stress2, s, x1, x2, dldu(3), dl(2), stress_penalty

    do i = 1, size(this%index)
      n1 = this%index(i)

      do d = 1,3
        x(d) = dot_product(this%normal_d(:,i), F(:,d,n1))
      end do
      diag(:,n1) = diag(:,n1) - this%normal_d(:,i) * x &
          &                   - this%penalty * this%normal_d(:,i)**2

      ! ! contact part
      ! n2 = this%linked_node(i)
      ! stress_penalty = this%penalty * stress_factor(n1)

      ! stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      ! stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      ! x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      ! x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      ! s = x2 - x1
      ! tn = - stress1 / this%area(i)
      ! l = contact_factor(s, tn, this%distance, this%normal_traction)

      ! stress1 = dot_product(this%align(:,i), ftot(:,n1))
      ! stress2 = dot_product(this%align(:,i), ftot(:,n2))
      ! x1 = dot_product(this%align(:,i), displ(:,n1))
      ! x2 = dot_product(this%align(:,i), displ(:,n2))
      ! v = stress2 + stress_penalty * this%alpha(i) * s

      ! dl = derivative_contact_factor(s, tn, this%distance, this%normal_traction)
      ! dldu = -dl(1)*this%normal_gap(:,i) - dl(2)*this%normal_gap(:,i)*diag(:,n1) / this%area(i)

      ! diag(:,n1) = diag(:,n1) - this%align(:,i)**2 * l * this%penalty * this%alpha(i)
      ! diag(:,n1) = diag(:,n1) + this%align(:,i) * v * dldu / stress_factor(n1)
    end do

  end subroutine compute_deriv_diag


  subroutine compute_deriv_full(this, time, displ, ftot, stress_factor, Aforce, A)

    use pcsr_matrix_type

    class(sm_bc_c1d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    type(pcsr_matrix), intent(in) :: Aforce
    type(pcsr_matrix), intent(inout) :: A

    integer :: i, n, d, ii, jj, n1, n2, n3
    real(r8) :: s, x1, x2, tn, l, stress1, stress2
    real(r8), pointer :: A1(:) => null(), A2(:) => null(), A3(:) => null()
    integer, pointer :: indices(:) => null()

    do i = 1, size(this%index)
      n = this%index(i)
      n1 = 3*(n-1) + 1
      n2 = 3*(n-1) + 2
      n3 = 3*(n-1) + 3

      ! It is assumed that the indices for each row here are identical. This *should* be the case.
      call A%get_row_view(n1, A1, indices)
      call A%get_row_view(n2, A2, indices)
      call A%get_row_view(n3, A3, indices)

      ! project out displacement direction
      this%dot = A1 * this%normal_d(1,i) + A2 * this%normal_d(2,i) + A3 * this%normal_d(3,i)
      A1 = A1 - this%normal_d(1,i) * this%dot
      A2 = A2 - this%normal_d(2,i) * this%dot
      A3 = A3 - this%normal_d(3,i) * this%dot

      ! displacement part
      do ii = 1, 3
        do jj = 1, 3
          call A%add_to(3*(n-1) + ii, 3*(n-1) + jj, &
              -this%penalty * this%normal_d(ii,i) * this%normal_d(jj,i))
        end do
      end do

      ! contact part
      n1 = this%index(i)
      n2 = this%linked_node(i)

      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      s = x2 - x1
      tn = - stress1 / this%area(i)
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      do ii = 1, 3
        do jj = 1, 3
          call A%add_to(3*(n1-1) + ii, 3*(n1-1) + jj, &
              -this%penalty * l * this%align(ii,i) * this%align(jj,i) * this%alpha(i))
          call A%add_to(3*(n1-1) + ii, 3*(n2-1) + jj, &
              this%penalty * l * this%align(ii,i) * this%align(jj,i) * this%alpha(i))
        end do
      end do

      ! linked stress part
      block
        real(r8), pointer :: A1(:) => null(), A2(:) => null(), A3(:) => null()
        integer, pointer :: indices(:) => null()
        integer :: ix

        ! These indices at the linked node are expected to be identical to each other,
        ! but not identical to the indices for the local node.
        ! NB: The stress factor is already taken into account in Aforce.
        call Aforce%get_row_view(3*(n2-1)+1, A1, indices)
        call Aforce%get_row_view(3*(n2-1)+2, A2, indices)
        call Aforce%get_row_view(3*(n2-1)+3, A3, indices)
        this%dot = this%align(1,i) * A1 + this%align(2,i) * A2 + this%align(3,i) * A3

        do ix = 1, size(indices)
          do ii = 1, 3
            call A%add_to(3*(n1-1) + ii, indices(ix), l * this%align(ii,i) * this%dot(ix))
          end do
        end do
      end block

      ! ! dldu term. Only active when the surfaces are barely in contact, so should be negligible.
      ! block
      !   use sm_bc_utilities, only: derivative_contact_factor
      !   real(r8) :: dl(2), v
      !   real(r8), pointer :: A1(:) => null(), A2(:) => null(), A3(:) => null()
      !   integer, pointer :: indices(:) => null()

      !   dl = derivative_contact_factor(s, tn, this%distance, this%normal_traction)
      !   stress1 = dot_product(this%align(:,i), ftot(:,n1))
      !   stress2 = dot_product(this%align(:,i), ftot(:,n2))
      !   x1 = dot_product(this%align(:,i), displ(:,n1))
      !   x2 = dot_product(this%align(:,i), displ(:,n2))
      !   s = x2 - x1
      !   v = stress2 / stress_factor(n1) + this%penalty * s

      !   ! It is assumed that the indices for each row here are identical. This *should* be the case.
      !   call A%get_row_view(3*(n1-1)+1, A1, indices)
      !   call A%get_row_view(3*(n1-1)+2, A2, indices)
      !   call A%get_row_view(3*(n1-1)+3, A3, indices)

      !   this%dot = A1 * this%normal_gap(1,i) + A2 * this%normal_gap(2,i) + A3 * this%normal_gap(3,i)
      !   A1 = A1 - this%align(1,i) * v * (dl(1) + dl(2) / this%area(i) * this%dot)
      !   A2 = A2 - this%align(2,i) * v * (dl(1) + dl(2) / this%area(i) * this%dot)
      !   A3 = A3 - this%align(3,i) * v * (dl(1) + dl(2) / this%area(i) * this%dot)
      ! end block
    end do

  end subroutine compute_deriv_full


  subroutine add_graph_links(this, gforce, g)

    use pcsr_matrix_type
    use sm_bc_utilities, only: alloc_at_least

    class(sm_bc_c1d1), intent(in) :: this
    type(pcsr_graph), intent(in) :: gforce
    type(pcsr_graph), intent(inout) :: g

    integer :: i, j, n1, n2, clique(6), s
    integer, allocatable :: cliques(:)

    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(i)

      ! contact-displacement part
      do j = 1, 3
        clique(j) = 3*(n1-1) + j
        clique(3+j) = 3*(n2-1) + j
      end do
      call g%add_clique(clique)

      ! contact-stress part
      s = 3 + gforce%xadj(3*(n2-1) + 2) - gforce%xadj(3*(n2-1) + 1)
      call alloc_at_least(cliques, s)
      cliques(1) = 3*(n1 - 1) + 1
      cliques(2) = 3*(n1 - 1) + 2
      cliques(3) = 3*(n1 - 1) + 3
      cliques(4:s) = gforce%adjncy(gforce%xadj(3*(n2-1) + 1):gforce%xadj(3*(n2-1) + 2)-1)
      call g%add_clique(cliques(:s))
    end do

  end subroutine add_graph_links

end module sm_bc_c1d1_type
