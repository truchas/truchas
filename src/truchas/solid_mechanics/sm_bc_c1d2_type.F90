!!
!! SM_BC_C1D2_TYPE
!!
!! This module implements a type for applying boundary conditions with one
!! gap-contact and two Dirichlet displacements.
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
  use sm_bc_utilities, only: contact_factor
  use sm_bc_class
  implicit none
  private

  !! One contact and two displacements
  type, extends(sm_bc), public :: sm_bc_c1d2
    private
    real(r8), allocatable, public :: value(:,:), dvalue(:,:), tangent(:,:)

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    integer, allocatable :: linked_node(:)
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal_d(:,:,:), normal_gap(:,:), alpha(:), dot(:)
    type(scalar_func_ptr), allocatable :: displf(:,:)
  contains
    procedure :: init
    procedure :: add_graph_links
    procedure :: apply
    procedure :: compute_deriv_diag
    procedure :: compute_deriv_full
  end type sm_bc_c1d2

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use cell_geometry, only: cross_product, normalized
    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c1d2), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, d, icontact(1), idispl(2)
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
    allocate(this%normal_d(3,2,nnode), this%normal_gap(3,nnode), this%tangent(3,nnode))
    allocate(this%linked_node(nnode), this%area(nnode), this%displf(2,nnode), this%alpha(nnode))
    allocate(this%displacement(nnode), this%traction(nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)
      this%linked_node(nnode) = nodebc%linked_node(icontact(1))

      do d = 1, 2
        this%displf(d,nnode)%f => bc%displacement(nodebc%bcid(idispl(d)))%f
        this%normal_d(:,d,nnode)  = nodebc%normal(:,idispl(d)) / norm2(nodebc%normal(:,idispl(d)))
      end do
      this%normal_gap(:,nnode) = nodebc%normal(:,icontact(1)) / norm2(nodebc%normal(:,icontact(1)))
      this%area(nnode) = norm2(nodebc%normal(:,icontact(1)))

      this%tangent(:,nnode) = normalized(cross_product(this%normal_d(:,1,nnode), this%normal_d(:,2,nnode)))
      this%alpha(nnode) = dot_product(this%tangent(:,nnode), this%normal_gap(:,nnode))**2
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C1D2 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c1d2), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, s, tn, l, v, stress_penalty
    real(r8) :: args(0:3), displbc(2), x(3)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      n2 = this%linked_node(i)
      ! the stress_factor gets divided out in compute_residual
      stress_penalty = this%penalty * stress_factor(n1)
      r(:,n1) = dot_product(r(:,n1), this%tangent(:,i)) * this%tangent(:,i)

      ! displacement part
      args(1:) = this%mesh%x(:,n1)
      displbc(1) = this%displf(1,i)%eval(args) ! associated with this%normal_d(:,1,i)
      displbc(2) = this%displf(2,i)%eval(args) ! associated with this%normal_d(:,2,i)
      x = displ(:,n1) - displacement_vector(this%normal_d(:,:,i), displbc)
      x = x - dot_product(x, this%tangent(:,i)) * this%tangent(:,i)
      r(:,n1) = r(:,n1) - stress_penalty * x

      ! contact part
      stress1 = dot_product(this%normal_gap(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal_gap(:,i), ftot(:,n2))
      x1 = dot_product(this%normal_gap(:,i), displ(:,n1))
      x2 = dot_product(this%normal_gap(:,i), displ(:,n2))
      s = x2 - x1
      tn = - stress1 / this%area(i)
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      stress1 = dot_product(this%tangent(:,i), ftot(:,n1))
      stress2 = dot_product(this%tangent(:,i), ftot(:,n2))
      x1 = dot_product(this%tangent(:,i), displ(:,n1))
      x2 = dot_product(this%tangent(:,i), displ(:,n2))
      v = stress2 + stress_penalty * this%alpha(i) * (x2 - x1)

      r(:,n1) = r(:,n1) + this%tangent(:,i) * l * v

      ! visualization storage
      this%displacement(i) = s
      this%traction(i) = tn
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
  subroutine compute_deriv_diag(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c1d2), intent(inout) :: this
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
          &       - this%penalty * (1 - this%tangent(:,i)**2)
    end do

  end subroutine compute_deriv_diag


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine compute_deriv_full(this, time, displ, ftot, stress_factor, Aforce, A)

    use pcsr_matrix_type

    class(sm_bc_c1d2), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    type(pcsr_matrix), intent(in) :: Aforce
    type(pcsr_matrix), intent(inout) :: A

    integer :: i, n, d, ii, jj, n1, n2, n3
    real(r8) :: s, x1, x2, tn, l, stress1, stress2
    real(r8), pointer :: A1(:) => null(), A2(:) => null(), A3(:) => null()
    integer, pointer :: indices(:) => null()

    do i = 1, size(this%index)
      n = this%index(i)
      if (n > this%mesh%nnode_onP) cycle
      n1 = 3*(n-1) + 1
      n2 = 3*(n-1) + 2
      n3 = 3*(n-1) + 3

      ! It is assumed that the indices for each row here are identical.
      ! This *should* be the case.
      call A%get_row_view(n1, A1, indices)
      call A%get_row_view(n2, A2, indices)
      call A%get_row_view(n3, A3, indices)

      ! project out displacement directions
      this%dot = A1 * this%tangent(1,i) + A2 * this%tangent(2,i) + A3 * this%tangent(3,i)
      A1 = this%tangent(1,i) * this%dot
      A2 = this%tangent(2,i) * this%dot
      A3 = this%tangent(3,i) * this%dot

      ! displacement part
      do ii = 1, 3
        call A%add_to(3*(n-1) + ii, 3*(n-1) + ii, -this%penalty)
        do jj = 1, 3
          call A%add_to(3*(n-1) + ii, 3*(n-1) + jj, &
              this%penalty * this%tangent(ii,i) * this%tangent(jj,i))
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
              -this%penalty * l * this%tangent(ii,i) * this%tangent(jj,i) * this%alpha(i))
          call A%add_to(3*(n1-1) + ii, 3*(n2-1) + jj, &
              this%penalty * l * this%tangent(ii,i) * this%tangent(jj,i) * this%alpha(i))
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
        this%dot = this%tangent(1,i) * A1 + this%tangent(2,i) * A2 + this%tangent(3,i) * A3

        do ix = 1, size(indices)
          do ii = 1, 3
            call A%add_to(3*(n1-1) + ii, indices(ix), l * this%tangent(ii,i) * this%dot(ix))
          end do
        end do
      end block

      ! Currently neglecting dldu term (only active when the surfaces are barely in contact)
    end do

  end subroutine compute_deriv_full


  subroutine add_graph_links(this, gforce, g)

    use pcsr_matrix_type
    use sm_bc_utilities, only: alloc_at_least

    class(sm_bc_c1d2), intent(in) :: this
    type(pcsr_graph), intent(in) :: gforce
    type(pcsr_graph), intent(inout) :: g

    integer :: i, j, n1, n2, s, clique(6)
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

end module sm_bc_c1d2_type
