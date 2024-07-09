!!
!! SM_BC_C0D1_TYPE
!!
!! This module implements a type for applying boundary conditions with one
!! Dirichlet displacement.
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

module sm_bc_c0d1_type

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

  !! Zero contact and one displacement
  type, extends(sm_bc), public :: sm_bc_c0d1
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: penalty
    real(r8), allocatable :: normal_d(:,:), dot(:)
    type(scalar_func_ptr), allocatable :: displf(:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: compute_deriv_diag
    procedure :: compute_deriv_full
  end type sm_bc_c0d1

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
  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c0d1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, icontact(0), idispl(1)
    logical :: matching_node

    this%mesh => mesh
    this%penalty = penalty

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (matching_node) nnode = nnode + 1
    end do
    allocate(this%index(nnode), this%normal_d(3,nnode), this%displf(nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      call check_if_matching_node(ni, nodebc%node(ni), nodebc%bcid, nodebc%xbcid, &
          mesh%nnode_onP, bc%xcontact, icontact, idispl, matching_node)
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      this%displf(nnode)%f => bc%displacement(nodebc%bcid(idispl(1)))%f
      this%normal_d(:,nnode)  = nodebc%normal(:,idispl(1)) / norm2(nodebc%normal(:,idispl(1)))
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C0D1 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  ! Normal displacement BCs enforce a given displacement in the normal
  ! direction, and apply zerotraction in tangential directions. To do this
  ! we must rotate the residual components components to a surface-aligned
  ! coordinate space, and set the normal component of the residual to the
  ! appropriate value. The node normal vector is computed from an
  ! area-weighted average of the surrounding integration points. After
  ! setting the appropriate coordinate-aligned value, it's important to
  ! rotate back into the usual coordinate space, for consistency with the
  ! preconditioner.
  subroutine apply(this, time, displ, ftot, stress_factor, r)

    class(sm_bc_c0d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1
    real(r8) :: args(0:3), stress_penalty, v, x(3)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      ! the stress_factor gets divided out in compute_residual
      stress_penalty = this%penalty * stress_factor(n1)
      r(:,n1) = r(:,n1) - dot_product(r(:,n1), this%normal_d(:,i)) * this%normal_d(:,i)

      ! displacement part
      args(1:) = this%mesh%x(:,n1)
      v = dot_product(this%normal_d(:,i), displ(:,n1)) - this%displf(i)%eval(args)
      x = this%normal_d(:,i) * v
      r(:,n1) = r(:,n1) - stress_penalty * x
    end do

  end subroutine apply


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine compute_deriv_diag(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c0d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, n, d
    real(r8) :: x(3)

    do i = 1, size(this%index)
      n = this%index(i)
      do d = 1,3
        x(d) = dot_product(this%normal_d(:,i), F(:,d,n))
      end do
      diag(:,n) = diag(:,n) - this%normal_d(:,i) * x &
          &                 - this%penalty * this%normal_d(:,i)**2
    end do

  end subroutine compute_deriv_diag


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine compute_deriv_full(this, time, displ, ftot, stress_factor, Aforce, A)

    use pcsr_matrix_type

    class(sm_bc_c0d1), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    type(pcsr_matrix), intent(in) :: Aforce
    type(pcsr_matrix), intent(inout) :: A

    integer :: i, n, d, ii, jj, n1, n2, n3
    real(r8) :: stress_penalty
    real(r8), pointer :: A1(:) => null(), A2(:) => null(), A3(:) => null()
    integer, pointer :: indices1(:) => null(), indices2(:) => null(), indices3(:) => null()

    do i = 1, size(this%index)
      n = this%index(i)
      n1 = 3*(n-1) + 1
      n2 = 3*(n-1) + 2
      n3 = 3*(n-1) + 3
      stress_penalty = this%penalty !/ stress_factor(n)

      ! It is assumed that the indices for each row here are identical. This *should* be the case.
      call A%get_row_view(n1, A1, indices1)
      call A%get_row_view(n2, A2, indices2)
      call A%get_row_view(n3, A3, indices3)
      ASSERT(all(indices1 == indices2) .and. all(indices2 == indices3))

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
    end do

  end subroutine compute_deriv_full

end module sm_bc_c0d1_type
