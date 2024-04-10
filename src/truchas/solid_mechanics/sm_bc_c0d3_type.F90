!!
!! SM_BC_C0D3_TYPE
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

module sm_bc_c0d3_type

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

  !! Zero contact and three displacements
  type, extends(sm_bc), public :: sm_bc_c0d3
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: penalty
    real(r8), allocatable :: normal_d(:,:,:)
    type(scalar_func_ptr), allocatable :: displf(:,:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_deriv
  end type sm_bc_c0d3

contains

  subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)

    use sm_bc_utilities, only: check_if_matching_node

    class(sm_bc_c0d3), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, ni, ncontact, d, ndispl
    logical :: matching_node

    integer :: xd, ibc, xibc
    logical :: degenerate

    this%mesh => mesh
    this%penalty = penalty

    nnode = 0
    do ni = 1, size(nodebc%node)
      matching_node = .false.
      block
        ndispl = 0
        do xibc = nodebc%xbcid(ni), nodebc%xbcid(ni+1)-1
          ibc = nodebc%bcid(xibc)
          if (ibc < bc%xcontact) ndispl = ndispl + 1
        end do
        matching_node = ndispl >= 3
      end block
      if (matching_node) nnode = nnode + 1
    end do
    allocate(this%index(nnode), this%normal_d(3,3,nnode), this%displf(3,nnode))

    nnode = 0
    do ni = 1, size(nodebc%node)
      ! For the 3-displacement BC, we'll also consider nodes with any
      ! number of contact BCs. These nodes are overconstrained, so
      ! contact BCs are neglected.
      matching_node = .false.
      block
        ndispl = 0
        do xibc = nodebc%xbcid(ni), nodebc%xbcid(ni+1)-1
          ibc = nodebc%bcid(xibc)
          if (ibc < bc%xcontact) ndispl = ndispl + 1
        end do
        matching_node = ndispl >= 3
      end block
      if (.not.matching_node) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)

      d = 1
      do xibc = nodebc%xbcid(ni), nodebc%xbcid(ni+1)-1
        ibc = nodebc%bcid(xibc)
        if (ibc >= bc%xcontact) cycle
        this%displf(d,nnode)%f => bc%displacement(ibc)%f
        this%normal_d(:,d,nnode) = nodebc%normal(:,xibc) / norm2(nodebc%normal(:,xibc))
        if (d > 1) then
          ! TODO: only add BCs orthogonal to already added BCs.
          ! Can probably do this more efficiently with LAPACK.
          degenerate = .false.
          do xd = 1, d-1
            if (dot_product(this%normal_d(:,d,nnode), this%normal_d(:,xd,nnode)) > 1-1d-2) then
              degenerate = .true.
              exit
            end if
          end do
          if (degenerate) cycle
        end if
        d = d + 1
        if (d > 3) exit
      end do
      INSIST(d > 3)
    end do

    nnode = count(this%index <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    this%enabled = nnode > 0
    if (this%enabled) then
      write(msg,"('SM-C0D3 nodes: ',i6)") nnode
      call TLS_info(trim(msg))
    end if

  end subroutine init


  !! When given three displacements (d1,d2,d3) associated with linearly
  !! independent normal directions, we compute a displacement vector a
  !! such that:
  !!   dot(n1, a) = d1
  !!   dot(n2, a) = d2
  !!   dot(n3, a) = d3
  subroutine apply(this, time, displ, ftot, stress_factor, r)

    external dgesv ! LAPACK

    class(sm_bc_c0d3), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
    real(r8), intent(inout) :: r(:,:)

    integer :: i, n1, stat, ipiv(3)
    real(r8) :: args(0:3), stress_penalty, y(3)

    args(0) = time
    do i = 1, size(this%index)
      n1 = this%index(i)
      stress_penalty = this%penalty * stress_factor(n1)

      args(1:) = this%mesh%x(:,n1)

      y(1) = this%displf(1,i)%eval(args) ! associated with this%normal(:,1,i)
      y(2) = this%displf(2,i)%eval(args) ! associated with this%normal(:,2,i)
      y(3) = this%displf(3,i)%eval(args) ! associated with this%normal(:,3,i)

      call dgesv(3, 1, this%normal_d(:,:,i), 3, ipiv, y, 3, stat)
      INSIST(stat == 0)

      r(:,n1) = - stress_penalty * (displ(:,n1) - y)
    end do

  end subroutine apply


  !! Only the displacement part is currently implemented in the preconditioner.
  subroutine apply_deriv(this, time, displ, ftot, stress_factor, F, diag)

    class(sm_bc_c0d3), intent(inout) :: this
    real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, n

    do i = 1, size(this%index)
      n = this%index(i)
      diag(:,n) = - this%penalty * stress_factor(n)
    end do

  end subroutine apply_deriv

end module sm_bc_c0d3_type
