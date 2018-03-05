!!
!! FISCHER_TYPE
!!
!! Implements the Fischer intial guess technique from
!!
!! P. F. Fischer, Projection Techniques for iterative solution of Ax=b with
!!  successive right-hand sides.  Comput. Methods Appl. Mech. Engrg. v. 163
!!  pp. 193--204, 1998
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type FISCHER_GUESS that encapsulates the data
!!  and procedures for generating an solution to the pressure poisson system
!!
!!  Interaction with this object occurs through its methods, generally via the flow_driver
!!
!!  INIT(FLOW_MESH MESH) - Allocates the data members for this object and initializes the
!!    internal solver
!!
!!  GUESS(REAL(R8) GUESSED_SOLUTION(:), REAL(R8) RHS(:))
!!    write GUESSED_SOLUTION based on RHS
!!
!!  UPDATE(REAL(R8) SOLUTION(:), REAL(R8) RHS(:), PCSR_MATRIX LHS)
!!    update internals with the computed solution
!!
!! USAGE (after init)
!!
!!    call fischer_guess%guess( x', b)
!!    solve Ax=b using x'
!!    call fischer_guess%update( x, b, A)
!!

#include "f90_assert.fpp"

module fischer_guess_type

  use kinds, only: r8
  use flow_mesh_type
  use index_partitioning
  use parallel_communication
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  type, public :: fischer_guess
    private
    type(flow_mesh), pointer :: mesh ! unowned reference
    real(r8), allocatable :: b_tilde(:,:), x_tilde(:,:), x_guess(:)
    integer :: size, max_size
    integer :: lri ! least-recently-inserted
  contains
    procedure :: read_params
    procedure :: init
    procedure :: guess
    procedure :: update
  end type fischer_guess

contains

  subroutine read_params(this, p)
    class(fischer_guess), intent(inout) :: this
    type(parameter_list), intent(inout) :: p

    call p%get('fischer_history', this%max_size, 6)

  end subroutine read_params


  subroutine init (this, mesh)
    class(fischer_guess), intent(inout) :: this
    type(flow_mesh), pointer, intent(in)  :: mesh
    !-
    integer :: n

    if (this%max_size < 1) return

    this%mesh => mesh
    this%size = 0
    this%lri = 1
    n = mesh%mesh%ncell

    allocate(this%b_tilde(n,this%max_size+1), this%x_tilde(n,this%max_size+1), &
        this%x_guess(n))

    this%x_guess = 0.0_r8
    this%b_tilde = 0.0_r8
    this%x_tilde = 0.0_r8

  end subroutine init


  ! Given the new right hand side b, this computes a new guess for
  ! the solution x based on previous stored solution values
  subroutine guess (this, b, x)
    class(fischer_guess), intent(inout) :: this
    real(r8),       intent(out)   :: x(:)
    real(r8),       intent(in)    :: b(:)

    real(r8) :: alpha(this%size)
    integer :: i, j, nop

    x = 0.0_r8
    if (this%max_size < 1) return

    nop = this%mesh%mesh%ncell_onP
    this%x_guess = 0.0_r8

    do i = 1, this%size
      alpha(i) = global_dot_product(b(1:nop), this%b_tilde(1:nop,i))
    end do

    do i = 1, this%size
      do j = 1, nop
        this%x_guess(j) = this%x_guess(j)+alpha(i)*this%x_tilde(j,i)
      end do
    end do
    x = this%x_guess

  end subroutine guess

  ! After a PPE solution, the set of projection vectors
  ! (i.e. previous solutions and right hand sides) is updated.
  subroutine update (this, b, x_soln, lhs)
    class(fischer_guess), intent(inout) :: this
    real(r8), intent(in) :: x_soln(:), b(:)
    type(pcsr_matrix), pointer, intent(in) :: lhs

    integer :: i, j, idx, nop
    real(r8) :: b_norm, alpha(this%max_size)
    real(r8), pointer :: mval(:)
    integer, pointer :: midx(:)

    if (this%max_size < 1) return

    nop = this%mesh%mesh%ncell_onP
    idx = this%size+1

    this%x_tilde(:,idx) = x_soln - this%x_guess
    call gather_boundary(this%mesh%mesh%cell_ip, this%x_tilde(:,idx))

    do i = 1, nop
      call lhs%get_row_view(i, mval, midx)
      this%b_tilde(i, idx) = dot_product(mval,this%x_tilde(midx,idx))
    end do

    do i = 1, this%size
      alpha(i) = global_dot_product(this%b_tilde(1:nop,idx), this%b_tilde(1:nop,i))
      do j = 1, nop
        this%b_tilde(j,idx) = this%b_tilde(j,idx) - alpha(i)*this%b_tilde(j,i)
        this%x_tilde(j,idx) = this%x_tilde(j,idx) - alpha(i)*this%x_tilde(j,i)
      end do
    end do

    b_norm = sqrt(global_dot_product(this%b_tilde(1:nop,idx), this%b_tilde(1:nop,idx)))

    if (b_norm > epsilon(1.0_r8)) then
      if (this%size == this%max_size) then
        ! reuse the least-recently-inserted entry
        this%b_tilde(:,this%lri) = this%b_tilde(:,idx)/b_norm
        this%x_tilde(:,this%lri) = this%x_tilde(:,idx)/b_norm
        this%lri = 1+mod(this%lri, this%max_size)
      else
        this%b_tilde(:,idx) = this%b_tilde(:,idx)/b_norm
        this%x_tilde(:,idx) = this%x_tilde(:,idx)/b_norm
        this%size = idx
      end if
    end if
  end subroutine update

end module fischer_guess_type
