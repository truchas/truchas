!!
!! FLOW_2D_MOMENTUM_SOLVER_TYPE
!!
!! This module defines FLOW_2D_MOMENTUM_SOLVER, the linear-solver adapter for
!! FLOW_2D_MOMENTUM. It owns a BLOCK_PCSR_BRIDGE and HYPRE_HYBRID while the
!! momentum operator remains stored and assembled as a 2-by-2 block matrix.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_momentum_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use pcsr_matrix_type
  use block_pcsr_bridge_type
  use hypre_hybrid_type
  use flow_2d_momentum_type
  implicit none
  private

  type, public :: flow_2d_momentum_solver
    private
    type(flow_2d_momentum), pointer :: momentum => null()  ! unowned reference
    type(parameter_list), pointer :: params => null()  ! unowned reference
    type(block_pcsr_bridge), pointer :: bridge => null()
    type(hypre_hybrid) :: solver
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: get_metrics
    procedure :: metrics_string
    final :: delete
  end type

contains

  subroutine init(this, momentum, params)
    class(flow_2d_momentum_solver), intent(out) :: this
    type(flow_2d_momentum), target, intent(in) :: momentum
    type(parameter_list), target, intent(in) :: params

    this%momentum => momentum
    this%params => params
    allocate(this%bridge)
  end subroutine


  subroutine delete(this)
    type(flow_2d_momentum_solver), intent(inout) :: this

    if (associated(this%bridge)) deallocate(this%bridge)
  end subroutine


  !! Update the scalar solver representation after the block momentum matrix
  !! has been assembled.
  subroutine setup(this)
    class(flow_2d_momentum_solver), intent(inout) :: this

    type(pcsr_matrix), pointer :: matrix

    ASSERT(associated(this%momentum))
    ASSERT(associated(this%params))
    ASSERT(associated(this%bridge))
    call this%bridge%update(this%momentum%matrix())
    matrix => this%bridge%matrix()
    call this%solver%init(matrix, this%params)
    call this%solver%setup()
  end subroutine


  !! Solve for a cell-centered two-component velocity. RHS and VELOCITY use
  !! the block-interleaved storage expected by BLOCK_PCSR_BRIDGE.
  subroutine solve(this, rhs, velocity, stat)
    class(flow_2d_momentum_solver), intent(inout) :: this
    real(r8), target, contiguous, intent(in) :: rhs(:,:)
    real(r8), target, contiguous, intent(inout) :: velocity(:,:)
    integer, intent(out) :: stat

    real(r8), pointer :: b(:), x(:)

    ASSERT(size(rhs,1) == 2)
    ASSERT(size(velocity,1) == 2)
    ASSERT(size(rhs,2) == size(velocity,2))
    b(1:size(rhs)) => rhs
    x(1:size(velocity)) => velocity
    call this%solver%solve(b, x, stat)
  end subroutine


  subroutine get_metrics(this, num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
    class(flow_2d_momentum_solver), intent(in) :: this
    integer, intent(out), optional :: num_itr, num_dscg_itr, num_pcg_itr
    real(r8), intent(out), optional :: rel_res_norm

    call this%solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
  end subroutine


  function metrics_string(this) result(string)
    class(flow_2d_momentum_solver), intent(in) :: this
    character(:), allocatable :: string

    string = this%solver%metrics_string()
  end function

end module flow_2d_momentum_solver_type
