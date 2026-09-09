!!
!! T2D_FLOW_PROJECTION_SOLVER_TYPE
!!
!! This module defines T2D_FLOW_PROJECTION_SOLVER, the linear-solver adapter
!! for the scalar pressure system stored by T2D_FLOW_PROJECTION. It owns the
!! HYPRE_HYBRID solver and retains unowned references to the projection
!! operator and its solver parameter list.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flow_projection_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use t2d_flow_projection_type
  use hypre_hybrid_type
  implicit none
  private

  type, public :: t2d_flow_projection_solver
    private
    type(t2d_flow_projection), pointer :: projection => null()  ! unowned reference
    type(parameter_list), pointer :: params => null()  ! unowned reference
    type(hypre_hybrid) :: solver
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: get_metrics
    procedure :: metrics_string
  end type

contains

  subroutine init(this, projection, params, stat, errmsg)
    class(t2d_flow_projection_solver), intent(out) :: this
    type(t2d_flow_projection), target, intent(in) :: projection
    type(parameter_list), target, intent(in) :: params
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    call this%solver%validate_params(params, stat, errmsg)
    if (present(stat)) then
      if (stat /= 0) then
        if (present(errmsg)) errmsg = 'processing ' // params%path() // ': ' // errmsg
        return
      end if
    end if
    this%projection => projection
    this%params => params
  end subroutine


  !! Update the HYPRE solver after the pressure matrix has been assembled.
  subroutine setup(this)
    class(t2d_flow_projection_solver), intent(inout) :: this

    ASSERT(associated(this%projection))
    ASSERT(associated(this%params))
    call this%solver%init(this%projection%matrix(), this%params)
    call this%solver%setup()
  end subroutine


  subroutine solve(this, rhs, pressure, stat)
    class(t2d_flow_projection_solver), intent(inout) :: this
    real(r8), intent(in) :: rhs(:)
    real(r8), intent(inout) :: pressure(:)
    integer, intent(out) :: stat

    call this%solver%solve(rhs, pressure, stat)
  end subroutine


  subroutine get_metrics(this, num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
    class(t2d_flow_projection_solver), intent(in) :: this
    integer, intent(out), optional :: num_itr, num_dscg_itr, num_pcg_itr
    real(r8), intent(out), optional :: rel_res_norm

    call this%solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, rel_res_norm)
  end subroutine


  function metrics_string(this) result(string)
    class(t2d_flow_projection_solver), intent(in) :: this
    character(:), allocatable :: string

    string = this%solver%metrics_string()
  end function

end module t2d_flow_projection_solver_type
