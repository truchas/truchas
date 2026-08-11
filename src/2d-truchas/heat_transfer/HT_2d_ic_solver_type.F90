!! HT_2D_IC_SOLVER_TYPE
!!
!! Constructs a consistent initial state for the 2D heat-transfer DAE.

#include "f90_assert.fpp"

module HT_2d_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use HT_2d_model_type
  implicit none
  private

  type, public :: HT_2d_ic_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(HT_2d_model), pointer :: model => null()
  contains
    procedure :: init
    procedure :: compute
  end type HT_2d_ic_solver

contains

  subroutine init(this, model)
    class(HT_2d_ic_solver), intent(out) :: this
    type(HT_2d_model), intent(in), target :: model
    this%model => model
    this%mesh => model%mesh
  end subroutine init


  subroutine compute(this, t, dt, temp, u, udot, rel_tol, max_itr, stat, errmsg)
    class(HT_2d_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, temp(:), rel_tol
    real(r8), intent(out), target :: u(:), udot(:)
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: state(:,:), Hcell(:)
    real(r8), pointer :: ucell(:), uface(:)

    stat = 0
    errmsg = ''

    ASSERT(size(temp) == this%mesh%ncell_onP)
    ASSERT(size(u) == this%model%num_dof())
    ASSERT(size(udot) == size(u))

    call this%model%get_cell_temp_view(u, ucell)
    call this%model%get_face_temp_view(u, uface)
    ucell = temp

    allocate(Hcell(this%mesh%ncell))
    call this%model%new_state_array(u, state)
    call this%model%H_of_T%compute_value(state, Hcell)
    call this%model%set_cell_heat(Hcell, u)
    deallocate(state, Hcell)

    !! The face solve obtains the consistent algebraic state. The zero face
    !! field is only an initial guess for that solve.
    uface = 0.0_r8
    call this%model%compute_face_temp(t, u, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) return

    call this%model%compute_udot(t, dt, u, udot, rel_tol, max_itr, stat, errmsg)
  end subroutine compute

end module HT_2d_ic_solver_type
