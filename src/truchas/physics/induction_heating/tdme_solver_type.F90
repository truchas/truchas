!!
!! TDME_SOLVER_TYPE
!!
!! This module provides a derived type that encapsulates time stepping of
!! the discrete time-domain Maxwell equations.
!!
!! Neil Carlson
!! Refactored February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module tdme_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use tdme_model_type
  use tdme_cg_solver_type
  use state_history_type
  use parallel_communication
  use truchas_logging_services
  implicit none
  private

  type, public :: tdme_solver
    type(simpl_mesh), pointer :: mesh  => null() ! unowned reference
    type(tdme_model), pointer :: model => null() ! unowned reference
    real(r8) :: t, dt
    type(state_history) :: ehist, bhist
    real(r8), allocatable :: r(:), de(:) ! persistent workspace
    type(tdme_cg_solver) :: eslv
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: step
  end type

contains

  subroutine init(this, mesh, model, params, stat, errmsg)
    use parameter_list_type
    class(tdme_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(tdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%mesh => mesh
    this%model => model
    allocate(this%r(this%mesh%nedge), this%de(this%mesh%nedge))
    call this%eslv%init(this%model, params, stat, errmsg)
  end subroutine

  subroutine set_initial_state(this, t, e, b)
    class(tdme_solver), intent(inout) :: this
    real(r8), intent(in) :: t, e(:), b(:)
    ASSERT(size(e) == this%mesh%nedge)
    ASSERT(size(b) == this%mesh%nface)
    this%t = t
    this%dt = this%model%dt
    call this%ehist%init(3, t, e)
    call this%bhist%init(1, t, b)  !NB: will need mvec > 1 for interpolated output
  end subroutine


  subroutine step(this, t, e, b, stat, errmsg)

    class(tdme_solver), intent(inout), target :: this
    real(r8), intent(out) :: t, e(:), b(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), pointer :: e0(:), b0(:)
    character(80) :: message

    ASSERT(size(e) == this%mesh%nedge)
    ASSERT(size(b) == this%mesh%nface)

    t = this%t + this%dt

    call this%ehist%get_last_state_view(e0)
    call this%bhist%get_last_state_view(b0)

    call this%ehist%interp_state(t, e) ! initial guess: quadratic extrapolation
    call this%model%trap_res(this%t, e0, b0, e, this%r)

    this%de = 0.0_r8 ! initial guess for the correction
    call this%eslv%solve(this%r, this%de, this%mesh%nedge_onP, stat, errmsg)
    if (stat /= 0) return
    call this%mesh%edge_imap%gather_offp(this%de)
    e = e + this%de

    call this%model%advance_bfield(e0, b0, e, b)

    this%t = t
    call this%ehist%record_state(t, e)
    call this%bhist%record_state(t, b)

    if (this%eslv%output_level >= 3) then
      write(message,fmt='(t6,a,es10.3)') 'step |de|_max/|e|_max=', &
          global_maxval(abs(this%de))/global_maxval(abs(e))
      call TLS_info(trim(message))
    end if

  end subroutine step

end module tdme_solver_type
