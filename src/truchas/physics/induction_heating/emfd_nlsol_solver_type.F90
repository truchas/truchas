!!
!! emfd_nlsol_solver_type
!!
!! This module provides a type for solving the frequency-domain Maxwell
!! Equations.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module emfd_nlsol_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fdme_model_type
  use fdme_minres_solver2_type
  use fdme_mixed_minres_solver_type
  use fdme_mumps_solver_type
  use simpl_mesh_type
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: emfd_nlsol_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference

    type(fdme_model), pointer :: model => null()
    type(fdme_minres_solver2), allocatable :: minres
    type(fdme_mixed_minres_solver), allocatable :: mixed_minres
    type(fdme_mumps_solver), allocatable :: mumps

    integer :: print_level
  contains
    procedure :: init
    procedure :: solve
  end type emfd_nlsol_solver

contains

  subroutine init(this, model, params, stat, errmsg)

    class(emfd_nlsol_solver), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    integer :: ierr
    character(:), allocatable :: solver_type, choice

    call params%get('print-level', this%print_level, stat, errmsg, default=0)
    if (stat /= 0) return
    call params%set('print-level', max(0,this%print_level-1))

    call params%get('solver-type', solver_type, stat, errmsg)
    if (stat /= 0) return

    select case (solver_type)
    case ('minres')
      if (model%use_mixed_form) then
        allocate(this%mixed_minres)
      else
        allocate(this%minres)
      end if
    case ('mumps')
      allocate(this%mumps)
    case default
      stat = 1
      errmsg = 'invalid solver-type value: ' // solver_type
      return
    end select

    this%model => model
    this%mesh => model%mesh

    if (allocated(this%minres)) then
      call this%minres%init(this%model, params, ierr, errmsg)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%init(this%model, params, ierr, errmsg)
    else if (allocated(this%mumps)) then
      call this%mumps%init(this%model, params, ierr, errmsg)
    else
      INSIST(.false.)
    end if
    if (ierr /= 0) call tls_fatal("EMFD_NLSOL INIT: " // errmsg)

  end subroutine init


  subroutine solve(this, efield, stat, errmsg)

    class(emfd_nlsol_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    character(72) :: message

    call start_timer("solve")
    if (allocated(this%minres)) then
      call this%minres%solve(efield, stat)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%solve(efield, stat)
    else if (allocated(this%mumps)) then
      call this%mumps%solve(efield, stat)
    else
      INSIST(.false.)
    end if
    if (stat /= 0) errmsg = 'FDME convergence failure'

    block
      use parallel_communication, only: global_norm2
      complex(r8) :: r(this%mesh%nedge), div_efield(this%mesh%nnode)
      real(r8) :: r_norm2, d_norm2
      character(128) :: msg
      call this%model%residual(efield, r)
      call this%model%compute_div(efield, div_efield)
      r_norm2 = global_norm2(r(:this%mesh%nedge_onP))
      d_norm2 = global_norm2(div_efield(:this%mesh%nnode_onP))
      write (msg,'(a,2es14.4)') "EMFD solve complete. Residual, divE: ", r_norm2, d_norm2
      call tls_info(msg)
      !INSIST(.false.)
    end block

    call stop_timer("solve")
  end subroutine

end module emfd_nlsol_solver_type
