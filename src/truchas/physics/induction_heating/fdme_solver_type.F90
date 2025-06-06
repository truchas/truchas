!!
!! FDME_SOLVER_TYPE
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

module fdme_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use fdme_model_type
  use fdme_minres_solver_type
  use fdme_mixed_minres_solver_type
  use fdme_mumps_solver_type
  implicit none
  private

  type, public :: fdme_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(fdme_model), pointer :: model => null()
    type(fdme_minres_solver), allocatable :: minres
    type(fdme_mixed_minres_solver), allocatable :: mixed_minres
    type(fdme_mumps_solver), allocatable :: mumps
    integer :: print_level
  contains
    procedure :: init
    procedure :: solve
  end type fdme_solver

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(fdme_solver), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
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
      call this%minres%init(this%model, params, stat, errmsg)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%init(this%model, params, stat, errmsg)
    else if (allocated(this%mumps)) then
      call this%mumps%init(this%model, params, stat, errmsg)
    else
      INSIST(.false.)
    end if

  end subroutine init


  subroutine solve(this, efield, stat, errmsg)

    use string_utilities, only: i_to_c

    class(fdme_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    if (allocated(this%minres)) then
      call this%minres%solve(efield, stat, errmsg)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%solve(efield, stat, errmsg)
    else if (allocated(this%mumps)) then
      call this%mumps%solve(efield, stat)
      if (stat /= 0) errmsg = 'MUMPS solve failed: stat=' // i_to_c(stat)
    else
      INSIST(.false.)
    end if

    block
      use parallel_communication, only: global_norm2
      use truchas_logging_services
      complex(r8) :: r(this%mesh%nedge), div_efield(this%mesh%nnode)
      real(r8) :: rnorm, bnorm, dnorm
      character(128) :: msg
      call this%model%residual(efield, r)
      call this%model%compute_div(efield, div_efield)
      rnorm = global_norm2(r(:this%mesh%nedge_onP))
      bnorm = global_norm2(this%model%rhs(:this%mesh%nedge_onP))
      dnorm = global_norm2(div_efield(:this%mesh%nnode_onP))
      write(msg,'(*(a,es8.2))') 'FDME_SOLVE: |r|=', rnorm, &
          ', |b|=', bnorm, ', |r|/|b|=', rnorm/bnorm, ', |Div E|=', dnorm
      call TLS_info(msg)
    end block

  end subroutine

end module fdme_solver_type
