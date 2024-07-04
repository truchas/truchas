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
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use fdme_vector_type
  use fdme_model_type
  use fdme_precon_class
  use fdme_ams_precon_type
  use fdme_hiptmair_precon_type
  use fdme_nlk_solver_type
  use fdme_gmres_solver_type
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
    class(fdme_precon), pointer :: precon => null()
    type(fdme_gmres_solver), allocatable :: gmres
    type(fdme_nlk_solver),   allocatable :: nlk

    type(fdme_vector) :: efield ! electric field (real and imaginary parts)
  contains
    procedure :: init
    procedure :: solve
    final :: emfd_nlsol_solver_delete
  end type emfd_nlsol_solver

contains

  !! Final subroutine for emfd_nlsol_solver type objects.
  subroutine emfd_nlsol_solver_delete(this)
    type(emfd_nlsol_solver), intent(inout) :: this
    if (associated(this%precon)) deallocate(this%precon)
  end subroutine emfd_nlsol_solver_delete


  subroutine init(this, model, params, stat, errmsg)

    class(emfd_nlsol_solver), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    integer :: ierr
    character(:), allocatable :: solver_type, choice

    call params%get('solver-type', solver_type, stat, errmsg)
    if (stat /= 0) return

    select case (solver_type)
    case ('gmres')
      allocate(this%gmres)
    case ('nlk')
      allocate(this%nlk)
    case default
      stat = 1
      errmsg = 'invalid solver-type value: ' // solver_type
      return
    end select

    this%model => model

    this%mesh => model%mesh
    call this%efield%init(this%mesh)

    plist => params%sublist('precon')
    call plist%get('type', choice)
    select case (choice)
    case ('ams', 'AMS')
      allocate(fdme_ams_precon :: this%precon)
    case ('hiptmair')
      allocate(fdme_hiptmair_precon :: this%precon)
    case default
      call tls_fatal('unknown preconditioner type: ' // choice)
    end select
    call this%precon%init(this%model, plist)

    if (allocated(this%gmres)) then
      call this%gmres%init(this%efield, this%model, this%precon, params, ierr, errmsg)
    else if (allocated(this%nlk)) then
      call this%nlk%init(this%efield, this%model, this%precon, params, ierr, errmsg)
    else
      INSIST(.false.)
    end if
    if (ierr /= 0) call tls_fatal("EMFD_NLSOL INIT: " // errmsg)

  end subroutine init


  subroutine solve(this, efield, stat, errmsg)
    class(emfd_nlsol_solver), intent(inout) :: this
    type(fdme_vector), intent(inout) :: efield
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%precon%setup ! assume that the model has changed
    if (allocated(this%gmres)) then
      call this%gmres%solve(efield, stat)
    else if (allocated(this%nlk)) then
      call this%nlk%solve(efield, stat)
    else
      INSIST(.false.)
    end if
    if (stat /= 0) errmsg = 'FDME convergence failure'
  end subroutine

end module emfd_nlsol_solver_type
