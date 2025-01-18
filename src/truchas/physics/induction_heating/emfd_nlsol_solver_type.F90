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
  use fdme_precon_class
  use fdme_ams_precon_type
  use fdme_hiptmair_precon_type
  use fdme_nlk_solver_type
  use fdme_gmres_solver_type
  use fdme_minres_solver_type
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
    class(fdme_precon), pointer :: precon => null()
    type(fdme_gmres_solver), allocatable :: gmres
    !type(fdme_minres_solver), allocatable :: minres
    type(fdme_minres_solver2), allocatable :: minres
    type(fdme_mixed_minres_solver), allocatable :: mixed_minres
    type(fdme_nlk_solver),   allocatable :: nlk
    type(fdme_mumps_solver), allocatable :: mumps

    integer :: print_level
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

    call params%get('print-level', this%print_level, stat, errmsg, default=0)
    if (stat /= 0) return
    call params%set('print-level', max(0,this%print_level-1))

    call params%get('solver-type', solver_type, stat, errmsg)
    if (stat /= 0) return

    select case (solver_type)
    case ('gmres')
      allocate(this%gmres)
    case ('minres')
      if (model%use_mixed_form) then
        allocate(this%mixed_minres)
      else
        allocate(this%minres)
      end if
    case ('nlk')
      allocate(this%nlk)
    case ('mumps')
      allocate(this%mumps)
    case default
      stat = 1
      errmsg = 'invalid solver-type value: ' // solver_type
      return
    end select

    this%model => model
    this%mesh => model%mesh

    !NB: precon is ignored by minres and mixed_minres (but can't be null)
    if (.not.allocated(this%mumps)) then
      plist => params%sublist('precon')
      call plist%get('type', choice)
      select case (choice)
      case ('ams', 'AMS')
        allocate(fdme_ams_precon :: this%precon)
        call this%precon%init(this%model, plist)
      case ('hiptmair')
        allocate(fdme_hiptmair_precon :: this%precon)
        call this%precon%init(this%model, plist)
      case ('gs', 'boomer') ! do nothing -- for minres
      case default
        call tls_fatal('unknown preconditioner type: ' // choice)
      end select
    end if

    if (allocated(this%gmres)) then
      call this%gmres%init(this%model, this%precon, params, ierr, errmsg)
    else if (allocated(this%minres)) then
      call this%minres%init(this%model, params, ierr, errmsg)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%init(this%model, this%precon, params, ierr, errmsg)
    else if (allocated(this%nlk)) then
      call this%nlk%init(this%model, this%precon, params, ierr, errmsg)
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
    if (associated(this%precon)) call this%precon%setup ! assume that the model has changed
    if (allocated(this%gmres)) then
      call this%gmres%solve(efield, stat)
      if (this%print_level > 0) then
        write(message,'(t8,a,i4,*(a,es10.3))') 'gmres solve:', this%gmres%gmres%num_iter, &
            ' iterations, |r|/|r0|=', this%gmres%gmres%rel_rnorm, &
            ', |r0|=', this%gmres%gmres%r0norm
        call TLS_info(message)
      end if
    else if (allocated(this%minres)) then
      call this%minres%solve(efield, stat)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%solve(efield, stat)
    else if (allocated(this%nlk)) then
      call this%nlk%solve(efield, stat)
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
