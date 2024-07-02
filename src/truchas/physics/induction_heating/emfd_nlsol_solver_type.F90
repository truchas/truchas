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

!#define SOLVER fdme_gmres_solver
#define SOLVER fdme_nlk_solver

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
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: emfd_nlsol_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference

    class(fdme_precon), pointer :: precon => null()
    type(SOLVER) :: solver
    type(fdme_model), pointer :: model => null()

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

    use em_bc_factory_type
    use parameter_list_type

    class(emfd_nlsol_solver), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    !type(simpl_mesh), intent(in), target :: mesh
    !type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    integer :: ierr
    real(r8) :: atol, rtol, ftol
    character(:), allocatable :: choice

    ! set defaults
    !if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 0.0_r8)
    !if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 1d-1)
    !if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 5d-4) ! gmres-hiptmair
    !if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 1d-8) ! gmres-hiptmair HF-ND
    if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 1d-13) ! gmres-hiptmair HF-ND
    !if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 1d-16) ! gmres-hiptmair HF-ND
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 2d-6)
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 3d-7) ! gmres-ams
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 3d-7) ! nlsol-ams
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 1d-8) ! gmres-hiptmair
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 1d-11) ! gmres-hiptmair HF
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 0.0_r8) ! gmres-hiptmair HF
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 0.0_r8) ! gmres-hiptmair HF
    if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 1d-8) ! nlsol-hiptmair
    if (.not.params%is_parameter("res-tol")) call params%set("res-tol", huge(1.0_r8))
    !if (.not.params%is_parameter("res-tol")) call params%set("res-tol", 1d-1) ! gmres-hiptmair HF
    !if (.not.params%is_parameter("res-tol")) call params%set("res-tol", 1d0) ! nlsol-hiptmair
    !if (.not.params%is_parameter("nlk-tol")) call params%set("nlk-tol", 1d-1)
    if (.not.params%is_parameter("nlk-tol")) call params%set("nlk-tol", 1d-8)
    if (.not.params%is_parameter("nlk-max-iter")) call params%set("nlk-max-iter", 4000)
    if (.not.params%is_parameter("verbosity")) call params%set("verbosity", 2)

    ! use nlk_solver defaults for the present
    !if (.not.params%is_parameter('nlk-abs-tol')) call params%set('nlk-abs-tol', 0.0_r8)
    !if (.not.params%is_parameter('nlk-rel-tol')) call params%set('nlk-rel-tol', 0.0_r8)

    !call params%set("nlk-max-vec", 0)
    !call params%set("gmres-krylov-dim", 20)
    !call params%set("gmres-krylov-dim", 2)

    call params%set("max-iter", 2000)

    ! nlsol parameters
    if (.not.params%is_parameter("nlk-max-vec")) call params%set("nlk-max-vec", 20)
    !if (.not.params%is_parameter("nlk-vec-tol")) call params%set("nlk-vec-tol", 1d-2)

    this%model => model
    !allocate(this%model)
    !call this%model%init(mesh, bc_fac, params, stat, errmsg)
    !if (stat /= 0) return

    this%mesh => model%mesh
    call this%efield%init(this%mesh)

    call params%get("abs-tol", atol)
    call params%get("rel-tol", rtol)
    call params%get("res-tol", ftol)

    plist => params%sublist("precon")
    call plist%set('type', 'hiptmair')
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

    call this%solver%init(this%efield, this%model, this%precon, params, ierr, errmsg)
    if (ierr /= 0) call tls_fatal("EMFD_NLSOL INIT: " // errmsg)

  end subroutine init


  subroutine solve(this, efield)

    class(emfd_nlsol_solver), intent(inout) :: this
    type(fdme_vector), intent(inout) :: efield

    integer :: ierr

    print '(a,2es13.3)', "max |rhs| = ", maxval(abs(this%model%rhs)), maxval(abs(this%model%hbc%value))

    ASSERT(all(ieee_is_finite(this%efield%array)))
    ASSERT(all(ieee_is_finite(this%model%rhs)))

    call start_timer("solve")
    call this%solver%solve(this%efield, ierr)
    call stop_timer("solve")
    !call tls_info('  EMFD solve: ' // this%solver%metrics_string())
    print *, "ierr: ", ierr
    if (ierr /= 0) call tls_error("EMFD solve unsuccessful")

    print *, "max |rhs| = ", maxval(abs(this%model%rhs)), maxval(abs(this%model%hbc%value))
    print *, "max |er| = ", maxval(abs(this%efield%array(1,:)))
    print *, "max |ei| = ", maxval(abs(this%efield%array(2,:)))
    print *, "max |e-rhs| = ", maxval(abs(this%efield%array - this%model%rhs))

    call efield%copy(this%efield)
  end subroutine solve

end module emfd_nlsol_solver_type
