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

#define PRECON ams_precon
#define SOLVER gmres_left_solver
!#define PRECON hiptmair_precon
!#define SOLVER nlk_solver

module emfd_nlsol_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use vector_class
  use fdme_vector_type
  use fdme_model_type
  use ams_precon_type
  use hiptmair_precon_type
  use nlk_solver_type
  use gmres_left_solver_type
  use pcsr_matrix_type
  use simpl_mesh_type
  use index_map_type
  use bndry_func1_class
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, extends(nlk_solver_model) :: emfd_nlsol_model
    private
    type(fdme_model), pointer :: model
    real(r8) :: atol, rtol, ftol
    real(r8), allocatable :: b(:,:)
    !type(pcsr_matrix), pointer :: A(:,:) => null() ! unowned reference
    type(PRECON), pointer :: precon => null() ! unowned reference
  contains
    procedure :: init => emfd_nlsol_model_init
    !!!! INHERITED !!!!
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: is_converged
  end type emfd_nlsol_model


  type, public :: emfd_nlsol_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference

    type(PRECON), pointer :: precon => null()
    type(SOLVER) :: solver
    type(emfd_nlsol_model), pointer :: model => null()
    type(fdme_model), pointer :: newmodel => null()
    type(pcsr_matrix), pointer :: Ap => null()
    !type(pcsr_matrix), pointer :: A(:,:) => null(), Ap => null()

    type(fdme_vector) :: efield ! electric field (real and imaginary parts)
    !real(r8), allocatable :: rhs(:) ! rhs (set by BCs)
    !real(r8), allocatable :: rhs_r(:), rhs_i(:)

    !real(r8), allocatable :: epsr(:)  ! real part of the permittivity
    !real(r8), allocatable :: epsi(:)  ! imaginary part of the permittivity
    !real(r8), allocatable :: mu(:)    ! permeability
    !real(r8), allocatable :: sigma(:) ! conductivity
    !real(r8) :: omega                 ! driving angular frequency

    ! Non-dimensionalization parameters
    real(r8) :: Z0 = 376.730313668 ! Ohms -- vacuum impedance
    real(r8) :: c0 = 299792458.0_r8 ! speed of light
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: compute_heat_source
    final :: emfd_nlsol_solver_delete
  end type emfd_nlsol_solver

contains

  !! Final subroutine for emfd_nlsol_solver type objects.
  subroutine emfd_nlsol_solver_delete(this)
    type(emfd_nlsol_solver), intent(inout) :: this
    if (associated(this%Ap)) deallocate(this%Ap)
    if (associated(this%precon)) deallocate(this%precon)
    if (associated(this%model)) deallocate(this%model)
  end subroutine emfd_nlsol_solver_delete


  subroutine init(this, mesh, bc_fac, params, stat, errmsg)

    use em_bc_factory_type
    use parameter_list_type

    class(emfd_nlsol_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    integer :: ierr
    real(r8) :: atol, rtol, ftol

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
    if (.not.params%is_parameter("nlk-tol")) call params%set("nlk-tol", 1d-1)
    if (.not.params%is_parameter("nlk-max-iter")) call params%set("nlk-max-iter", 4000)
    if (.not.params%is_parameter("verbosity")) call params%set("verbosity", 2)

    !call params%set("nlk-max-vec", 0)
    !call params%set("gmres-krylov-dim", 20)
    !call params%set("gmres-krylov-dim", 2)

    call params%set("max-iter", 2000)

    ! nlsol parameters
    if (.not.params%is_parameter("nlk-max-vec")) call params%set("nlk-max-vec", 20)
    !if (.not.params%is_parameter("nlk-vec-tol")) call params%set("nlk-vec-tol", 1d-2)

    this%mesh => mesh
    ! this%epsr = epsr
    ! this%epsi = epsi
    ! this%mu = mu
    ! this%sigma = sigma
    ! this%omega = omega
    !this%emask = emask
    call this%efield%init(mesh)

    allocate(this%newmodel)
    call this%newmodel%init(mesh, bc_fac, params, stat, errmsg)
    if (stat /= 0) return

!    allocate(this%rhs_r(mesh%nedge), this%rhs_i(mesh%nedge), source=0.0_r8)
!
!    block ! full system in block form
!      type(pcsr_graph), pointer :: g
!      allocate(g)
!      call g%init(this%mesh%edge_imap)
!      call g%add_clique(this%mesh%cedge)
!      call g%add_complete
!      allocate(this%A(2,2), this%Ap)
!      call this%A(1,1)%init(g, take_graph=.true.)
!      call this%A(1,2)%init(mold=this%A(1,1))
!      call this%A(2,1)%init(mold=this%A(1,1))
!      call this%A(2,2)%init(mold=this%A(1,1))
!      call this%Ap%init(mold=this%A(1,1))
!    end block

!    allocate(this%Ap)
!    call this%Ap%init(mold=this%newmodel%A(1,1))

    call params%get("abs-tol", atol)
    call params%get("rel-tol", rtol)
    call params%get("res-tol", ftol)

    allocate(this%precon, this%model)
    plist => params%sublist("precon")
    call this%precon%init(this%newmodel, plist)
    call this%model%init(this%newmodel, this%precon, atol, rtol, ftol)
    call this%solver%init(this%model, params, ierr, errmsg)
    if (ierr /= 0) call tls_fatal("EMFD_NLSOL INIT: " // errmsg)

  end subroutine init


  subroutine setup(this, t, epsr, epsi, mu, sigma, omega)

    class(emfd_nlsol_solver), intent(inout) :: this
    real(r8), intent(in) :: t, epsr(:), epsi(:), mu(:), sigma(:), omega

    call this%newmodel%setup(t, epsr, epsi, mu, sigma, omega)

    select type (rhs => this%model%rhs)
    type is (fdme_vector)
      rhs%array(:,:) = this%newmodel%rhs
    end select

!    ASSERT(all(ieee_is_finite(this%model%rhs)))

!    ! Preconditioner
!    this%Ap%values = this%newmodel%A(1,1)%values - this%newmodel%A(1,2)%values
!    if (allocated(this%newmodel%ebc)) then
!      do j = 1, size(this%newmodel%ebc%index)
!        n = this%newmodel%ebc%index(j)
!        call this%Ap%set(n, n, 1.0_r8)
!      end do
!    end if
    call this%precon%setup ! this implicitly uses this%Ap

  end subroutine setup


  subroutine solve(this)

    class(emfd_nlsol_solver), intent(inout) :: this

    integer :: ierr

    print '(a,2es13.3)', "max |rhs| = ", maxval(abs(this%newmodel%rhs)), maxval(abs(this%newmodel%hbc%value))

    ASSERT(all(ieee_is_finite(this%efield%array)))
    ASSERT(all(ieee_is_finite(this%newmodel%rhs)))

    call start_timer("solve")
    call this%solver%solve(this%efield, ierr)
    call stop_timer("solve")
    !call tls_info('  EMFD solve: ' // this%solver%metrics_string())
    print *, "ierr: ", ierr
    if (ierr /= 0) call tls_error("EMFD solve unsuccessful")

    print *, "max |rhs| = ", maxval(abs(this%newmodel%rhs)), maxval(abs(this%newmodel%hbc%value))
    print *, "max |er| = ", maxval(abs(this%efield%array(1,:)))
    print *, "max |ei| = ", maxval(abs(this%efield%array(2,:)))
    print *, "max |e-rhs| = ", maxval(abs(this%efield%array - this%newmodel%rhs))

  end subroutine solve


  subroutine compute_heat_source(this, q)

    class(emfd_nlsol_solver), intent(in), target :: this
    real(r8), intent(out) :: q(:)

    call start_timer("heat source")

    call this%newmodel%compute_heat_source(this%efield%array, q)

    call stop_timer("heat source")

  end subroutine compute_heat_source


!!!! MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine emfd_nlsol_model_init(this, model, precon, atol, rtol, ftol)
    class(emfd_nlsol_model), intent(out) :: this
    !type(pcsr_matrix), pointer :: A(:,:)
    type(fdme_model), pointer :: model
    type(PRECON), pointer :: precon
    real(r8), intent(in) :: atol, rtol, ftol
    !this%A => A
    this%model => model
    this%precon => precon
    this%atol = atol
    this%rtol = rtol
    this%ftol = ftol
    block
      type(fdme_vector), allocatable :: tmp
      allocate(tmp)
      call tmp%init(model%mesh)
      call move_alloc(tmp, this%rhs)
    end block
  end subroutine emfd_nlsol_model_init

  subroutine alloc_vector(this, u)
    class(emfd_nlsol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: u
    type(fdme_vector), allocatable :: tmp
    allocate(tmp)
    call this%model%init_vector(tmp)
    call move_alloc(tmp, u)
  end subroutine

  subroutine compute_f(this, u, f, ax)
    class(emfd_nlsol_model) :: this
    class(vector), intent(in) :: u
    class(vector), intent(inout) :: f
    logical, intent(in), optional :: ax
    select type (u)
    type is (fdme_vector)
      select type (f)
      type is (fdme_vector)
        call this%model%compute_f(u%array, f%array, ax)
      end select
    end select
  end subroutine


  subroutine apply_precon(this, u, f)
    class(emfd_nlsol_model) :: this
    class(vector), intent(in) :: u  ! not used
    class(vector), intent(inout) :: f
    select type (u)
    type is (fdme_vector)
      select type (f)
      type is (fdme_vector)
        this%b = f%array  !this allocates b
        f%array = 0.0_r8
        ASSERT(all(ieee_is_finite(this%b)))
        call this%precon%apply(this%b, f%array)
        ASSERT(all(ieee_is_finite(f%array)))
      end select
    end select
  end subroutine


  subroutine compute_precon(this, u)
    class(emfd_nlsol_model) :: this
    class(vector), intent(in) :: u
    ! compute precon currently handled above the nlsol solver.
    !call this%precon%setup(interior_nodes)
  end subroutine


  real(r8) function du_norm(this, u, du)

    use parallel_communication, only: global_maxval

    class(emfd_nlsol_model) :: this
    class(vector), intent(in) :: u, du

    integer :: i, j
    real(r8) :: err

    select type (u)
    type is (fdme_vector)
      select type (du)
      type is (fdme_vector)

        du_norm = 0.0_r8
        do j = 1, this%model%mesh%nedge_onP
          do i = 1, 2
            if (this%atol == 0 .and. abs(u%array(i,j)) == 0) then
              err = huge(1.0_r8)
            else
              err = abs(du%array(i,j)) / (this%atol + this%rtol*abs(u%array(i,j)))
            end if
            du_norm = max(du_norm, err)
          end do
        end do
        du_norm = global_maxval(du_norm)

      end select
    end select
  end function


  logical function is_converged(this, itr, u, du, f_lnorm, tol)
    class(emfd_nlsol_model) :: this
    integer, intent(in) :: itr
    class(vector), intent(in) :: u, du
    real(r8), intent(in) :: tol, f_lnorm(:)
    is_converged = this%du_norm(u, du) < tol .and. f_lnorm(3) < this%ftol
  end function

end module emfd_nlsol_solver_type
