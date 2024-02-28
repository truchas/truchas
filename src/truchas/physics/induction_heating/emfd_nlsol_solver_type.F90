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

!#define PRECON hiptmair_precon
#define PRECON ams_precon
#define SOLVER gmres_left_solver
!#define SOLVER nlsol

module emfd_nlsol_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  !use,intrinsic :: iso_fortran_env, only: r8 => real128
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ams_precon_type
  use hiptmair_precon_type
  use nlsol_type
  use gmres_left_solver_type
  use pcsr_matrix_type
  use simpl_mesh_type
  use index_map_type
  use em_bc_type
  use bndry_func1_class
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, extends(nlsol_model) :: emfd_nlsol_model
    private
    real(r8) :: atol, rtol, ftol
    real(r8), allocatable :: b(:)
    type(pcsr_matrix), pointer :: A => null() ! unowned reference
    type(PRECON), pointer :: precon => null() ! unowned reference
  contains
    procedure :: init => emfd_nlsol_model_init
    !!!! INHERITED !!!!
    procedure :: size => emfd_nlsol_model_size
    procedure :: compute_f => emfd_nlsol_model_compute_f
    procedure :: apply_precon => emfd_nlsol_model_apply_precon
    procedure :: compute_precon => emfd_nlsol_model_compute_precon
    procedure :: du_norm => emfd_nlsol_model_du_norm
    procedure :: is_converged => emfd_nlsol_model_is_converged
  end type emfd_nlsol_model


  type, public :: emfd_nlsol_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(simpl_mesh) :: meshr

    type(PRECON), pointer :: precon => null()
    type(SOLVER) :: solver
    type(emfd_nlsol_model), pointer :: model => null()
    type(index_map), pointer :: imap => null()
    type(pcsr_matrix), pointer :: A => null(), Ap => null()

    type(em_bc) :: bc

    real(r8), allocatable :: efield(:) ! electric field (real and imaginary parts)
    real(r8), allocatable :: rhs(:) ! rhs (set by BCs)

    real(r8), allocatable :: epsr(:)  ! real part of the permittivity
    real(r8), allocatable :: epsi(:)  ! imaginary part of the permittivity
    real(r8), allocatable :: mu(:)    ! permeability
    real(r8), allocatable :: sigma(:) ! conductivity
    real(r8) :: omega                 ! driving angular frequency

    ! Non-dimensionalization parameters
    real(r8) :: Z0 = 376.730313668 ! Ohms -- vacuum impedance
    real(r8) :: c0 = 299792458.0_r8 ! speed of light
    real(r8) :: L0, H0
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
    if (associated(this%imap)) deallocate(this%imap)
    if (associated(this%A)) deallocate(this%A)
    if (associated(this%Ap)) deallocate(this%Ap)
    if (associated(this%precon)) deallocate(this%precon)
    if (associated(this%model)) deallocate(this%model)
  end subroutine emfd_nlsol_solver_delete


  subroutine init(this, mesh, bc_fac, params, stat, errmsg)

    use mimetic_discretization
    use em_bc_factory_type
    use parameter_list_type
    use parallel_communication, only: is_IOP

    class(emfd_nlsol_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    type(index_map), pointer :: row_imap
    integer, allocatable :: nvars(:)
    integer, allocatable :: ebedge(:)
    integer :: e, c, er, ei
    integer :: e1, e1x, e1r, e1i, e2, e2x, e2r, e2i
    integer :: ierr
    real(r8) :: atol, rtol, ftol

    block
      class(bndry_func1), allocatable :: ebc, hbc
      call bc_fac%alloc_nxE_bc(ebc, stat, errmsg)
      if (stat /= 0) return
      call bc_fac%alloc_fd_nxH_bc(hbc, stat, errmsg, omit_edge_list=ebc%index)
      if (stat /= 0) return
      call this%bc%init(mesh, ebc, hbc)
    end block

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
    if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 0.0_r8) ! gmres-hiptmair HF
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 0.0_r8) ! gmres-hiptmair HF
    !if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 1d-8) ! nlsol-hiptmair
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
    allocate(this%efield(2*mesh%nedge), this%rhs(2*mesh%nedge))
    this%efield = 0
    this%rhs = 0

    call init_full_system(this%A)
    call init_precon_system(this%Ap)

    call params%get("abs-tol", atol)
    call params%get("rel-tol", rtol)
    call params%get("res-tol", ftol)

    ! non-dimensionalization
    this%L0 = 1
    this%H0 = 1 !maxval(abs(this%bc%hsource))
    this%meshr = this%mesh
    this%meshr%x = this%meshr%x / this%L0
    this%meshr%length = this%meshr%length / this%L0
    this%meshr%volume = this%meshr%volume / this%L0**3

    allocate(this%precon, this%model)
    plist => params%sublist("precon")
    call this%precon%init(plist, this%meshr, this%Ap, this%bc)
    call this%model%init(this%A, this%precon, atol, rtol, ftol)
    call this%solver%init(this%model, params, ierr, errmsg)
    if (ierr /= 0) call tls_fatal("EMFD_NLSOL INIT: " // errmsg)

  contains

    subroutine init_full_system(A)

      type(pcsr_matrix), pointer :: A
      type(pcsr_graph), pointer :: g

      allocate(g, this%imap, A)
      row_imap => mesh%edge_imap
      allocate(nvars(merge(this%mesh%edge_imap%global_size, 0, is_IOP)))
      nvars = 2
      call this%imap%init(row_imap, nvars)

      call g%init(this%imap)
      do c = 1, mesh%ncell
        ! add curl terms
        call g%add_clique(2*(mesh%cedge(:,c)-1) + 1) ! real components
        call g%add_clique(2*(mesh%cedge(:,c)-1) + 2) ! imaginary components

        ! add cross-diagonal terms
        do e1x = 1, size(mesh%cedge(:,c))
          e1 = mesh%cedge(e1x,c)
          e1r = 2*(e1-1) + 1
          e1i = 2*(e1-1) + 2
          do e2x = 1, e1x
            e2 = mesh%cedge(e2x,c)
            e2r = 2*(e2-1) + 1
            e2i = 2*(e2-1) + 2
            call g%add_edge(e1r, e2i)
            call g%add_edge(e1i, e2r)
          end do
        end do
      end do
      call g%add_complete
      call A%init(g, take_graph=.true.)

    end subroutine init_full_system

    subroutine init_precon_system(A)
      type(pcsr_matrix), pointer :: A
      type(pcsr_graph), pointer :: g
      allocate(g, A)
      call g%init(this%mesh%edge_imap)
      do c = 1, mesh%ncell
        call g%add_clique(mesh%cedge(:,c))
      end do
      call g%add_complete
      call A%init(g, take_graph=.true.)
    end subroutine init_precon_system

  end subroutine init


  subroutine setup(this, t, epsr, epsi, mu, sigma, omega)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we
    use upper_packed_matrix_procs, only: sym_matmul

    class(emfd_nlsol_solver), intent(inout) :: this
    real(r8), intent(in) :: t, epsr(:), epsi(:), mu(:), sigma(:), omega

    real(r8), parameter :: curl(4,6) = reshape(source=[&
        0.0_r8,  0.0_r8,  1.0_r8,  1.0_r8, &
        0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
        0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
        1.0_r8,  0.0_r8,  0.0_r8,  1.0_r8, &
        -1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, &
        1.0_r8,  1.0_r8,  0.0_r8,  0.0_r8], &
        shape=shape(curl))

    integer :: i, j, k, l, m, n, e1, e1x, e1r, e1i, e2, e2x, e2r, e2i
    real(r8) :: ctm2(4), ctm2c(6,6), a0(4,4), m1(21), m2(10), mtr1, mtr2, mtr3, omegar

    ASSERT(size(epsr) == this%mesh%ncell)
    ASSERT(size(epsi) == this%mesh%ncell)
    ASSERT(size(mu) == this%mesh%ncell)
    ASSERT(size(sigma) == this%mesh%ncell)

    print *, "epsr: ", epsr(1)
    print *, "epsi: ", epsi(1)
    print *, "mu: ", mu(1)
    print *, "sigma: ", sigma(1)
    print *, "omega: ", omega

    call start_timer("setup")

    this%sigma = sigma
    this%epsi = epsi

    ! non-dimensionalization
    omegar = omega * this%L0 / this%c0
    print *, "omegar: ", omegar

    call this%bc%compute(t)
    call this%A%set_all(0.0_r8)
    call this%Ap%set_all(0.0_r8)
    this%rhs = 0

    ! print *, "setup1"

    do j = 1, this%mesh%ncell_onP
      ctm2c = matmul(transpose(curl), sym_matmul(W2_matrix_WE(this%meshr, j), curl))

      !! NB: we traverse the elements of the upper triangular part column by column.
      m1 = W1_matrix_WE(this%meshr, j)
      l = 0
      do e1x = 1, 6
        e1 = this%mesh%cedge(e1x,j)
        e1r = 2*(e1-1) + 1
        e1i = 2*(e1-1) + 2

        ! solve dummy equations on the Dirichlet edges
        if (this%bc%is_ebc_edge(e1)) then
          call this%A%set(e1r, e1r, 1.0_r8)
          call this%A%set(e1i, e1i, 1.0_r8)
          call this%Ap%set(e1, e1, 1.0_r8)
          this%rhs(e1r) = this%bc%efield(e1)
          this%rhs(e1i) = 0
          l = l + e1x
          cycle
        end if

        ! set up the system
        do e2x = 1, e1x
          l = l + 1
          e2 = this%mesh%cedge(e2x,j)
          e2r = 2*(e2-1) + 1
          e2i = 2*(e2-1) + 2

          ! mtr1 = ctm2c(e1x, e2x) / (omega * mu(j))
          ! mtr2 = omega * epsr(j) * m1(l)
          ! mtr3 = (omega * epsi(j) - sigma(j)) * m1(l)

          ! non-dimensionalized
          mtr1 = ctm2c(e1x, e2x) / mu(j)
          mtr2 = omegar**2  * epsr(j) * m1(l)
          mtr3 = (omegar**2 * epsi(j) - omegar * sigma(j) * this%Z0 * this%L0) * m1(l)
          !mtr1 = 0

          if (this%bc%is_ebc_edge(e2)) then
            ! Apply Dirichlet BCs to the RHS
            this%rhs(e1r) = this%rhs(e1r) - (mtr1 - mtr2) * this%bc%efield(e2)
            this%rhs(e1i) = this%rhs(e1i) - mtr3 * this%bc%efield(e2)
          else
            call this%A%add_to(e1r, e2r,  mtr1 - mtr2)
            call this%A%add_to(e1i, e2i, -mtr1 + mtr2)
            call this%A%add_to(e1r, e2i,  mtr3)
            call this%A%add_to(e1i, e2r,  mtr3)
            call this%Ap%add_to(e1, e2, mtr1 - mtr2 - mtr3) ! preconditioner

            if (e1x /= e2x) then
              call this%A%add_to(e2r, e1r,  mtr1 - mtr2)
              call this%A%add_to(e2i, e1i, -mtr1 + mtr2)
              call this%A%add_to(e2i, e1r,  mtr3)
              call this%A%add_to(e2r, e1i,  mtr3)
              call this%Ap%add_to(e2, e1, mtr1 - mtr2 - mtr3) ! preconditioner
            end if
          end if
        end do
      end do
    end do

    ! add the source to the rhs
    do e1 = 1, this%mesh%nedge
      e1r = 2*(e1-1) + 1
      e1i = 2*(e1-1) + 2
      !this%rhs(e1r) = this%rhs(e1r) + this%bc%hsource(e1)
      this%rhs(e1r) = this%rhs(e1r) + omegar * this%bc%hsource(e1) / this%H0

      !this%rhs(e1i) = this%rhs(e1i) + omegar * this%bc%hsource(e1) / this%H0
    end do

    call this%precon%setup(mu, epsr, epsi, sigma, omegar) ! this implicitly uses this%Ap
    this%model%rhs = this%rhs

    ASSERT(all(ieee_is_finite(this%model%rhs)))

    call stop_timer("setup")

  end subroutine setup


  subroutine solve(this)

    class(emfd_nlsol_solver), intent(inout) :: this

    integer :: ierr

    print '(a,2es13.3)', "max |rhs| = ", maxval(abs(this%rhs)), maxval(abs(this%bc%hsource))

    ASSERT(all(ieee_is_finite(this%efield)))
    ASSERT(all(ieee_is_finite(this%rhs)))

    call start_timer("solve")
    call this%solver%solve(0.0_r8, 0.0_r8, this%efield, this%efield, ierr)
    call stop_timer("solve")
    !call tls_info('  EMFD solve: ' // this%solver%metrics_string())
    print *, "ierr: ", ierr
    if (ierr /= 0) call tls_error("EMFD solve unsuccessful")

    print *, "max |rhs| = ", maxval(abs(this%rhs)), maxval(abs(this%bc%hsource))
    print *, "max |er| = ", maxval(abs(this%efield(::2)))
    print *, "max |ei| = ", maxval(abs(this%efield(2::2)))
    print *, "max |e-rhs| = ", maxval(abs(this%efield - this%rhs))

  end subroutine solve


  subroutine compute_heat_source(this, q)

    use mimetic_discretization, only: w1_matrix_we
    use parallel_communication

    class(emfd_nlsol_solver), intent(in) :: this
    real(r8), intent(out) :: q(:)

    real(r8) :: efieldt, efield2, q_joule, q_dielectric, m1(21)
    integer :: i, j, k, l1, l2, edges_r(6), edges_i(6)
    character(256) :: string

    ASSERT(size(q) == this%mesh%ncell)

    call start_timer("heat source")

    do j = 1, this%mesh%ncell
      if (this%sigma(j) == 0 .and. this%epsi(j) == 0) then
        q(j) = 0
        cycle
      end if

      ! compute |E|^2 on cell j
      m1 = W1_matrix_WE(this%mesh, j)
      edges_r = 2 * (this%mesh%cedge(:,j) - 1) + 1
      edges_i = 2 * (this%mesh%cedge(:,j) - 1) + 2
      associate(efield_r => this%efield(edges_r), efield_i => this%efield(edges_i))
        efield2 = 0

        ! real part
        l1 = 1
        do k = 1, 6
          efieldt = 0
          l2 = l1
          do i = 1, k-1
            efieldt = efieldt + m1(l2) * efield_r(i)
            l2 = l2 + 1
          end do
          do i = k, 6
            efieldt = efieldt + m1(l2) * efield_r(i)
            l2 = l2 + i
          end do
          efield2 = efield2 + efieldt * efield_r(k)
          l1 = l1 + k
        end do

        ! imaginary part
        l1 = 1
        do k = 1, 6
          efieldt = 0
          l2 = l1
          do i = 1, k-1
            efieldt = efieldt + m1(l2) * efield_i(i)
            l2 = l2 + 1
          end do
          do i = k, 6
            efieldt = efieldt + m1(l2) * efield_i(i)
            l2 = l2 + i
          end do
          efield2 = efield2 + efieldt * efield_i(k)
          l1 = l1 + k
        end do
      end associate

      ! compute heat sources
      efield2 = efield2 * (this%Z0 * this%H0)**2
      ! q_joule = efield2
      ! q_dielectric = 0
      q_joule = this%sigma(j) * efield2 / 2
      q_dielectric = this%omega * this%epsi(j) * efield2 / 2

      ! want a source density
      q(j) = (q_joule + q_dielectric) / abs(this%mesh%volume(j))
    end do

    write(string,fmt='(2(a,es11.4))') '|Q|_max=', global_maxval(q), ', Q_total=', &
        global_dot_product(q(:this%mesh%ncell_onP), abs(this%mesh%volume(:this%mesh%ncell_onP)))
    call tls_info(trim(string))

    call stop_timer("heat source")

  end subroutine compute_heat_source


!!!! MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine emfd_nlsol_model_init(this, A, precon, atol, rtol, ftol)
    class(emfd_nlsol_model), intent(out) :: this
    type(pcsr_matrix), pointer :: A
    type(PRECON), pointer :: precon
    real(r8), intent(in) :: atol, rtol, ftol
    this%A => A
    this%precon => precon
    this%atol = atol
    this%rtol = rtol
    this%ftol = ftol
  end subroutine emfd_nlsol_model_init


  integer function emfd_nlsol_model_size(this)
    class(emfd_nlsol_model), intent(in) :: this
    emfd_nlsol_model_size = this%A%nrow_onP
  end function emfd_nlsol_model_size


  subroutine emfd_nlsol_model_compute_f(this, t, u, udot, f, ax)
    class(emfd_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), udot(:)
    real(r8), intent(out), contiguous, target :: f(:)
    logical, intent(in), optional :: ax
    ASSERT(all(ieee_is_finite(u)))
    call this%A%matvec(u, f)
    if (present(ax)) then
      if (ax) return
    end if
    f = this%rhs - f
  end subroutine emfd_nlsol_model_compute_f


  subroutine emfd_nlsol_model_apply_precon(this, t, u, f)
    class(emfd_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:)
    real(r8), intent(inout), contiguous, target :: f(:)
    integer :: ierr
    this%b = f
    f = 0
    call this%precon%apply(this%b, f, ierr)
    INSIST(ierr == 0)
  end subroutine emfd_nlsol_model_apply_precon


  subroutine emfd_nlsol_model_compute_precon(this, t, u, dt)
    class(emfd_nlsol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), contiguous, target :: u(:)
    ! compute precon currently handled above the nlsol solver.
    !call this%precon%setup(interior_nodes)
  end subroutine emfd_nlsol_model_compute_precon


  real(r8) function emfd_nlsol_model_du_norm(this, t, u, du)

    use parallel_communication, only: global_maxval

    class(emfd_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), du(:)

    integer :: i
    real(r8) :: err

    emfd_nlsol_model_du_norm = 0
    do i = 1, size(u)
      if (this%atol == 0 .and. abs(u(i)) == 0) then
        err = huge(1.0_r8)
      else
        err = abs(du(i)) / (this%atol + this%rtol*abs(u(i)))
      end if
      emfd_nlsol_model_du_norm = max(emfd_nlsol_model_du_norm, err)
    end do
    emfd_nlsol_model_du_norm = global_maxval(emfd_nlsol_model_du_norm)

  end function emfd_nlsol_model_du_norm


  logical function emfd_nlsol_model_is_converged(this, itr, t, u, du, f_lnorm, tol)
    class(emfd_nlsol_model) :: this
    integer, intent(in) :: itr
    real(r8), intent(in) :: t, tol
    real(r8), intent(in), contiguous, target :: u(:), du(:), f_lnorm(:)
    emfd_nlsol_model_is_converged = this%du_norm(t, u, du) < tol .and. f_lnorm(3) < this%ftol
  end function emfd_nlsol_model_is_converged

end module emfd_nlsol_solver_type
