#include "f90_assert.fpp"

module fdme_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  !use pcsr_matrix_type
  use msr_matrix_type
  use bndry_func1_class
  use truchas_timers
  implicit none
  private

  type, public :: fdme_model
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(msr_matrix) :: A(2,2)
    real(r8), allocatable :: rhs(:)
    real(r8) :: omega
    real(r8), allocatable :: epsi(:), epsr(:), mu(:), sigma(:)
    class(bndry_func1), allocatable :: ebc  ! tangential E condition (nxE)
    class(bndry_func1), allocatable :: hbc  ! tangential H condition (nxH)
    ! Non-dimensionalization parameters
    real(r8) :: Z0 = 376.730313668 ! Ohms -- vacuum impedance
    real(r8) :: c0 = 299792458.0_r8 ! speed of light
  contains
    procedure :: init
    procedure :: setup
    procedure :: compute_f
    procedure :: compute_heat_source
  end type

contains

  subroutine init(this, mesh, bc_fac, params, stat, errmsg)

    use em_bc_factory_type
    use parameter_list_type

    class(fdme_model), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n

    this%mesh => mesh

    !! Boundary condition data
    call bc_fac%alloc_nxE_bc(this%ebc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(this%ebc)) then
      call bc_fac%alloc_fd_nxH_bc(this%hbc, stat, errmsg, omit_edge_list=this%ebc%index)
    else
      call bc_fac%alloc_fd_nxH_bc(this%hbc, stat, errmsg)
    end if
    if (stat /= 0) return

    block ! full system in block form
      !type(pcsr_graph), pointer :: g
      type(msr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nedge)
      call g%add_clique(this%mesh%cedge)
      call g%add_complete
      call this%A(1,1)%init(g, take_graph=.true.)
      call this%A(1,2)%init(mold=this%A(1,1))
      call this%A(2,1)%init(mold=this%A(1,1))
      call this%A(2,2)%init(mold=this%A(1,1))
    end block

    n = this%mesh%ncell
    allocate(this%epsr(n), this%epsi(n), this%mu(n), this%sigma(n), source=0.0_r8)

    allocate(this%rhs(2*this%mesh%nedge), source=0.0_r8)

  end subroutine init

  subroutine setup(this, t, epsr, epsi, mu, sigma, omega)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we, cell_curl
    use upper_packed_matrix_procs, only: upm_cong_prod

    class(fdme_model), intent(inout) :: this
    real(r8), intent(in) :: t, epsr(:), epsi(:), mu(:), sigma(:), omega

    integer :: j, n
    real(r8) :: m1(21), m2(10), ctm2c(21), a(21), omegar
    real(r8) :: rhs_r(this%mesh%nedge), rhs_i(this%mesh%nedge)

    ASSERT(size(epsr) == this%mesh%ncell)
    ASSERT(size(epsi) == this%mesh%ncell)
    ASSERT(size(mu) == this%mesh%ncell)
    ASSERT(size(sigma) == this%mesh%ncell)

    call start_timer("setup")

    this%epsr(:) = epsr
    this%epsi(:) = epsi
    this%mu(:) = mu
    this%sigma(:) = sigma

    ! non-dimensionalization
    this%omega = omega
    omegar = omega / this%c0

    call this%A(1,1)%set_all(0.0_r8)
    call this%A(1,2)%set_all(0.0_r8)
    call this%A(2,1)%set_all(0.0_r8)
    call this%A(2,2)%set_all(0.0_r8)

    ! Raw system matrix ignoring boundary condtions
    do j = 1, this%mesh%ncell_onP !TODO: this should run over ALL cells
      m1 = W1_matrix_WE(this%mesh, j)
      m2 = W2_matrix_WE(this%mesh, j)
      ctm2c = upm_cong_prod(4, 6, m2, cell_curl)

      a = (1.0_r8/mu(j)) * ctm2c - (omegar**2 * epsr(j)) * m1
      call this%A(1,1)%add_to(this%mesh%cedge(:,j), a)
      call this%A(2,2)%add_to(this%mesh%cedge(:,j), -a)

      a = (omegar**2 * epsi(j) - omegar * sigma(j) * this%Z0) * m1
      call this%A(1,2)%add_to(this%mesh%cedge(:,j), a)
      call this%A(2,1)%add_to(this%mesh%cedge(:,j), a)
    end do

    ! RHS contribution from nxE boundary conditions
    if (allocated(this%ebc)) then
      block
        real(r8), allocatable :: efield_r(:)
        call this%ebc%compute(t)
        allocate(efield_r(this%mesh%nedge), source=0.0_r8)
        do j = 1, size(this%ebc%index)
          efield_r(this%ebc%index(j)) = this%ebc%value(j)
        end do
        rhs_r = efield_r - this%A(1,1)%matvec(efield_r)
        rhs_i = - this%A(2,1)%matvec(efield_r)
      end block
    end if

    ! RHS contribution from nxH source boundary conditions
    if (allocated(this%hbc)) then
      call this%hbc%compute(t)
      do j = 1, size(this%hbc%index)
        n = this%hbc%index(j)
        rhs_r(n) = rhs_r(n) + omegar * this%hbc%value(j)
      end do
    end if

    this%rhs(1:2*this%mesh%nedge:2) = rhs_r
    this%rhs(2:2*this%mesh%nedge:2) = rhs_i

    ! Apply nxE boundary conditions to the system matrix
    if (allocated(this%ebc)) then
      do j = 1, size(this%ebc%index)
        n = this%ebc%index(j)
        call this%A(1,1)%project_out(n)
        call this%A(1,2)%project_out(n)
        call this%A(2,1)%project_out(n)
        call this%A(2,2)%project_out(n)
        call this%A(1,1)%set(n, n, 1.0_r8)
        call this%A(2,2)%set(n, n, 1.0_r8)
      end do
    end if

    call stop_timer("setup")

  end subroutine setup


  subroutine compute_f(this, u, f, ax)
    class(fdme_model) :: this
    real(r8), intent(in) :: u(:)
    real(r8), intent(out) :: f(:)
    logical, intent(in), optional :: ax
    f(1::2) = this%A(1,1)%matvec(u(1::2)) + this%A(1,2)%matvec(u(2::2))
    f(2::2) = this%A(2,1)%matvec(u(1::2)) + this%A(2,2)%matvec(u(2::2))
    !f(2::2) = -f(2::2)
    if (present(ax)) then
      if (ax) return
    end if
    f = this%rhs - f
  end subroutine

  subroutine compute_heat_source(this, e, q)

    use mimetic_discretization, only: w1_matrix_we
    use upper_packed_matrix_procs, only: upm_quad_form
    use parallel_communication
    use truchas_logging_services

    class(fdme_model), intent(in) :: this
    real(r8), intent(in), target :: e(:)
    real(r8), intent(out) :: q(:)

    real(r8) :: efield2, q_joule, q_dielectric, m1(21)
    integer :: j
    character(256) :: string
    real(r8), pointer :: efield(:,:)

    ASSERT(size(q) == this%mesh%ncell)

    efield(1:2,1:this%mesh%nedge) => e

    do j = 1, this%mesh%ncell_onP
      if (this%sigma(j) == 0 .and. this%epsi(j) == 0) then
        q(j) = 0
        cycle
      end if

      ! compute |E|^2 on cell j
      associate(efield_r => efield(1,this%mesh%cedge(:,j)), efield_i => efield(2,this%mesh%cedge(:,j)))
        m1 = W1_matrix_WE(this%mesh, j)
        efield2 = this%Z0**2 * (upm_quad_form(m1, efield_r) + upm_quad_form(m1, efield_i))
      end associate

      ! compute heat sources
      q_joule = this%sigma(j) * efield2 / 2
      q_dielectric = (this%omega/this%c0) * this%epsi(j) * efield2 / 2

      ! want a source density
      q(j) = (q_joule + q_dielectric) / abs(this%mesh%volume(j))
    end do
    call this%mesh%cell_imap%gather_offp(q)

    write(string,fmt='(2(a,es11.4))') '|Q|_max=', global_maxval(q), ', Q_total=', &
        global_dot_product(q(:this%mesh%ncell_onP), abs(this%mesh%volume(:this%mesh%ncell_onP)))
    call tls_info(trim(string))

  end subroutine compute_heat_source


end module fdme_model_type
