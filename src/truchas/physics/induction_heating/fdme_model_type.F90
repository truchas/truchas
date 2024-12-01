#include "f90_assert.fpp"

module fdme_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use complex_pcsr_matrix_type
  use pbsr_matrix_type
  use bndry_func1_class
  use bndry_vfunc_class
  use bndry_cfunc1_class
  use truchas_timers
  implicit none
  private

  type, public :: fdme_model
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(complex_pcsr_matrix) :: A, B, BT
    complex(r8), allocatable :: rhs(:)
    type(pbsr_matrix) :: A2
    real(r8) :: omega
    real(r8), allocatable :: epsi(:), epsr(:), mu(:), sigma(:)
    class(bndry_func1), allocatable :: ebc  ! tangential E condition (nxE)
    class(bndry_func1), allocatable :: hbc  ! tangential H condition (nxH)
    class(bndry_cfunc1), allocatable :: robin_lhs, robin_rhs
    ! Non-dimensionalization parameters
    real(r8) :: Z0 ! vacuum impedance
    real(r8) :: c0 ! speed of light
    logical :: use_mixed_form = .false.
  contains
    procedure :: init
    procedure :: setup
    procedure :: residual
    procedure :: compute_bfield
    procedure :: compute_heat_source
    procedure :: compute_div
    procedure, private :: init_div_matrix
    procedure, private :: setup_div_matrix
  end type

contains

  subroutine init(this, mesh, bc_fac, params, stat, errmsg)

    use em_bc_factory_type
    use parameter_list_type
    use physical_constants, only: vacuum_permittivity, vacuum_permeability

    class(fdme_model), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(em_bc_factory), intent(in) :: bc_fac
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n

    this%Z0 = sqrt(vacuum_permeability/vacuum_permittivity)
    this%c0 = 1.0_r8/sqrt(vacuum_permittivity*vacuum_permeability)
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

    if (allocated(this%ebc)) then
      call bc_fac%alloc_robin_bc(this%robin_lhs, this%robin_rhs, stat, errmsg, omit_edge_list=this%ebc%index)
    else
      call bc_fac%alloc_robin_bc(this%robin_lhs, this%robin_rhs, stat, errmsg)
    end if
    if (stat /= 0) return

    block
      type(pcsr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%edge_imap)
      call g%add_clique(this%mesh%cedge)
      call g%add_complete
      call this%A%init(g, take_graph=.true.)
    end block

    n = this%mesh%ncell
    allocate(this%epsr(n), this%epsi(n), this%mu(n), this%sigma(n), source=0.0_r8)

    allocate(this%rhs(this%mesh%nedge))

    call this%init_div_matrix ! currently used for diagnostics regardless of solver

  end subroutine init

  subroutine init_div_matrix(this)

    class(fdme_model), intent(inout) :: this

    integer :: j, xn, n, xe, e
    type(pcsr_graph), pointer :: g1, g2

    call start_timer("div")

    allocate(g1, g2)
    call g1%init(this%mesh%node_imap, this%mesh%edge_imap)
    call g2%init(this%mesh%edge_imap, this%mesh%node_imap)
    do j = 1, this%mesh%ncell
      do xn = 1, size(this%mesh%cnode(:,j))
        n = this%mesh%cnode(xn,j)
        do xe = 1, size(this%mesh%cedge(:,j))
          e = this%mesh%cedge(xe,j)
          call g1%add_edge(n, e)
          call g2%add_edge(e, n)
        end do
      end do
    end do
    call g1%add_complete
    call g2%add_complete
    call this%B%init(g1, take_graph=.true.)
    call this%BT%init(g2, take_graph=.true.)

    call stop_timer("div")

  end subroutine init_div_matrix

  subroutine setup(this, t, epsr, epsi, mu, sigma, omega)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we, cell_curl, w1_face_matrix
    use upper_packed_matrix_procs, only: upm_cong_prod

    class(fdme_model), intent(inout) :: this
    real(r8), intent(in) :: t, epsr(:), epsi(:), mu(:), sigma(:), omega

    integer :: j, n
    real(r8) :: m1(21), m2(10), ctm2c(21), k0
    complex(r8) :: Aj(21)

    ASSERT(size(epsr) == this%mesh%ncell)
    ASSERT(size(epsi) == this%mesh%ncell)
    ASSERT(size(mu) == this%mesh%ncell)
    ASSERT(size(sigma) == this%mesh%ncell)

    call start_timer("setup")

    k0 = omega / this%c0 ! free-space angular wave number
    this%omega = omega

    this%mu(:) = mu
    this%epsr(:) = epsr
    this%epsi(:) = epsi + sigma*(this%Z0/k0) ! fold cond into complex perm
    this%sigma(:) = 0.0_r8 ! it's been folded into the complex permittivity

    ! Raw system matrix ignoring boundary condtions
    do j = 1, this%mesh%ncell
      m1 = W1_matrix_WE(this%mesh, j)
      m2 = W2_matrix_WE(this%mesh, j)
      ctm2c = upm_cong_prod(4, 6, m2, cell_curl)
      Aj%re = (1.0_r8/mu(j)) * ctm2c - (epsr(j)*k0**2) * m1
      Aj%im = (-this%epsi(j)*k0**2) * m1
      call this%A%add_to(this%mesh%cedge(:,j), Aj)
    end do

    ! LHS contribution from Robin boundary conditions
    !FIXME: ONLY CORRECT FOR MU=1
    if (allocated(this%robin_lhs)) then
      block
        complex(r8) ::a(6)
        call this%robin_lhs%compute(t)
        do j = 1, size(this%robin_lhs%index)
          n = this%robin_lhs%index(j)
          a = -this%robin_lhs%value(j) * w1_face_matrix(this%mesh, n)
          call this%A%add_to(this%mesh%fedge(:,n), a)
        end do
      end block
    end if

    this%rhs = 0.0_r8

    ! RHS contribution from nxE boundary conditions
    if (allocated(this%ebc)) then
      block
        complex(r8) :: efield(this%mesh%nedge)
        efield = 0.0_r8
        do j = 1, size(this%ebc%index)
          efield(this%ebc%index(j))%re = this%ebc%value(j)
        end do
        call this%A%matvec(efield, this%rhs)
        this%rhs = efield - this%rhs
        call this%mesh%edge_imap%gather_offp(this%rhs)
      end block
    end if

    ! RHS contribution from Robin boundary conditions
    !FIXME: only correct for uniform mu = 1 (relative). For other mu, it needs to
    !be incorporated into the computation of robin_rhs.
    if (allocated(this%robin_rhs)) then
      call this%robin_rhs%compute(t)
      do j = 1, size(this%robin_rhs%index)
        n = this%robin_rhs%index(j)
        this%rhs(n) = this%rhs(n) - this%robin_rhs%value(j)
      end do
      call this%mesh%edge_imap%gather_offp(this%rhs) ! necessary?
    end if

    !! Apply the nxE boundary conditions to the system matrix
    if (allocated(this%ebc)) then
      do j = 1, size(this%ebc%index)
        n = this%ebc%index(j)
        call this%A%project_out(n)
        call this%A%set(n, n, cmplx(1,0,kind=r8))
      end do
    end if

    ! RHS contribution from nxH source boundary conditions
    if (allocated(this%hbc)) then
      call this%hbc%compute(t)
      do j = 1, size(this%hbc%index)
        n = this%hbc%index(j)
        this%rhs(n)%im = this%rhs(n)%im - k0 * this%Z0 * this%hbc%value(j)
      end do
      call this%mesh%edge_imap%gather_offp(this%rhs)
    end if

    !! Copy the complex CSR matrix to an equivalent real 2x2 block CSR matrix
    call this%A2%init(2, this%A%graph, take_graph=.false.)
    this%A2%values(1,1,:) =  this%A%values%re
    this%A2%values(2,1,:) =  this%A%values%im
    this%A2%values(1,2,:) = -this%A%values%im
    this%A2%values(2,2,:) =  this%A%values%re

    call this%setup_div_matrix(epsr, epsi) ! currently used for diagnostics regardless of system

    call stop_timer("setup")

  end subroutine setup


  subroutine setup_div_matrix(this, epsr, epsi)

    use mimetic_discretization, only: w1_matrix_we, cell_grad
    use upper_packed_matrix_procs, only: sym_matmul

    class(fdme_model), intent(inout) :: this
    real(r8), intent(in) :: epsr(:), epsi(:)

    real(r8) :: bT(6, 4), b(4, 6), m1(21)
    integer :: j, xn, n, xe, e
    complex(r8) :: c

    call this%B%set_all(cmplx(0, 0, kind=r8))

    do j = 1, this%mesh%ncell
      m1 = W1_matrix_WE(this%mesh, j)
      bT = sym_matmul(m1, cell_grad)
      b = transpose(bT)
      c = cmplx(epsr(j), epsi(j), kind=r8)
      do xn = 1, size(this%mesh%cnode(:,j))
        n = this%mesh%cnode(xn,j)
        do xe = 1, size(this%mesh%cedge(:,j))
          e = this%mesh%cedge(xe,j)
          call this%B%add_to(n, e, c*b(xn, xe))
          call this%BT%add_to(e, n, c*bT(xe, xn))
        end do
      end do
    end do

    ! zero Lagrange multipliers on nxE BCs
    block
      integer :: xj, nb, xnb, Nn
      if (allocated(this%ebc)) then
        do xj = 1, size(this%ebc%index)
          j = this%ebc%index(xj)
          where (j == this%B%graph%adjncy) this%B%values = 0
          do xnb = 1, 2
            nb = this%mesh%enode(xnb, j)
            this%B%values(this%B%graph%xadj(nb):this%B%graph%xadj(nb+1)-1) = 0
          end do
        end do
      end if
    end block

  end subroutine setup_div_matrix


  subroutine residual(this, e, r)
    class(fdme_model), intent(in) :: this
    complex(r8), intent(inout) :: e(:), r(:)
    call this%mesh%edge_imap%gather_offp(e)
    call this%A%matvec(e, r)
    r = this%rhs - r
    !call this%mesh%edge_imap%gather_offp(r) ! not necessary?
  end subroutine

  subroutine compute_bfield(this, efield, bfield)
    use mimetic_discretization, only: curl
    class(fdme_model), intent(in) :: this
    complex(r8), intent(in)  :: efield(:)
    complex(r8), intent(out) :: bfield(:)
    ASSERT(size(efield) == this%mesh%nedge)
    ASSERT(size(bfield) == this%mesh%nface)
    bfield%re = (-1.0/this%omega) * curl(this%mesh,efield%im)
    bfield%im =  (1.0/this%omega) * curl(this%mesh,efield%re)
  end subroutine

  subroutine compute_heat_source(this, efield, q)

    use mimetic_discretization, only: w1_matrix_we
    use upper_packed_matrix_procs, only: upm_quad_form

    class(fdme_model), intent(in) :: this
    complex(r8), intent(in) :: efield(:)
    real(r8), intent(out) :: q(:)

    real(r8) :: efield2, q_joule, q_dielectric
    integer :: j

    ASSERT(size(efield) == this%mesh%nedge)
    ASSERT(size(q) <= this%mesh%ncell)

    do j = 1, size(q)
      if (this%sigma(j) == 0 .and. this%epsi(j) == 0) then
        q(j) = 0
        cycle
      end if

      ! compute |E|^2 on cell j
      efield2 = upm_quad_form(W1_matrix_WE(this%mesh,j), efield(this%mesh%cedge(:,j)))

      ! compute heat sources
      q_joule = this%sigma(j) * efield2 / 2
      q_dielectric = (this%omega/this%c0) * this%epsi(j) * efield2 / 2

      ! want a source density
      q(j) = (q_joule + q_dielectric) / abs(this%mesh%volume(j))
    end do

  end subroutine compute_heat_source

  subroutine compute_div(this, efield, div_efield)
    class(fdme_model), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:) ! inout for gather_offp
    complex(r8), intent(out) :: div_efield(:)
    ASSERT(size(efield) == this%mesh%nedge)
    ASSERT(size(div_efield) == this%mesh%nnode)
    call this%mesh%edge_imap%gather_offp(efield)
    call this%B%matvec(efield, div_efield)
    call this%mesh%node_imap%gather_offp(div_efield)
  end subroutine compute_div

end module fdme_model_type
