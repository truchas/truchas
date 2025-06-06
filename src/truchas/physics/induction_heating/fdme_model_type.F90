!!
!! FDME_MODEL_TYPE
!!
!! This module provides a derived type that encapsulates data and methods
!! for a discrete form of the frequency-domain Maxwell equations. The
!! discretization uses the mimetic Whitney finite element complex on a 3D
!! simplicial mesh.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fdme_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use complex_pcsr_matrix_type
  use pcsr_matrix_type
  use pbsr_matrix_type
  use bndry_func1_class
  use bndry_vfunc_class
  use bndry_cfunc1_class
  use truchas_timers
  use index_map_type
  use index_corrector_type
  implicit none
  private

  type, public :: fdme_model
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(complex_pcsr_matrix) :: A, B, BT
    complex(r8), allocatable :: rhs(:)
    type(pbsr_matrix) :: A2
    real(r8) :: omega
    real(r8), allocatable :: mu(:), sigma(:), epsi(:)
    complex(r8), allocatable :: k(:)
    class(bndry_func1), allocatable :: ebc  ! tangential E condition (nxE)
    class(bndry_func1), allocatable :: hbc  ! tangential H condition (nxH)
    class(bndry_cfunc1), allocatable :: robin_lhs, robin_rhs
    ! Non-dimensionalization parameters
    real(r8) :: Z0 ! vacuum impedance
    real(r8) :: c0 ! speed of light

    logical :: use_mixed_form
    type(pcsr_matrix) :: Am
    type(complex_pcsr_matrix) :: cAm
    type(index_corrector) :: icr, icc
    !FOR MINRES PC EXPERIMENTATION
    type(complex_pcsr_matrix) :: M
  contains
    procedure :: init
    procedure :: setup
    procedure :: residual
    procedure :: compute_bfield
    procedure :: compute_heat_source
    procedure :: compute_div
    procedure, private :: init_div_matrix
    procedure, private :: setup_div_matrix
    procedure, private :: init_mixed_matrix
    procedure, private :: setup_mixed_matrix
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
      call this%A2%init(2, this%A%graph, take_graph=.false.) ! used by AMS initialization
      !FOR MINRES PC EXPERIMENTATION
      call this%M%init(this%A%graph, take_graph=.false.)
    end block

    n = this%mesh%ncell
    allocate(this%epsi(n), this%sigma(n), this%mu(n), source=0.0_r8)
    allocate(this%k(n), source=cmplx(0,0,kind=r8))

    allocate(this%rhs(this%mesh%nedge))

    call params%get('use-mixed-form', this%use_mixed_form, stat, errmsg, default=.false.)
    if (stat /= 0) return
    call this%init_div_matrix ! currently used for diagnostics regardless of solver
    if (this%use_mixed_form) call this%init_mixed_matrix

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
    complex(r8), allocatable :: eps(:), k(:)

    ASSERT(size(epsr) == this%mesh%ncell)
    ASSERT(size(epsi) == this%mesh%ncell)
    ASSERT(size(mu) == this%mesh%ncell)
    ASSERT(size(sigma) == this%mesh%ncell)

    call start_timer("setup")

    k0 = omega / this%c0 ! free-space angular wave number
    this%omega = omega

    ! stored for preconditioner use
    this%mu(:) = mu
    this%k(:)%re = k0**2 * epsr
    this%k(:)%im = k0**2 * epsi + sigma*(k0*this%Z0)

    ! stored for heat calculation
    this%epsi(:) = epsi
    this%sigma(:) = sigma

    ! Raw system matrix ignoring boundary condtions
    do j = 1, this%mesh%ncell
      m1 = W1_matrix_WE(this%mesh, j)
      m2 = W2_matrix_WE(this%mesh, j)
      ctm2c = upm_cong_prod(4, 6, m2, cell_curl)
      Aj = (1.0_r8/mu(j)) * ctm2c - this%k(j) * m1
      call this%A%add_to(this%mesh%cedge(:,j), Aj)
      !FOR MINRES PC EXPERIMENTATION
      call this%M%add_to(this%mesh%cedge(:,j), this%k(j)*m1)
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
          !FOR MINRES PC EXPERIMENTATION
          call this%M%add_to(this%mesh%fedge(:,n), -a)
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
        !FOR MINRES PC EXPERIMENTATION
        call this%M%project_out(n)
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
    this%A2%values(1,1,:) =  this%A%values%re
    this%A2%values(2,1,:) =  this%A%values%im
    this%A2%values(1,2,:) = -this%A%values%im
    this%A2%values(2,2,:) =  this%A%values%re

    call this%setup_div_matrix(epsr, epsi) ! currently used for diagnostics regardless of system
    if (this%use_mixed_form) call this%setup_mixed_matrix

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
    r(:this%mesh%nedge_onP) = this%rhs(:this%mesh%nedge_onP) - r(:this%mesh%nedge_onP)
    !call this%mesh%edge_imap%gather_offp(r) ! not necessary?
  end subroutine

  subroutine compute_bfield(this, efield, bfield)
    use mimetic_discretization, only: curl
    class(fdme_model), intent(in) :: this
    complex(r8), intent(in)  :: efield(:)
    complex(r8), intent(out) :: bfield(:)
    ASSERT(size(efield) == this%mesh%nedge)
    ASSERT(size(bfield) == this%mesh%nface)
    bfield%re =  (1.0/this%omega) * curl(this%mesh,efield%im)
    bfield%im = (-1.0/this%omega) * curl(this%mesh,efield%re)
  end subroutine

  subroutine compute_heat_source(this, efield, q)

    use mimetic_discretization, only: w1_matrix_we
    use upper_packed_matrix_procs, only: upm_quad_form
    use physical_constants, only: vacuum_permittivity

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
      q_dielectric = this%omega*vacuum_permittivity*this%epsi(j) * efield2 / 2

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


  !!! MIXED FORMULATION ROUTINES !!!!!!!!!!!!!!!!!
  subroutine init_mixed_matrix(this)

    class(fdme_model), intent(inout) :: this

    integer :: i, xj, j, xnb, nb
    type(pcsr_graph), pointer :: gr => null(), gc => null()

    call start_timer("mixed")

    !! Set up a new imap for the 4-variable system: 2 edge-centered variables
    !! (real and imaginary E-field) and 2 node-centered variables (real and
    !! imaginary Lagrange multiplier).
    call this%icr%init([this%mesh%edge_imap, this%mesh%edge_imap, &
        this%mesh%node_imap, this%mesh%node_imap])
    call this%icc%init([this%mesh%edge_imap, this%mesh%node_imap])

    allocate(gr, gc)
    call gr%init(this%icr%imap)
    call gc%init(this%icc%imap)

    do i = 1, this%A%nrow
      do xj = this%A%graph%xadj(i), this%A%graph%xadj(i+1)-1
        j = this%A%graph%adjncy(xj)

        call gr%add_edge(this%icr%eval(i,1), this%icr%eval(j,1))
        call gr%add_edge(this%icr%eval(i,1), this%icr%eval(j,2))
        call gr%add_edge(this%icr%eval(i,2), this%icr%eval(j,1))
        call gr%add_edge(this%icr%eval(i,2), this%icr%eval(j,2))

        call gc%add_edge(this%icc%eval(i,1), this%icc%eval(j,1))
      end do
    end do

    do i = 1, this%B%nrow ! nnode
      do xj = this%B%graph%xadj(i), this%B%graph%xadj(i+1)-1
        j = this%B%graph%adjncy(xj) ! nedge

        call gr%add_edge(this%icr%eval(i,3), this%icr%eval(j,1))
        call gr%add_edge(this%icr%eval(i,4), this%icr%eval(j,2))
        call gr%add_edge(this%icr%eval(j,1), this%icr%eval(i,3))
        call gr%add_edge(this%icr%eval(j,2), this%icr%eval(i,4))

        call gr%add_edge(this%icr%eval(i,3), this%icr%eval(j,2))
        call gr%add_edge(this%icr%eval(i,4), this%icr%eval(j,1))
        call gr%add_edge(this%icr%eval(j,1), this%icr%eval(i,4))
        call gr%add_edge(this%icr%eval(j,2), this%icr%eval(i,3))

        call gc%add_edge(this%icc%eval(i,2), this%icc%eval(j,1))
        call gc%add_edge(this%icc%eval(j,1), this%icc%eval(i,2))
      end do
      call gc%add_edge(this%icc%eval(i,2), this%icc%eval(i,2)) ! used only for kdiag in minres
    end do

    ! for BCs
    do xj = 1, size(this%ebc%index)
      j = this%ebc%index(xj)
      do xnb = 1, 2
        nb = this%mesh%enode(xnb, j)
        call gr%add_edge(this%icr%eval(nb,3), this%icr%eval(nb,3))
        call gr%add_edge(this%icr%eval(nb,4), this%icr%eval(nb,4))
        call gc%add_edge(this%icc%eval(nb,2), this%icc%eval(nb,2))
      end do
    end do

    call gr%add_complete
    call gc%add_complete
    call this%Am%init(gr, take_graph=.true.)
    call this%cAm%init(gc, take_graph=.true.)

    call stop_timer("mixed")

  end subroutine init_mixed_matrix


  subroutine setup_mixed_matrix(this)

    class(fdme_model), intent(inout) :: this

    integer :: i, xj, j, xnb, nb

    call this%Am%set_all(0.0_r8)
    call this%cAm%set_all(cmplx(0, 0, kind=r8))

    do i = 1, this%A%nrow ! nedge
      do xj = this%A%graph%xadj(i), this%A%graph%xadj(i+1)-1
        j = this%A%graph%adjncy(xj) ! nedge
        call this%Am%set(this%icr%eval(i,1), this%icr%eval(j,1),  this%A%values(xj)%re)
        call this%Am%set(this%icr%eval(i,1), this%icr%eval(j,2), -this%A%values(xj)%im)
        call this%Am%set(this%icr%eval(i,2), this%icr%eval(j,1),  this%A%values(xj)%im)
        call this%Am%set(this%icr%eval(i,2), this%icr%eval(j,2),  this%A%values(xj)%re)
        call this%cAm%set(this%icc%eval(i,1), this%icc%eval(j,1),  this%A%values(xj))
      end do
    end do

    do i = 1, this%B%nrow ! nnode
      do xj = this%B%graph%xadj(i), this%B%graph%xadj(i+1)-1
        j = this%B%graph%adjncy(xj) ! nedge

        call this%Am%set(this%icr%eval(i,3), this%icr%eval(j,1),  this%B%values(xj)%re)
        call this%Am%set(this%icr%eval(i,4), this%icr%eval(j,2), -this%B%values(xj)%re)
        call this%Am%set(this%icr%eval(j,1), this%icr%eval(i,3),  this%B%values(xj)%re) ! B^T
        call this%Am%set(this%icr%eval(j,2), this%icr%eval(i,4), -this%B%values(xj)%re) ! B^T

        call this%Am%set(this%icr%eval(i,3), this%icr%eval(j,2), -this%B%values(xj)%im)
        call this%Am%set(this%icr%eval(i,4), this%icr%eval(j,1), -this%B%values(xj)%im)
        call this%Am%set(this%icr%eval(j,1), this%icr%eval(i,4), -this%B%values(xj)%im) ! B^T
        call this%Am%set(this%icr%eval(j,2), this%icr%eval(i,3), -this%B%values(xj)%im) ! B^T

        call this%cAm%set(this%icc%eval(i,2), this%icc%eval(j,1),  this%B%values(xj))
        call this%cAm%set(this%icc%eval(j,1), this%icc%eval(i,2),  this%B%values(xj)) ! B^T
      end do
    end do

    ! Apply nxE boundary conditions to the system matrix
    if (allocated(this%ebc)) then
      do xj = 1, size(this%ebc%index)
        j = this%ebc%index(xj)
        do xnb = 1, 2
          nb = this%mesh%enode(xnb, j)
          call this%Am%set(this%icr%eval(nb,3), this%icr%eval(nb,3), 1.0_r8)
          call this%Am%set(this%icr%eval(nb,4), this%icr%eval(nb,4), 1.0_r8)
          call this%cAm%set(this%icc%eval(nb,2), this%icc%eval(nb,2), cmplx(1,1,kind=r8))
        end do
      end do
    end if

  end subroutine setup_mixed_matrix


end module fdme_model_type
