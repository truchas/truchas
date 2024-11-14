#include "f90_assert.fpp"

module fdme_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use vector_class
  use fdme_vector_type
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
    type(complex_pcsr_matrix) :: AA
    type(pbsr_matrix) :: A
    type(fdme_vector) :: rhs
    complex(r8), allocatable :: crhs(:)
    real(r8) :: omega
    real(r8), allocatable :: epsi(:), epsr(:), mu(:), sigma(:)
    class(bndry_func1), allocatable :: ebc  ! tangential E condition (nxE)
    class(bndry_func1), allocatable :: hbc  ! tangential H condition (nxH)
    class(bndry_cfunc1), allocatable :: robin_lhs, robin_rhs
    ! Non-dimensionalization parameters
    real(r8) :: Z0 ! vacuum impedance
    real(r8) :: c0 ! speed of light
  contains
    procedure :: init
    procedure :: setup
    procedure :: matvec
    procedure :: matvec2
    procedure :: residual
    procedure :: compute_heat_source
    procedure :: compute_b
  end type

contains

  subroutine compute_b(this, e, b)
    use mimetic_discretization, only: curl
    class(fdme_model), intent(in) :: this
    real(r8), intent(in)  :: e(:,:)
    real(r8), intent(out) :: b(:,:)
    ASSERT(size(e,1) == 2)
    ASSERT(size(e,2) == this%mesh%nedge)
    ASSERT(size(b,1) == 2)
    ASSERT(size(b,2) == this%mesh%nface)
    b(1,:) = (-1.0/this%omega) * curl(this%mesh,e(2,:))
    b(2,:) =  (1.0/this%omega) * curl(this%mesh,e(1,:))
    call this%mesh%face_imap%gather_offp(b)
  end subroutine

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

    call bc_fac%alloc_robin_bc(this%robin_lhs, this%robin_rhs, stat, errmsg)
    if (stat /= 0) return

    block
      type(pcsr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%edge_imap)
      call g%add_clique(this%mesh%cedge)
      call g%add_complete
      call this%AA%init(g, take_graph=.true.)
    end block

    n = this%mesh%ncell
    allocate(this%epsr(n), this%epsi(n), this%mu(n), this%sigma(n), source=0.0_r8)

    allocate(this%crhs(this%mesh%nedge))
    call this%rhs%init(this%mesh)

  end subroutine init

  subroutine setup(this, t, epsr, epsi, mu, sigma, omega)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we, cell_curl, w1_face_matrix
    use upper_packed_matrix_procs, only: upm_cong_prod

    class(fdme_model), intent(inout) :: this
    real(r8), intent(in) :: t, epsr(:), epsi(:), mu(:), sigma(:), omega

    integer :: j, n
    real(r8) :: m1(21), m2(10), ctm2c(21), k0
    complex(r8) :: AA(21)

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
      AA%re = (1.0_r8/mu(j)) * ctm2c - (epsr(j)*k0**2) * m1
      AA%im = (-this%epsi(j)*k0**2) * m1
      call this%AA%add_to(this%mesh%cedge(:,j), AA)
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
          call this%AA%add_to(this%mesh%fedge(:,n), a)
        end do
      end block
    end if

    this%crhs = 0.0_r8

    ! RHS contribution from nxE boundary conditions
    if (allocated(this%ebc)) then
      block
        complex(r8) :: efield(this%mesh%nedge)
        efield = 0.0_r8
        do j = 1, size(this%ebc%index)
          efield(this%ebc%index(j))%re = this%ebc%value(j)
        end do
        call this%AA%matvec(efield, this%crhs)
        this%crhs = efield - this%crhs
        call this%mesh%edge_imap%gather_offp(this%crhs)
      end block
    end if

    ! RHS contribution from Robin boundary conditions
    !FIXME: only correct for uniform mu = 1 (relative). For other mu, it needs to
    !be incorporated into the computation of robin_rhs.
    if (allocated(this%robin_rhs)) then
      call this%robin_rhs%compute(t)
      do j = 1, size(this%robin_rhs%index)
        n = this%robin_rhs%index(j)
        this%crhs(n) = this%crhs(n) - this%robin_rhs%value(j)
      end do
      call this%mesh%edge_imap%gather_offp(this%crhs) ! necessary?
    end if

    !! Apply the nxE boundary conditions to the system matrix
    if (allocated(this%ebc)) then
      do j = 1, size(this%ebc%index)
        n = this%ebc%index(j)
        call this%AA%project_out(n)
        call this%AA%set(n, n, cmplx(1,0,kind=r8))
      end do
    end if

    ! RHS contribution from nxH source boundary conditions
    if (allocated(this%hbc)) then
      call this%hbc%compute(t)
      do j = 1, size(this%hbc%index)
        n = this%hbc%index(j)
        this%crhs(n)%im = this%crhs(n)%im - k0 * this%Z0 * this%hbc%value(j)
      end do
      call this%mesh%edge_imap%gather_offp(this%crhs)
    end if

    !! Copy the complex CSR matrix to an equivalent real 2x2 block CSR matrix
    call this%A%init(2, this%AA%graph, take_graph=.false.)
    this%A%values(1,1,:) =  this%AA%values%re
    this%A%values(2,1,:) =  this%AA%values%im
    this%A%values(1,2,:) = -this%AA%values%im
    this%A%values(2,2,:) =  this%AA%values%re

    !! Copy the complex RHS vector to an equivalent real RHS vector
    this%rhs%array(1,:) = this%crhs%re
    this%rhs%array(2,:) = this%crhs%im

    call stop_timer("setup")

  end subroutine setup


  subroutine residual(this, e, r)
    class(fdme_model), intent(in) :: this
    class(fdme_vector), intent(inout) :: e, r
    call e%gather_offp
    call this%A%matvec(e%array, r%array)
    r%array = this%rhs%array - r%array
    !call r%gather_offp ! not necessary
  end subroutine

  subroutine matvec(this, x, ax)
    class(fdme_model), intent(in) :: this
    type(fdme_vector), intent(inout) :: x, ax
    call this%A%matvec(x%array, ax%array)
  end subroutine

  subroutine matvec2(this, x, ax)
    class(fdme_model), intent(in) :: this
    complex(r8), intent(in) :: x(:)
    complex(r8), intent(out) :: ax(:)
    call this%AA%matvec(x, ax)
  end subroutine

  subroutine compute_heat_source(this, efield, q)

    use mimetic_discretization, only: w1_matrix_we
    use upper_packed_matrix_procs, only: upm_quad_form

    class(fdme_model), intent(in) :: this
    type(fdme_vector), intent(in) :: efield
    real(r8), intent(out) :: q(:)

    real(r8) :: efield2, q_joule, q_dielectric, m1(21)
    integer :: j

    ASSERT(size(q) <= this%mesh%ncell)

    do j = 1, size(q)
      if (this%sigma(j) == 0 .and. this%epsi(j) == 0) then
        q(j) = 0
        cycle
      end if

      ! compute |E|^2 on cell j
      associate(efield_r => efield%array(1,this%mesh%cedge(:,j)), &
                efield_i => efield%array(2,this%mesh%cedge(:,j)))
        m1 = W1_matrix_WE(this%mesh, j)
        efield2 = upm_quad_form(m1, efield_r) + upm_quad_form(m1, efield_i)
      end associate

      ! compute heat sources
      q_joule = this%sigma(j) * efield2 / 2
      q_dielectric = (this%omega/this%c0) * this%epsi(j) * efield2 / 2

      ! want a source density
      q(j) = (q_joule + q_dielectric) / abs(this%mesh%volume(j))
    end do

  end subroutine compute_heat_source

end module fdme_model_type
