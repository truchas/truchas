!!
!! EM_FD_solver_type
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

module em_fd_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use hypre_hybrid_type
  use pcsr_matrix_type
  use simpl_mesh_type
  use index_map_type
  use em_bc_type
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: em_fd_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(simpl_mesh) :: meshr

    type(hypre_hybrid) :: solver
    type(index_map), pointer :: imap => null()
    type(pcsr_matrix), pointer :: A => null()

    type(em_bc), pointer :: bc => null()

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
    final :: em_fd_solver_delete
  end type em_fd_solver

contains

  !! Final subroutine for em_fd_solver type objects.
  subroutine em_fd_solver_delete(this)
    type(em_fd_solver), intent(inout) :: this
    if (associated(this%imap)) deallocate(this%imap)
    if (associated(this%A)) deallocate(this%A)
    if (associated(this%bc)) deallocate(this%bc)
  end subroutine em_fd_solver_delete


  subroutine init(this, mesh, bc, params)

    use mimetic_discretization
    use parameter_list_type
    use parallel_communication, only: is_IOP

    class(em_fd_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(em_bc), intent(inout), pointer :: bc
    type(parameter_list), intent(inout), pointer :: params

    type(pcsr_graph), pointer :: g
    type(index_map), pointer :: row_imap
    integer, allocatable :: nvars(:)
    integer, allocatable :: ebedge(:)
    integer :: e, c, er, ei
    integer :: e1, e1x, e1r, e1i, e2, e2x, e2r, e2i

    !ASSERT(size(emask) == mesh%nedge)
    ASSERT(associated(bc))

    this%mesh => mesh
    ! this%epsr = epsr
    ! this%epsi = epsi
    ! this%mu = mu
    ! this%sigma = sigma
    ! this%omega = omega
    !this%emask = emask
    this%bc => bc ! taking ownership
    allocate(this%efield(2*mesh%nedge), this%rhs(2*mesh%nedge))
    this%efield = 0
    this%rhs = 0

    ! non-dimensionalization
    this%L0 = 1
    this%H0 = 1 !maxval(abs(this%bc%hsource))
    this%meshr = this%mesh
    this%meshr%x = this%meshr%x / this%L0
    this%meshr%length = this%meshr%length / this%L0
    this%meshr%volume = this%meshr%volume / this%L0**3

    !allocate(this%mtr1(21,mesh%ncell), this%mtr2(21,mesh%ncell), this%mtr3(6,4,mesh%ncell))
    !allocate(this%w0mask(mesh%nnode))

    !! Set up solver
    allocate(g, this%imap, this%A)
    row_imap => mesh%edge_imap
    allocate(nvars(merge(this%mesh%edge_imap%global_size, 0, is_IOP)))
    nvars = 2
    call this%imap%init(row_imap, nvars)

    call g%init(this%imap)
    ! do e = 1, mesh%nedge_onP ! add cross-diagonal terms
    !   er = 2*(e-1) + 1 ! real component for this edge
    !   ei = 2*(e-1) + 2 ! imaginary component for this edge
    !   call g%add_edge(er, ei)
    !   call g%add_edge(ei, er)
    ! end do
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
          !call g%add_edge(e1r, e2r)
          !call g%add_edge(e1i, e2i)
          call g%add_edge(e1r, e2i)
          call g%add_edge(e1i, e2r)
        end do
      end do
    end do
    call g%add_complete

    call this%A%init(g, take_graph=.true.)

    ! set defaults
    if (.not.params%is_parameter("rel-tol")) call params%set("rel-tol", 0.0_r8)
    if (.not.params%is_parameter("abs-tol")) call params%set("abs-tol", 1e-9_r8)
    !if (.not.params%is_parameter("conv-rate-tol")) call params%set("conv-rate-tol", 1e-4_r8)
    if (.not.params%is_parameter("max-ds-iter")) call params%set("max-ds-iter", 100)
    if (.not.params%is_parameter("max-amg-iter")) call params%set("max-amg-iter", 1000)
    if (.not.params%is_parameter("krylov-method")) call params%set("krylov-method", "gmres")
    !if (.not.params%is_parameter("krylov-method")) call params%set("krylov-method", "bicgstab")
    if (.not.params%is_parameter("gmres-krylov-dim")) call params%set("gmres-krylov-dim", 10)
    call this%solver%init(this%A, params)

    ! !! Create the mask arrays for the spaces: any edge not tagged as a
    ! !! dirichlet (boundary) edge is masked (w1mask); any node on a unmasked
    ! !! edge is unmasked (w0mask).  'Masked' means a true value and marks
    ! !! a DOF that belongs to the space.
    ! this%edge_dirichlet = (emask /= 0)
    ! this%node_dirichlet = .false.
    ! do j = 1, mesh%nedge
    !   if (this%edge_dirichlet(j)) this%node_dirichlet(mesh%enode(:,j)) = .true.
    ! end do
    ! call this%mesh%node_imap%gather_offp(this%node_dirichlet)

    ! !! Assemble the edge-based coefficient matrix A1
    ! block
    !   type(msr_graph), pointer :: g
    !   allocate(g)
    !   call g%init(mesh%nedge)
    !   call g%add_clique(mesh%cedge)
    !   call g%add_complete
    !   call this%a1%init(g, take_graph=.true.)
    ! end block

    ! !! Project out the rows and columns corresponding to Dirichlet edges.
    ! ebedge = pack(array=[(j,j=1,mesh%nedge)], mask=(emask/=0))
    ! call remove_dof(this%a1, ebedge)
    ! !ASSERT(this%a1%is_symmetric())

    ! !! Projected System.  We form the projected system defined on the nullspace
    ! !! of the curl operator.  This corresponds to the range of the gradient
    ! !! operator, and thus the system is representable as node-based system
    ! !! Compute and assemble the node-based projected matrix A0.
    ! block
    !   type(msr_graph), pointer :: g
    !   allocate(g)
    !   call g%init(mesh%nnode)
    !   call g%add_clique(mesh%cnode)
    !   call g%add_complete
    !   call this%a0%init(g, take_graph=.true.)
    ! end block
    ! do j = 1, this%mesh%ncell
    !   a0 = (etasq*eps(j) + 0.5_r8*dt*sigma(j)) &
    !       * matmul(transpose(grad), sym_matmul(W1_matrix_WE(mesh, j), grad))
    !   do k = 1, 4
    !     do i = 1, 4
    !       call this%a0%add_to(mesh%cnode(i,j), mesh%cnode(k,j), a0(i,k))
    !     end do
    !   end do
    ! end do

    ! !! Project out the rows and columns corresponding to Dirichlet edges.
    ! ebedge = pack(array=[(j,j=1,mesh%nnode)], mask=.not.this%w0mask)
    ! call remove_dof(this%a0, ebedge)
    ! !ASSERT(this%a0%is_symmetric()) ! will generally fail due to order-of-operation differences

    ! call this%eslv%init(this, params)

  ! contains

  !   subroutine remove_dof(this, dof)
  !     type(msr_matrix), intent(inout) :: this
  !     integer, intent(in) :: dof(:)
  !     integer :: j
  !     ASSERT(this%nrow == this%ncol)
  !     do j = 1, size(dof)
  !       call this%project_out(dof(j))
  !       this%diag(dof(j)) = 1.0_r8
  !     end do
  !   end subroutine remove_dof

  end subroutine init


  subroutine setup(this, t, epsr, epsi, mu, sigma, omega)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we
    use upper_packed_matrix_procs, only: sym_matmul

    class(em_fd_solver), intent(inout) :: this
    real(r8), intent(in) :: t, epsr(:), epsi(:), mu(:), sigma(:), omega

    real(r8), parameter :: curl(4,6) = reshape(source=[&
        0.0_r8,  0.0_r8,  1.0_r8,  1.0_r8, &
        0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
        0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
        1.0_r8,  0.0_r8,  0.0_r8,  1.0_r8, &
        -1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, &
        1.0_r8,  1.0_r8,  0.0_r8,  0.0_r8], &
        shape=shape(curl))

    ! real(r8), parameter :: grad(6,4) = reshape(source=[& !! Local gradient matrix
    !     -1.0_r8, -1.0_r8, -1.0_r8,  0.0_r8,  0.0_r8,  0.0_r8, &
    !     1.0_r8,  0.0_r8,  0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
    !     0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
    !     0.0_r8,  0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  1.0_r8], &
    !     shape=shape(grad))

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

    call start_timer("setup")

    this%sigma = sigma
    this%epsi = epsi

    ! non-dimensionalization
    omegar = omega * this%L0 / this%c0
    print *, "omegar: ", omegar

    print *, "setup0"

    call this%bc%compute(t)
    call this%A%set_all(0.0_r8)
    this%rhs = 0

    print *, "setup1"

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
          ! ! mtr1 = 0
          ! ! mtr2 = 0
          ! ! mtr3 = 0
          ! ! mtr1 = merge(1.0_r8, 0.0_r8, e2x==e1x)

          ! non-dimensionalized
          mtr1 = ctm2c(e1x, e2x) / mu(j)
          mtr2 = omegar**2  * epsr(j) * m1(l)
          mtr3 = (omegar**2 * epsi(j) - omegar * sigma(j) * this%Z0 * this%L0) * m1(l)
          !mtr1 = 0

          if (this%bc%is_ebc_edge(e2)) then
            ! Apply Dirichlet BCs to the RHS
            this%rhs(e1r) = this%rhs(e1r) - (mtr1 + mtr2) * this%bc%efield(e2)
            this%rhs(e1i) = this%rhs(e1i) - mtr3 * this%bc%efield(e2)
            !this%rhs(e1r) = this%rhs(e1r) - mtr1 * this%bc%efield(e2r)
            !this%rhs(e1r) = this%rhs(e1r) + mtr1 * this%bc%efield(e2i)
          else
            call this%A%add_to(e1r, e2r,  mtr1 - mtr2)
            call this%A%add_to(e1i, e2i, -mtr1 + mtr2)
            call this%A%add_to(e1r, e2i,  mtr3)
            call this%A%add_to(e1i, e2r,  mtr3)
            if (e1x /= e2x) then
              call this%A%add_to(e2r, e1r,  mtr1 - mtr2)
              call this%A%add_to(e2i, e1i, -mtr1 + mtr2)
              call this%A%add_to(e2i, e1r,  mtr3)
              call this%A%add_to(e2r, e1i,  mtr3)
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
      !this%rhs(e1r) = this%rhs(e1r) + omega * this%bc%hsource(e1)
      !this%rhs(e1i) = this%rhs(e1i) - omega * this%bc%hsource(e1)
    end do

    call this%solver%setup
    call stop_timer("setup")

  end subroutine setup


  subroutine solve(this)

    class(em_fd_solver), intent(inout) :: this

    integer :: ierr

    print *, "max |rhs| = ", maxval(abs(this%rhs)), maxval(abs(this%bc%hsource))

    call start_timer("solve")
    call this%solver%solve(this%rhs, this%efield, ierr)
    !this%efield = this%rhs
    call stop_timer("solve")
    call tls_info('  EMFD solve: ' // this%solver%metrics_string())
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

    class(em_fd_solver), intent(in) :: this
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

end module em_fd_solver_type
