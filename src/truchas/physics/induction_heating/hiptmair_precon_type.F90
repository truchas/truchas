!!
!! HIPTMAIR_PRECON_TYPE
!!
!! This module implements a preconditioner [1] modified for the block
!! time-harmonic Maxwell equations [2].
!!
!! Zach Jibben <zjibben@lanl.gov>
!!
!! References:
!!
!! 1. Hiptmair, Multigrid Method For Maxwell's Equations, 1998.
!!
!! 2. Grayver and Kolev, Large-scale 3D geoelectromagnetic modeling using
!! parallel adaptive high-order finite element method, 2015.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module hiptmair_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  !use,intrinsic :: iso_fortran_env, only: r8 => real128
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use parameter_list_type
  use truchas_timers
  use simpl_mesh_type
  use index_map_type
  use pcsr_matrix_type
  use em_bc_type
  implicit none
  private

  type, public :: hiptmair_precon
    private
    integer, public :: niter = 0
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(index_map), pointer :: node_imap => null(), edge_imap => null()
    type(pcsr_matrix) :: An, Ae
    real(r8), allocatable :: un(:), rn(:), r(:), b(:), grad(:,:), curl(:,:)
    logical, allocatable :: is_ebc_edge(:), is_ebc_node(:)

    ! Non-dimensionalization parameters
    real(r8) :: Z0 = 376.730313668 ! Ohms -- vacuum impedance
    real(r8) :: c0 = 299792458.0_r8 ! speed of light
    real(r8) :: L0, H0
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply
  end type hiptmair_precon

contains

  subroutine hiptmair_precon_final(this)
    type(hiptmair_precon), intent(inout) :: this
    if (associated(this%node_imap)) deallocate(this%node_imap)
    if (associated(this%edge_imap)) deallocate(this%edge_imap)
  end subroutine hiptmair_precon_final


  subroutine init(this, params, mesh, A, bc)

    use parallel_communication, only: is_IOP

    class(hiptmair_precon), intent(out) :: this
    type(parameter_list), pointer, intent(in) :: params
    type(simpl_mesh), intent(in), target :: mesh
    type(pcsr_matrix), pointer, intent(in) :: A !! taking ownership
    type(em_bc), intent(in), target :: bc

    integer :: i

    this%mesh => mesh
    allocate(this%un(2*mesh%nnode), this%rn(2*mesh%nnode), this%r(2*mesh%nedge), &
        this%b(2*mesh%nedge), this%is_ebc_node(mesh%nnode))

    this%grad = reshape([-1.0_r8, -1.0_r8, -1.0_r8,  0.0_r8,  0.0_r8,  0.0_r8, &
        &                 1.0_r8,  0.0_r8,  0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
        &                 0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
        &                 0.0_r8,  0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  1.0_r8], &
        shape=[6, 4])

    this%curl = reshape([ 0.0_r8,  0.0_r8,  1.0_r8,  1.0_r8, &
        &                 0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
        &                 0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
        &                 1.0_r8,  0.0_r8,  0.0_r8,  1.0_r8, &
        &                -1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, &
        &                 1.0_r8,  1.0_r8,  0.0_r8,  0.0_r8], &
        shape=[4, 6])

    ! mark Dirichlet BC edges & nodes
    this%is_ebc_edge = bc%is_ebc_edge
    this%is_ebc_node = .false.
    do i = 1, mesh%nedge
      if (bc%is_ebc_edge(i)) this%is_ebc_node(mesh%enode(:,i)) = .true.
    end do

    ! prepare full edge matrix & node-based projection matrix
    call init_matrix(this%Ae, this%edge_imap, mesh%edge_imap, mesh%cedge)
    call init_matrix(this%An, this%node_imap, mesh%node_imap, mesh%cnode)

  contains

    subroutine init_matrix(this_A, this_imap, mesh_imap, adjcncy)

      type(pcsr_matrix), intent(out) :: this_A
      type(index_map), pointer :: this_imap
      type(index_map), intent(in), pointer :: mesh_imap
      integer, intent(in) :: adjcncy(:,:)

      integer :: n1, n1x, n1r, n1i
      integer :: n2, n2x, n2r, n2i
      integer, allocatable :: nvars(:)
      type(index_map), pointer :: row_imap => null()
      type(pcsr_graph), pointer :: g => null()

      allocate(g, this_imap)
      row_imap => mesh_imap
      allocate(nvars(merge(mesh_imap%global_size, 0, is_IOP)))
      nvars = 2
      call this_imap%init(row_imap, nvars)

      call g%init(this_imap)
      do i = 1, this%mesh%ncell
        ! add curl terms
        call g%add_clique(2*(adjcncy(:,i)-1) + 1) ! real components
        call g%add_clique(2*(adjcncy(:,i)-1) + 2) ! imaginary components

        ! add cross-diagonal terms
        do n1x = 1, size(adjcncy(:,i))
          n1 = adjcncy(n1x,i)
          n1r = 2*(n1-1) + 1
          n1i = 2*(n1-1) + 2
          do n2x = 1, n1x
            n2 = adjcncy(n2x,i)
            n2r = 2*(n2-1) + 1
            n2i = 2*(n2-1) + 2
            call g%add_edge(n1r, n2i)
            call g%add_edge(n1i, n2r)
          end do
        end do
      end do
      call g%add_complete
      call this_A%init(g, take_graph=.true.)

    end subroutine init_matrix

  end subroutine init


  !! State-dependent setup: compute and assemble the node-based projected matrix An
  subroutine setup(this, mu, epsr, epsi, sigma, omega)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we
    use upper_packed_matrix_procs, only: sym_matmul

    class(hiptmair_precon), intent(inout) :: this
    real(r8), intent(in) :: mu(:), epsr(:), epsi(:), sigma(:), omega

    real(r8), parameter :: relaxation = 1.0_r8
    integer :: ierr, j, l
    integer :: n1, n1x, n1r, n1i
    integer :: n2, n2x, n2r, n2i
    integer :: e1, e1x, e1r, e1i
    integer :: e2, e2x, e2r, e2i
    real(r8) :: ct_m2_c(6,6), gt_m1_g(4,4), m1(21), mtr1, mtr2, mtr3
    real(r8), pointer :: values => null()
    integer, pointer :: indices => null()

    call start_timer("precon")

    this%L0 = 1

    call this%Ae%set_all(0.0_r8)
    call this%An%set_all(0.0_r8)
    do j = 1, this%mesh%ncell_onP
      ! TODO: precompute these lines, since they're only dependent on geometry?
      m1 = W1_matrix_WE(this%mesh, j)
      ct_m2_c = matmul(transpose(this%curl), sym_matmul(W2_matrix_WE(this%mesh, j), this%curl))
      gt_m1_g = matmul(transpose(this%grad), sym_matmul(W1_matrix_WE(this%mesh, j), this%grad))

      ! non-dimensionalized
      mtr1 = mu(j)
      mtr2 = omega**2 * epsr(j)
      mtr3 = omega**2 * epsi(j) - omega * sigma(j) * this%Z0 * this%L0
      mtr2 = mtr2 + relaxation * mtr3
      mtr3 = mtr2
      !mtr1 = huge(1.0_r8)

      ! full edge-based matrix
      l = 0
      do e2x = 1, 6
        e2 = this%mesh%cedge(e2x,j)
        e2r = 2*(e2-1) + 1
        e2i = 2*(e2-1) + 2

        ! project out the rows and columns corresponding to Dirichlet edges
        if (this%is_ebc_edge(e2)) then
          call this%Ae%set(e2r, e2r, 1.0_r8)
          call this%Ae%set(e2i, e2i, 1.0_r8)
          l = l + e2x
          cycle
        end if

        do e1x = 1, e2x
          l = l + 1
          e1 = this%mesh%cedge(e1x,j)
          e1r = 2*(e1-1) + 1
          e1i = 2*(e1-1) + 2
          if (this%is_ebc_edge(e1)) cycle

          call this%Ae%add_to(e1r, e2r,  ct_m2_c(e1x,e2x) / mtr1 - mtr2 * m1(l))
          call this%Ae%add_to(e1i, e2i, -ct_m2_c(e1x,e2x) / mtr1 + mtr2 * m1(l))
          call this%Ae%add_to(e1r, e2i,  mtr3 * m1(l))
          call this%Ae%add_to(e1i, e2r,  mtr3 * m1(l))
          if (e1x /= e2x) then
            call this%Ae%add_to(e2r, e1r,  ct_m2_c(e1x,e2x) / mtr1 - mtr2 * m1(l))
            call this%Ae%add_to(e2i, e1i, -ct_m2_c(e1x,e2x) / mtr1 + mtr2 * m1(l))
            call this%Ae%add_to(e2i, e1r,  mtr3 * m1(l))
            call this%Ae%add_to(e2r, e1i,  mtr3 * m1(l))
          end if
        end do
      end do

      ! node-based projection matrix
      do n2x = 1, 4
        n2 = this%mesh%cnode(n2x,j)
        n2r = 2*(n2-1) + 1
        n2i = 2*(n2-1) + 2
        ! project out the rows and columns corresponding to Dirichlet edges
        if (this%is_ebc_node(n2)) then
          call this%An%set(n2r, n2r, 1.0_r8)
          call this%An%set(n2i, n2i, 1.0_r8)
          cycle
        end if

        do n1x = 1, 4
          n1 = this%mesh%cnode(n1x,j)
          n1r = 2*(n1-1) + 1
          n1i = 2*(n1-1) + 2
          if (this%is_ebc_node(n1)) cycle
          call this%An%add_to(n1r, n2r,  -mtr2 * gt_m1_g(n1x, n2x))
          call this%An%add_to(n1i, n2i,   mtr2 * gt_m1_g(n1x, n2x))
          call this%An%add_to(n1r, n2i,   mtr3 * gt_m1_g(n1x, n2x))
          call this%An%add_to(n1i, n2r,   mtr3 * gt_m1_g(n1x, n2x))
        end do
      end do
    end do

    ! diagonals used in Gauss-Seidel relaxation
    call this%Ae%update_diag
    call this%An%update_diag

    print *, "completed setup"

    call stop_timer("precon")

  end subroutine setup


  subroutine apply(this, b, x, stat)

    use mimetic_discretization, only: grad, grad_t

    class(hiptmair_precon), intent(inout) :: this
    real(r8), intent(in) :: b(:) ! residual
    real(r8), intent(inout) :: x(:) ! state
    integer, intent(out) :: stat

    integer :: n, nr, ni, nnode, nedge, nnode_onP, nedge_onP

    ASSERT(size(b) == 2*this%mesh%edge_imap%onp_size)
    ASSERT(size(x) == 2*this%mesh%edge_imap%onp_size)
    ASSERT(all(ieee_is_finite(x)))
    ASSERT(all(ieee_is_finite(b)))
    ASSERT(all(.not.this%is_ebc_edge .or. b(1:2*this%mesh%nedge-1:2) == 0))
    ASSERT(all(.not.this%is_ebc_edge .or. b(2:2*this%mesh%nedge:2) == 0))

    nnode = this%mesh%nnode
    nedge = this%mesh%nedge
    nnode_onP = this%mesh%nnode_onP
    nedge_onP = this%mesh%nedge_onP

    stat = 0
    this%niter = 1
    call start_timer("precon")

    !! Forward Gauss-Seidel relaxation on the on-process edge system.
    x = 0
    call gauss_seidel_relaxation(this%Ae, b, x, "f")
    ! print *, 'debug 0: ', maxval(abs(this%b))
    ! print *, 'debug 0: ', maxval(abs(x))

    !! Update the local residual and project it to the nodes.
    call this%Ae%matvec(x, this%r)
    !print *, 'debug 0: ', maxval(abs(this%r))
    this%r = b - this%r
    call grad_t(this%mesh, this%r(1:2*nedge-1:2), this%rn(1:2*nnode-1:2))
    call grad_t(this%mesh, this%r(2:2*nedge:2), this%rn(2:2*nnode:2))
    !print *, 'debug 0: ', maxval(abs(this%rn))
    do n = 1, nnode
      if (this%is_ebc_node(n)) then
        nr = 2*(n-1) + 1
        ni = 2*(n-1) + 2
        this%rn(nr) = 0
        this%rn(ni) = 0
      end if
    end do
    ASSERT(all(ieee_is_finite(this%rn)))

    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
    this%un = 0
    !print *, 'debug 1: ', maxval(abs(this%rn))
    call gauss_seidel_relaxation(this%An, this%rn, this%un, "fb", verbose=.true.)
    !call sor_relaxation(this%An, this%rn(:2*this%mesh%nnode_onP), this%un, "fb", 0.5_r8, verbose=.true.)
    !call jacobi_relaxation(this%An, this%rn(:2*this%mesh%nnode_onP), this%un, 10, verbose=.true.)
    !print *, "debug 1"
    ASSERT(all(ieee_is_finite(this%un)))

    !! Update the the solution with the node-based correction.
    block
      integer :: e, er, ei, n1, n1r, n1i, n2, n2r, n2i
      do e = 1, this%mesh%nedge
        n1 = this%mesh%enode(1,e)
        n2 = this%mesh%enode(2,e)
        er = 2*(e-1) + 1
        ei = 2*(e-1) + 2
        n1r = 2*(n1-1) + 1
        n1i = 2*(n1-1) + 2
        n2r = 2*(n2-1) + 1
        n2i = 2*(n2-1) + 2

        x(er) = x(er) + this%un(n2r) - this%un(n1r) ! grad
        x(ei) = x(ei) + this%un(n2i) - this%un(n1i) ! grad
      end do
    end block
    ! call grad(this%mesh, this%un(1:2*nnode-1:2), x(1:2*nedge-1:2), increment=.true.)
    ! call grad(this%mesh, this%un(2:2*nnode:2), x(2:2*nedge:2), increment=.true.)
    ASSERT(all(ieee_is_finite(x)))

    !! Backward Gauss-Seidel relaxation on the on-process edge system.
    call gauss_seidel_relaxation(this%Ae, b, x, "b")
    ASSERT(all(ieee_is_finite(x)))

    call this%mesh%edge_imap%scatter_offp_sum(x)
    call this%mesh%edge_imap%gather_offp(x)

    call stop_timer("precon")

    ASSERT(all(ieee_is_finite(x)))
    ASSERT(all(.not.this%is_ebc_edge .or. x(1:2*this%mesh%nedge-1:2) == 0))
    ASSERT(all(.not.this%is_ebc_edge .or. x(2:2*this%mesh%nedge:2) == 0))

  end subroutine apply


  subroutine gauss_seidel_relaxation(A, f, u, pattern, verbose)

    type(pcsr_matrix), intent(in) :: A
    real(r8), intent(in) :: f(:)
    real(r8), intent(inout) :: u(:)
    character(*), intent(in) :: pattern
    logical, intent(in), optional :: verbose

    logical :: verbose_
    integer :: i, j, k, imin, imax, di
    real(r8) :: s
    integer, pointer :: indices(:) => null()
    real(r8), pointer :: values(:) => null()

    call start_timer("gauss-seidel")

    !ASSERT(A%nrow == A%ncol)
    ASSERT(size(f) <= A%nrow)
    !ASSERT(size(u) >= A%ncol)
    ASSERT(size(u) >= A%nrow)
    ASSERT(all(ieee_is_finite(f)))
    ASSERT(all(ieee_is_finite(u)))

    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    do j = 1, len(pattern)
      call direction(pattern(j:j), A%nrow, imin, imax, di)
      do i = imin, imax, di
        ASSERT(A%diag(i) /= 0)
        call A%get_row_view(i, values, indices)
        ASSERT(all(ieee_is_finite(values)))
        s = f(i)
        do k = 1, size(indices)
          if (indices(k) /= i) s = s - values(k) * u(indices(k))
        end do
        u(i) = s / A%diag(i)
        ASSERT(ieee_is_finite(u(i)))
      end do
    end do

    call stop_timer("gauss-seidel")

  end subroutine gauss_seidel_relaxation


  pure subroutine direction(pattern, len, imin, imax, di)
    character(1), intent(in) :: pattern
    integer, intent(in) :: len
    integer, intent(out) :: imin, imax, di
    select case (pattern)
    case ('f', 'F') ! forward sweep
      imin = 1
      imax = len
      di = 1
    case ('b', 'B') ! backward sweep
      imin = len
      imax = 1
      di = -1
    end select
  end subroutine direction


  subroutine sor_relaxation(A, f, u, pattern, omega, verbose)

    type(pcsr_matrix), intent(in) :: A
    real(r8), intent(in) :: f(:)
    real(r8), intent(inout) :: u(:)
    character(*), intent(in) :: pattern
    real(r8), intent(in) :: omega
    logical, intent(in), optional :: verbose

    logical :: verbose_
    integer :: i, j, k, imin, imax, di
    real(r8) :: s
    integer, pointer :: indices(:) => null()
    real(r8), pointer :: values(:) => null()

    !ASSERT(A%nrow == A%ncol)
    ASSERT(size(f) <= A%nrow)
    !ASSERT(size(u) >= A%ncol)
    ASSERT(size(u) >= A%nrow)
    ASSERT(all(ieee_is_finite(f)))
    ASSERT(all(ieee_is_finite(u)))

    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    do j = 1, len(pattern)
      call direction(pattern(j:j), size(f), imin, imax, di)
      do i = imin, imax, di
        call A%get_row_view(i, values, indices)
        s = f(i)
        do k = 1, size(indices)
          if (indices(k) == i) cycle
          s = s - values(k) * u(indices(k))
        end do
        u(i) = (1 - omega) * u(i) + omega * s / A%diag(i)
        ASSERT(ieee_is_finite(u(i)))
      end do
    end do

  end subroutine sor_relaxation


  subroutine jacobi_relaxation(A, f, u, niter, verbose)

    type(pcsr_matrix), intent(in) :: A
    real(r8), intent(in) :: f(:)
    real(r8), intent(inout) :: u(:)
    integer, intent(in) :: niter
    logical, intent(in), optional :: verbose

    logical :: verbose_
    integer :: i, j, k, n
    real(r8) :: s, un(size(u))
    integer, pointer :: indices(:) => null()
    real(r8), pointer :: values(:) => null()

    !ASSERT(A%nrow == A%ncol)
    ASSERT(size(f) <= A%nrow)
    !ASSERT(size(u) >= A%ncol)
    ASSERT(size(u) >= A%nrow)
    ASSERT(all(ieee_is_finite(f)))
    ASSERT(all(ieee_is_finite(u)))

    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    do n = 1, niter
      do i = 1, size(f)
        ASSERT(A%diag(i) /= 0)
        call A%get_row_view(i, values, indices)
        ASSERT(all(ieee_is_finite(values)))
        s = f(i)
        do k = 1, size(indices)
          if (indices(k) == i) cycle
          s = s - values(k) * u(indices(k))
        end do
        un(i) = s / A%diag(i)
        ASSERT(ieee_is_finite(u(i)))
      end do
      u = un
    end do

  end subroutine jacobi_relaxation

end module hiptmair_precon_type
