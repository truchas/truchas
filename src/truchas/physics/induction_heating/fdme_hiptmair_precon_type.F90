!!
!! FDME_HIPTMAIR_PRECON_TYPE
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

#define ORIGINAL

#include "f90_assert.fpp"

module fdme_hiptmair_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fdme_precon_class
  use fdme_model_type
  use parameter_list_type
  use simpl_mesh_type
  use msr_matrix_type
  use bcsr_matrix_type
  use truchas_timers
  implicit none
  private

  type, extends(fdme_precon), public :: fdme_hiptmair_precon
    private
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(bcsr_matrix) :: An, Ae
    logical, allocatable :: is_ebc_edge(:), is_ebc_node(:)
    type(msr_matrix) :: Bn(2,2), Be(2,2)  ! alternative to An, Ae
    ! persistent workspace for apply
    real(r8), allocatable :: un(:,:), rn(:,:), r(:,:), b(:,:)
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply !=> alt_apply
  end type

contains

  subroutine init(this, model, params)

    class(fdme_hiptmair_precon), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params

    integer :: i

    this%model => model
    this%mesh => model%mesh
    allocate(this%un(2,this%mesh%nnode), this%rn(2,this%mesh%nnode))
    allocate(this%r(2,this%mesh%nedge), this%b(2,this%mesh%nedge))
    allocate(this%is_ebc_node(this%mesh%nnode))

    ! mark Dirichlet BC edges & nodes

    allocate(this%is_ebc_edge(this%mesh%nedge), source=.false.)
    if (allocated(this%model%ebc)) then
      do i = 1, size(this%model%ebc%index)
        this%is_ebc_edge(this%model%ebc%index(i)) = .true.
      end do
    end if
    this%is_ebc_node = .false.
    do i = 1, this%mesh%nedge
      if (this%is_ebc_edge(i)) this%is_ebc_node(this%mesh%enode(:,i)) = .true.
    end do

    !! 2x2 block CSR matrix storage
    block
      type(csr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nedge)
      call g%add_clique(this%mesh%cedge)
      call g%add_complete
      call this%Ae%init(2, g, take_graph=.true.)
    end block

    block
      type(csr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nnode)
      call g%add_clique(this%mesh%cnode)
      call g%add_complete
      call this%An%init(2, g, take_graph=.true.)
    end block

    !! Alternative real/imaginary partitioned CSR block storage.
    block
      type(msr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nedge)
      call g%add_clique(this%mesh%cedge)
      call g%add_complete
      call this%Be(1,1)%init(g, take_graph=.true.)
      call this%Be(2,1)%init(mold=this%Be(1,1))
      call this%Be(1,2)%init(mold=this%Be(1,1))
      call this%Be(2,2)%init(mold=this%Be(1,1))
    end block

    block
      type(msr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nnode)
      call g%add_clique(this%mesh%cnode)
      call g%add_complete
      call this%Bn(1,1)%init(g, take_graph=.true.)
      call this%Bn(2,1)%init(mold=this%Bn(1,1))
      call this%Bn(1,2)%init(mold=this%Bn(1,1))
      call this%Bn(2,2)%init(mold=this%Bn(1,1))
    end block

  end subroutine init


  !! State-dependent setup: compute and assemble the node-based projected matrix An
  subroutine setup(this)

    use mimetic_discretization, only: w1_matrix_we, w2_matrix_we, cell_grad, cell_curl
    use upper_packed_matrix_procs, only: upm_cong_prod

    class(fdme_hiptmair_precon), intent(inout) :: this

    integer :: j
    real(r8) :: omegar, c1, c2, c3
    real(r8) :: m1(21), ctm2c(21), etmp(2,2,21), gtm1g(10), ntmp(2,2,10)
    real(r8), parameter :: ID2(2,2) = reshape([1.0_r8, 0.0_r8, 0.0_r8, 1.0_r8], shape=[2,2])
    real(r8), parameter :: relaxation = 1.0_r8

    call start_timer("precon")

    omegar = this%model%omega/this%model%c0

    call this%Ae%set_all(0.0_r8)
    call this%An%set_all(0.0_r8)

    do j = 1, this%mesh%ncell
      ! non-dimensionalized
      c1 = 1.0_r8 / this%model%mu(j)
      c2 = this%model%epsr(j) * omegar**2
      c3 = this%model%epsi(j) * omegar**2 - omegar * this%model%sigma(j) * this%model%Z0
#ifdef ORIGINAL
      !NB: This results in a different matrix than the actual matrix of the edge-based system.
      !TODO: What is the rationale for the following modification?
      c2 = c2 + relaxation * c3
      c3 = c2
#endif

      m1 = W1_matrix_WE(this%mesh, j)
      ctm2c = upm_cong_prod(4, 6, W2_matrix_WE(this%mesh, j), cell_curl)
      gtm1g = upm_cong_prod(6, 4, m1, cell_grad)

      etmp(1,1,:) = (c1 * ctm2c - c2 * m1)
      etmp(2,2,:) = -(c1 * ctm2c - c2 * m1)
      etmp(1,2,:) = c3 * m1
      etmp(2,1,:) = c3 * m1
      ! Multiply imaginary equation by -1
      !etmp(2,1,:) = -etmp(2,1,:)
      !etmp(2,2,:) = -etmp(2,2,:)
      call this%Ae%add_to(this%mesh%cedge(:,j), etmp)

      ntmp(1,1,:) = -(c2 * gtm1g)
      ntmp(2,2,:) =  (c2 * gtm1g)
      ntmp(1,2,:) =  (c3 * gtm1g)
      ntmp(2,1,:) =  (c3 * gtm1g)
      ! Multiply imaginary equation by -1
      !ntmp(2,1,:) = -ntmp(2,1,:)
      !ntmp(2,2,:) = -ntmp(2,2,:)
      call this%An%add_to(this%mesh%cnode(:,j), ntmp)

      !! Alternative real/imaginary partitioned system
      call this%Be(1,1)%add_to(this%mesh%cedge(:,j), (c1 * ctm2c - c2 * m1))
      call this%Be(2,2)%add_to(this%mesh%cedge(:,j), -(c1 * ctm2c - c2 * m1))
      call this%Be(1,2)%add_to(this%mesh%cedge(:,j), c3*m1)
      call this%Be(2,1)%add_to(this%mesh%cedge(:,j), c3*m1)

      call this%Bn(1,1)%add_to(this%mesh%cnode(:,j), -c2*gtm1g)
      call this%Bn(2,2)%add_to(this%mesh%cnode(:,j), c2*gtm1g)
      call this%Bn(1,2)%add_to(this%mesh%cnode(:,j), c3*gtm1g)
      call this%Bn(2,1)%add_to(this%mesh%cnode(:,j), c3*gtm1g)
    end do

    do j = 1, this%mesh%nedge
      if (this%is_ebc_edge(j)) then
        call this%Ae%project_out(j)
        call this%Ae%set(j, j, ID2)
      end if
    end do

    do j = 1, this%mesh%nnode
      if (this%is_ebc_node(j)) then
        call this%An%project_out(j)
        call this%An%set(j, j, ID2)
      end if
    end do

    do j = 1, this%mesh%nedge
      if (this%is_ebc_edge(j)) then
        call this%Be(1,1)%project_out(j)
        call this%Be(2,1)%project_out(j)
        call this%Be(1,2)%project_out(j)
        call this%Be(2,2)%project_out(j)
        call this%Be(1,1)%set(j, j, 1.0_r8)
        call this%Be(2,2)%set(j, j, 1.0_r8)
      end if
    end do

    do j = 1, this%mesh%nnode
      if (this%is_ebc_node(j)) then
        call this%Bn(1,1)%project_out(j)
        call this%Bn(2,1)%project_out(j)
        call this%Bn(1,2)%project_out(j)
        call this%Bn(2,2)%project_out(j)
        call this%Bn(1,1)%set(j, j, 1.0_r8)
        call this%Bn(2,2)%set(j, j, 1.0_r8)
      end if
    end do

!    block ! test result -- Okay
!      real(r8), dimension(2,this%mesh%nedge) :: x, y1, y2
!      call random_number(x)
!      y1 = this%Ae%matvec(x)
!      y2(1,:) = this%Be(1,1)%matvec(x(1,:)) + this%Be(1,2)%matvec(x(2,:))
!      y2(2,:) = this%Be(2,1)%matvec(x(1,:)) + this%Be(2,2)%matvec(x(2,:))
!      print *, 'EDGE MATVEC ERROR:', count(abs(y1-y2) > 1d-12*abs(y1))
!    end block

!    block ! test result -- Okay
!      real(r8), dimension(2,this%mesh%nnode) :: x, y1, y2
!      call random_number(x)
!      y1 = this%An%matvec(x)
!      y2(1,:) = this%Bn(1,1)%matvec(x(1,:)) + this%Bn(1,2)%matvec(x(2,:))
!      y2(2,:) = this%Bn(2,1)%matvec(x(1,:)) + this%Bn(2,2)%matvec(x(2,:))
!      print *, 'NODE MATVEC ERROR:', count(abs(y1-y2) > 1d-12*abs(y1))
!    end block

    call stop_timer("precon")

  end subroutine setup

  subroutine apply(this, x)

    use mimetic_discretization, only: grad, grad_t
    use msr_matrix_type, only: gs_relaxation

    class(fdme_hiptmair_precon), intent(inout) :: this
    real(r8), intent(inout) :: x(:,:)

    integer :: j

    ASSERT(size(x,2) == this%mesh%edge_imap%onp_size) !FIXME?
    ASSERT(all(.not.this%is_ebc_edge .or. x(1,:) == 0))
    ASSERT(all(.not.this%is_ebc_edge .or. x(2,:) == 0))

    call start_timer('precon')

    this%b(:,:) = x
    x = 0.0_r8

    !! Forward Gauss-Seidel relaxation on the on-process edge system.
    call gs_relaxation(this%Ae, this%b, x, 'f')

    !! Update the local residual and project it to the nodes.
    this%r = this%b - this%Ae%matvec(x)
    call grad_t(this%mesh, this%r(1,:), this%rn(1,:))
    call grad_t(this%mesh, this%r(2,:), this%rn(2,:))
    do j = 1, this%mesh%nnode
      if (this%is_ebc_node(j)) this%rn(:,j) = 0.0_r8
    end do

    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
    this%un = 0.0_r8
    call gs_relaxation(this%An, this%rn, this%un, 'fb')

    !! Update the the solution with the node-based correction.
    call grad(this%mesh, this%un(1,:), x(1,:), increment=.true.)
    call grad(this%mesh, this%un(2,:), x(2,:), increment=.true.)

    !! Backward Gauss-Seidel relaxation on the on-process edge system.
    call gs_relaxation(this%Ae, this%b, x, 'b')

    !TODO: Fix this parallel step
    !call this%mesh%edge_imap%scatter_offp_sum(x)
    !call this%mesh%edge_imap%gather_offp(x)

    call stop_timer('precon')

    ASSERT(all(.not.this%is_ebc_edge .or. x(1,:) == 0))
    ASSERT(all(.not.this%is_ebc_edge .or. x(2,:) == 0))

  end subroutine apply

!  !! An alternate storage scheme for the preconditioner matrices stores them
!  !! as a 2x2 block system according to the real/imaginary part partitioning
!  !! of the unknowns. Each edge-based block has identical non-zero structure.
!  !! Morever there are only two distinct blocks in each matrix, with the other
!  !! two blocks being equal to or the negative of one of the first. This
!  !! allows the storage to be cut in half (not yet exploited). HOWEVER, the
!  !! Gauss-Seidel relaxations, which depend on the ordering of unknowns, are
!  !! fundamentally different in this scheme, and appear to be significantly
!  !! less effective.
!
!  subroutine alt_apply(this, b, x)
!
!    use mimetic_discretization, only: grad, grad_t
!    use msr_matrix_type, only: gs_relaxation
!
!    class(fdme_hiptmair_precon), intent(inout) :: this
!    real(r8), intent(in) :: b(:,:)
!    real(r8), intent(inout) :: x(:,:)
!
!    integer :: j
!
!    ASSERT(size(b,2) == this%mesh%edge_imap%onp_size) !FIXME?
!    ASSERT(size(x,2) == this%mesh%edge_imap%onp_size) !FIXME?
!    ASSERT(all(.not.this%is_ebc_edge .or. b(1,:) == 0))
!    ASSERT(all(.not.this%is_ebc_edge .or. b(2,:) == 0))
!
!    call start_timer('precon')
!
!    !! Forward Gauss-Seidel relaxation on the on-process edge system.
!    x = 0.0_r8
!    call gs_relaxation(this%Be(1,1), b(1,:), x(1,:), 'f')
!    this%r(2,:) = b(2,:) - this%Be(2,1)%matvec(x(1,:))
!    call gs_relaxation(this%Be(2,2), this%r(2,:), x(2,:), 'f')
!
!    !! Update the local residual and project it to the nodes.
!    this%r(1,:) = b(1,:) - this%Be(1,1)%matvec(x(1,:)) - this%Be(1,2)%matvec(x(2,:))
!    this%r(2,:) = b(2,:) - this%Be(2,1)%matvec(x(1,:)) - this%Be(2,2)%matvec(x(2,:))
!    call grad_t(this%mesh, this%r(1,:), this%rn(1,:))
!    call grad_t(this%mesh, this%r(2,:), this%rn(2,:))
!    do j = 1, this%mesh%nnode
!      if (this%is_ebc_node(j)) this%rn(:,j) = 0.0_r8
!    end do
!
!    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
!    this%un = 0.0_r8
!    call gs_relaxation(this%Bn(1,1), this%rn(1,:), this%un(1,:), 'f')
!    this%rn(2,:) = this%rn(2,:) - this%Bn(2,1)%matvec(this%un(1,:))
!    call gs_relaxation(this%Bn(2,2), this%rn(2,:), this%un(2,:), 'fb')
!    this%rn(1,:) = this%rn(1,:) - this%Bn(1,2)%matvec(this%un(2,:))
!    call gs_relaxation(this%Bn(1,1), this%rn(1,:), this%un(1,:), 'b')
!
!    !! Update the the solution with the node-based correction.
!    call grad(this%mesh, this%un(1,:), x(1,:), increment=.true.)
!    call grad(this%mesh, this%un(2,:), x(2,:), increment=.true.)
!
!    !! Backward Gauss-Seidel relaxation on the on-process edge system.
!    call gs_relaxation(this%Be(2,2), b(2,:), x(2,:), 'b')
!    this%r(1,:) = b(1,:) - this%Be(1,2)%matvec(x(2,:))
!    call gs_relaxation(this%Be(1,1), this%r(1,:), x(1,:), 'b')
!
!    !TODO: Fix this parallel step
!    !call this%mesh%edge_imap%scatter_offp_sum(x)
!    !call this%mesh%edge_imap%gather_offp(x)
!
!    call stop_timer('precon')
!
!    ASSERT(all(.not.this%is_ebc_edge .or. x(1,:) == 0))
!    ASSERT(all(.not.this%is_ebc_edge .or. x(2,:) == 0))
!
!  end subroutine alt_apply

end module fdme_hiptmair_precon_type
