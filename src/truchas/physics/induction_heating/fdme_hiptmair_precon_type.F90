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

#include "f90_assert.fpp"

module fdme_hiptmair_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fdme_precon_class
  use fdme_model_type
  use parameter_list_type
  use simpl_mesh_type
  use bsr_matrix_type
  use truchas_timers
  implicit none
  private

  type, extends(fdme_precon), public :: fdme_hiptmair_precon
    private
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(bsr_matrix) :: An, Ae
    logical, allocatable :: is_ebc_edge(:), is_ebc_node(:)
    ! persistent workspace for apply
    real(r8), allocatable :: un(:,:), rn(:,:), r(:,:), b(:,:)
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply
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

#define ORIGINAL
!#define SYMMETRIZE

    do j = 1, this%mesh%ncell
      ! non-dimensionalized
      c1 = 1.0_r8 / this%model%mu(j)
      c2 = this%model%epsr(j) * omegar**2
      c3 = -(this%model%epsi(j) * omegar**2 + omegar * this%model%sigma(j) * this%model%Z0)
print *, c2, -c3
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
      etmp(2,2,:) = etmp(1,1,:)
      etmp(1,2,:) = -c3 * m1 !+ c2 * m1 ! from Elmer doc
      etmp(2,1,:) = -etmp(1,2,:)
#ifdef SYMMETRIZE
      ! Symmetrize the matrix by multiplying the imaginary equation by -1
      etmp(2,1,:) = -etmp(2,1,:)
      etmp(2,2,:) = -etmp(2,2,:)
#endif
      call this%Ae%add_to(this%mesh%cedge(:,j), etmp)

      ntmp(1,1,:) = -(c2 * gtm1g)
      ntmp(2,2,:) = ntmp(1,1,:)
      ntmp(1,2,:) = -(c3 * gtm1g) !+ (c2 * gtm1g) ! from Elmer doc
      ntmp(2,1,:) = -ntmp(1,2,:)
#ifdef SYMMETRIZE
      ! Multiply imaginary equation by -1
      ntmp(2,1,:) = -ntmp(2,1,:)
      ntmp(2,2,:) = -ntmp(2,2,:)
#endif
      call this%An%add_to(this%mesh%cnode(:,j), ntmp)
    end do

    ! LHS contribution from Robin boundary conditions
    !FIXME? contribution to node-based subspace matrix is missing
    !FIXME: ONLY CORRECT FOR MU=1
    if (allocated(this%model%robin_lhs)) then
      block
        use mimetic_discretization, only: w1_face_matrix
        integer :: n
        real(r8) :: m(6), a(2,2,6), gtmg(6)
        real(r8), parameter :: g(3,3) = reshape([0,-1,1,0,-1,0,1,-1,1,0],shape=[3,3])
        !call this%model%robin_lhs%compute(t=0.0_r8)
        do j = 1, size(this%model%robin_lhs%index)
          n = this%model%robin_lhs%index(j)
          m = w1_face_matrix(this%mesh, n)
          a(1,1,:) = -(this%model%robin_lhs%value(j)%re) * m
          a(1,2,:) =  (this%model%robin_lhs%value(j)%im) * m
#ifdef SYMMETRIZE
          a(2,1,:) = a(1,2,:)
          a(2,2,:) = -a(1,1,:)
#else
          a(2,1,:) = -a(1,2,:)
          a(2,2,:) = a(1,1,:)
#endif
          call this%Ae%add_to(this%mesh%fedge(:,n), a)
          
          !NB: This needs to be checked for correctness.
          gtmg = upm_cong_prod(3, 3, m, g)
          a(1,1,:) = -(this%model%robin_lhs%value(j)%re) * gtmg
          a(1,2,:) =  (this%model%robin_lhs%value(j)%im) * gtmg
#ifdef SYMMETRIZE
          a(2,1,:) = a(1,2,:)
          a(2,2,:) = -a(1,1,:)
#else
          a(2,1,:) = -a(1,2,:)
          a(2,2,:) = a(1,1,:)
#endif
          call this%An%add_to(this%mesh%fnode(:,n), a)
        end do
      end block
    end if

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

    call stop_timer("precon")

  end subroutine setup

  subroutine apply(this, x)

    use mimetic_discretization, only: grad, grad_t
    use csr_matrix_type, only: gs_relaxation

    class(fdme_hiptmair_precon), intent(inout) :: this
    real(r8), intent(inout) :: x(:,:)

    integer :: j, nedge_onP, nnode_onP

    nedge_onP = this%mesh%nedge_onP
    nnode_onP = this%mesh%nnode_onP

    ASSERT(size(x,2) == this%mesh%nedge)
!    ASSERT(all(.not.this%is_ebc_edge .or. x(1,:) == 0))
!    ASSERT(all(.not.this%is_ebc_edge .or. x(2,:) == 0))

    call start_timer('precon')

    this%b(:,:) = x
    call this%mesh%edge_imap%gather_offp(this%b)

    !! Forward Gauss-Seidel relaxation on the on-process edge system.
    x = 0.0_r8
    call gs_relaxation(this%Ae, this%b(:,:nedge_onP), x, 'f')

    !! Update the local residual and project it to the nodes.
    this%r = this%b - this%Ae%matvec(x)
    call grad_t(this%mesh, this%r(1,:), this%rn(1,:))
    call grad_t(this%mesh, this%r(2,:), this%rn(2,:))
    do j = 1, this%mesh%nnode
      if (this%is_ebc_node(j)) this%rn(:,j) = 0.0_r8
    end do

    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
    this%un = 0.0_r8
    call gs_relaxation(this%An, this%rn(:,:nnode_onP), this%un, 'fb')

    !! Update the the solution with the node-based correction.
    call grad(this%mesh, this%un(1,:), x(1,:), increment=.true.)
    call grad(this%mesh, this%un(2,:), x(2,:), increment=.true.)

    !! Backward Gauss-Seidel relaxation on the on-process edge system.
    call gs_relaxation(this%Ae, this%b(:,:nedge_onP), x, 'b')

    !TODO: add missing rank-2 version of scatter_offp_sum
    call this%mesh%edge_imap%scatter_offp_sum(x(1,:))
    call this%mesh%edge_imap%scatter_offp_sum(x(2,:))
    call this%mesh%edge_imap%gather_offp(x)

    call stop_timer('precon')

!    ASSERT(all(.not.this%is_ebc_edge .or. x(1,:) == 0))
!    ASSERT(all(.not.this%is_ebc_edge .or. x(2,:) == 0))

  end subroutine apply

end module fdme_hiptmair_precon_type
