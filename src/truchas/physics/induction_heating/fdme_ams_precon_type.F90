!!
!! FDME_AMS_PRECON_TYPE
!!
!! This module defines an Auxilliary-space Maxwell Solver to be used as a
!! preconditioner for the flexible GMRES solver. It is specially designed for
!! the time-harmonic Maxwell equations, leading to a block system, see [1].
!!
!! Zach Jibben <zjibben@lanl.gov>
!!
!! References:
!!
!! 1. Grayver and Kolev, Large-scale 3D geoelectromagnetic modeling using
!! parallel adaptive high-order finite element method, 2015.
!!
!! TODO:
!!
!!   - We should provide the Aalpha and Abeta Poisson matrices to Hypre. This
!!     will allow Hypre to skip the step of internally backing them out from
!!     the full matrix, and save compute time in the setup stage.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fdme_ams_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use fdme_precon_class
  use parameter_list_type
  use truchas_timers
  use fhypre
  use simpl_mesh_type
  use pcsr_matrix_type
  use bndry_func1_class
  use fdme_model_type
  implicit none
  private

  type, extends(fdme_precon), public :: fdme_ams_precon
    private
    type(fdme_model), pointer :: model => null() ! unowned reference
    integer, public :: niter = 0
    integer :: nrows = 0, ilower = 0, iupper = 0
    integer, allocatable :: rows(:)
    type(hypre_obj) :: solver = hypre_null_obj ! HYPRE_Solver object handle
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference

    type(hypre_obj) :: bh = hypre_null_obj ! HYPRE_IJVector object handle
    type(hypre_obj) :: xh = hypre_null_obj ! HYPRE_IJVector object handle
    type(hypre_obj) :: lh = hypre_null_obj ! HYPRE_IJVector object handle

    type(hypre_obj) :: xnh = hypre_null_obj ! HYPRE_IJVector object handle
    type(hypre_obj) :: ynh = hypre_null_obj ! HYPRE_IJVector object handle
    type(hypre_obj) :: znh = hypre_null_obj ! HYPRE_IJVector object handle

    type(pcsr_matrix) :: A
    type(pcsr_matrix) :: grad
    type(hypre_obj) :: Ah = hypre_null_obj ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: gradh = hypre_null_obj ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: Aalpha = hypre_null_obj
    type(hypre_obj) :: Abeta = hypre_null_obj
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply
    final :: fdme_ams_precon_final
  end type fdme_ams_precon

contains

  subroutine fdme_ams_precon_final(this)
    type(fdme_ams_precon), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy(this%Ah, ierr)
    if (hypre_associated(this%gradh)) call fHYPRE_IJMatrixDestroy(this%gradh, ierr)
    if (hypre_associated(this%Aalpha)) call fHYPRE_IJMatrixDestroy(this%Aalpha, ierr)
    if (hypre_associated(this%Abeta)) call fHYPRE_IJMatrixDestroy(this%Abeta, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy(this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy(this%xh, ierr)
    if (hypre_associated(this%lh)) call fHYPRE_IJVectorDestroy(this%lh, ierr)
    if (hypre_associated(this%xnh)) call fHYPRE_IJVectorDestroy(this%xnh, ierr)
    if (hypre_associated(this%ynh)) call fHYPRE_IJVectorDestroy(this%ynh, ierr)
    if (hypre_associated(this%znh)) call fHYPRE_IJVectorDestroy(this%znh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_AMSDestroy(this%solver, ierr)
    INSIST(ierr == 0)
    call fHYPRE_ClearAllErrors
  end subroutine fdme_ams_precon_final


  subroutine init(this, model, params)

    class(fdme_ams_precon), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params

    integer :: ipar, i, ierr, comm
    real(r8) :: rpar

    this%model => model
    this%mesh => model%mesh
    call this%A%init(model%A%graph, take_graph=.false.)
    comm = this%mesh%edge_imap%comm
    this%nrows = this%mesh%edge_imap%onp_size
    this%ilower = this%mesh%edge_imap%first_gid
    this%iupper = this%mesh%edge_imap%last_gid
    this%rows = [ (i, i = this%ilower, this%iupper) ] ! global row indices for this process

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorCreate(comm, this%mesh%node_imap%first_gid, this%mesh%node_imap%last_gid, this%lh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%bh, 0, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%xh, 0, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%lh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_AMSCreate(this%solver, ierr)
    INSIST(ierr == 0)

    !! Parameters
    call params%get('max-iter', ipar, default=1)
    call fHYPRE_AMSSetMaxIter(this%solver, ipar, ierr)
    INSIST(ierr == 0)

    call params%get('tol', rpar, default=0.0_r8)
    call fHYPRE_AMSSetTol(this%solver, rpar, ierr)
    INSIST(ierr == 0)

    call params%get('print-level', ipar, default=0)
    call fHYPRE_AMSSetPrintLevel(this%solver, ipar, ierr)
    INSIST(ierr == 0)

    ! According to Hypre documentation, additive cycles are meant to be
    ! used when AMS is called as a preconditioner. Of the options, 8 is
    ! the only additive cycle in their list of fastest options.
    ! See https://hypre.readthedocs.io/en/latest/solvers-ams.html.
    !call fHYPRE_AMSSetCycleType(this%solver, 8, ierr)
    INSIST(ierr == 0)

    ! This can help, but not clear what it's doing in PC mode
    !call fHYPRE_AMSSetProjectionFrequency(this%solver, 20, ierr)
    !INSIST(ierr == 0)

    !! Geometry-dependent setup
    call set_gradient
    call set_coordinate_vectors

  contains

    subroutine set_gradient()

      integer :: e
      type(pcsr_graph), pointer :: g

      !! set up graph
      allocate(g)
      call g%init(this%mesh%edge_imap, this%mesh%node_imap)
      do e = 1, this%mesh%nedge_onP
        call g%add_edge(e, this%mesh%enode(:,e))
      end do
      call g%add_complete

      !! set up matrix
      call this%grad%init(g, take_graph=.true.)
      call this%grad%set_all(0.0_r8)
      do e = 1, this%mesh%nedge_onP
        call this%grad%set(e, this%mesh%enode(1,e), -1.0_r8)
        call this%grad%set(e, this%mesh%enode(2,e), +1.0_r8)
      end do

      !! copy to hypre
      call this%grad%copy_to_ijmatrix(this%gradh)
      call fHYPRE_AMSSetDiscreteGradient(this%solver, this%gradh, ierr)
      INSIST(ierr == 0)

    end subroutine set_gradient


    subroutine set_coordinate_vectors()

      integer :: nrows, ilower, iupper, i
      integer, allocatable :: rows(:)
      real(r8), allocatable :: x(:), y(:), z(:)

      nrows = this%mesh%node_imap%onp_size
      ilower = this%mesh%node_imap%first_gid
      iupper = this%mesh%node_imap%last_gid
      rows = [(i, i = ilower, iupper)] ! global row indices for this process

      x = this%mesh%x(1,:this%mesh%nnode_onP)
      y = this%mesh%x(2,:this%mesh%nnode_onP)
      z = this%mesh%x(3,:this%mesh%nnode_onP)

      call fHYPRE_IJVectorCreate(comm, ilower, iupper, this%xnh, ierr)
      call fHYPRE_IJVectorCreate(comm, ilower, iupper, this%ynh, ierr)
      call fHYPRE_IJVectorCreate(comm, ilower, iupper, this%znh, ierr)
      call fHYPRE_IJVectorSetMaxOffProcElmts(this%xnh, 0, ierr)
      call fHYPRE_IJVectorSetMaxOffProcElmts(this%ynh, 0, ierr)
      call fHYPRE_IJVectorSetMaxOffProcElmts(this%znh, 0, ierr)
      INSIST(ierr == 0)

      call fHYPRE_IJVectorInitialize(this%xnh, ierr)
      call fHYPRE_IJVectorInitialize(this%ynh, ierr)
      call fHYPRE_IJVectorInitialize(this%znh, ierr)
      call fHYPRE_IJVectorSetValues(this%xnh, nrows, rows, x, ierr)
      call fHYPRE_IJVectorSetValues(this%ynh, nrows, rows, y, ierr)
      call fHYPRE_IJVectorSetValues(this%znh, nrows, rows, z, ierr)
      call fHYPRE_IJVectorAssemble(this%xnh, ierr)
      call fHYPRE_IJVectorAssemble(this%ynh, ierr)
      call fHYPRE_IJVectorAssemble(this%znh, ierr)
      INSIST(ierr == 0)

      call fHYPRE_AMSSetCoordinateVectors(this%solver, this%xnh, this%ynh, this%znh, ierr)
      INSIST(ierr == 0)

    end subroutine set_coordinate_vectors

  end subroutine init


  !! State-dependent setup
  subroutine setup(this)

    class(fdme_ams_precon), intent(inout) :: this

    integer :: j, n, ierr
    real(r8) :: alpha(this%mesh%ncell), beta(this%mesh%ncell)

    !ASSERT(all(ieee_is_finite(this%A%values)))

    call start_timer("precon")

    !! Preconditioner matrix
    this%A%values(:) = this%model%A%values(1,1,:) - this%model%A%values(1,2,:)
    if (allocated(this%model%ebc)) then
      do j = 1, size(this%model%ebc%index)
        n = this%model%ebc%index(j)
        call this%A%set(n, n, 1.0_r8)
      end do
    end if

    !! Provide list of nodes inside the 0-dissipation region (tagged with 1)
    block
      integer :: nrow
      integer, allocatable :: rows(:)
      real(r8), allocatable :: interior_nodes(:)
      allocate(interior_nodes(this%mesh%nnode), source=1.0_r8)
      associate (epsi => this%model%epsi, sigma => this%model%sigma)
        do j = 1, this%mesh%ncell
          if (max(sigma(j), epsi(j)) > 0) interior_nodes(this%mesh%cnode(:,j)) = 0
        end do
      end associate
      nrow = this%mesh%node_imap%onp_size
      rows = [(j, j=this%mesh%node_imap%first_gid, this%mesh%node_imap%last_gid)]
      call fHYPRE_ClearAllErrors
      call fHYPRE_IJVectorInitialize(this%lh, ierr)
      call fHYPRE_IJVectorSetValues(this%lh, nrow, rows, interior_nodes(:nrow), ierr)
      call fHYPRE_IJVectorAssemble(this%lh, ierr)
      INSIST(ierr == 0)
      call fHYPRE_AMSSetInteriorNodes(this%solver, this%lh, ierr)
      INSIST(ierr == 0)
    end block

    block ! Optional Poisson matrix construction: Aalpha, Abeta
      integer :: tag(this%mesh%nnode)
      integer, allocatable :: ebc_nodes(:)

      ! The nxE BC are the essential BC for the edge-based A. The Hypre AMS
      ! documentation states, "It is assumed that the essential BC of A, Abeta
      ! and Aalpha are on the same part of the boundary." Now Aalpha and Abeta
      ! are node-based matrices, and so we choose their essential BC to be
      ! associated with those nodes that belong to one of the essential BC
      ! edges of A.
      tag = 0
      do j = 1, size(this%model%ebc%index) ! these are the nxE BC
        associate (enode => this%mesh%enode(:,this%model%ebc%index(j)))
          tag(enode) = enode
        end associate
      end do
      ebc_nodes = pack(tag, mask=(tag > 0))

      ! NB: Per the Hypre AMS documentation, these coefficients should match
      ! what goes into A above.
      associate (model => this%model, kappa => this%model%omega/this%model%c0)
        do j = 1, this%mesh%ncell
          alpha(j) = 1.0_r8/model%mu(j)
          beta(j) = model%sigma(j)*kappa*model%Z0 + model%epsi(j)*kappa**2 - model%epsr(j)*kappa**2
        end do
      end associate

      call create_poisson_matrices(this%mesh, ebc_nodes, alpha, beta, this%Aalpha, this%Abeta)
    end block
!    call fHYPRE_AMSSetAlphaPoissonMatrix(this%solver, this%Aalpha, ierr)
!    INSIST(ierr == 0)
!    call fHYPRE_AMSSetBetaPoissonMatrix(this%solver, this%Abeta, ierr)
!    INSIST(ierr == 0)

    call this%A%copy_to_ijmatrix(this%Ah)
    call fHYPRE_AMSSetup(this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

    call stop_timer("precon")

  end subroutine setup


  subroutine apply(this, x)

    class(fdme_ams_precon), intent(inout) :: this
    real(r8), intent(inout) :: x(:,:)

    integer :: ierr, stat

    ASSERT(size(x) == 2*this%nrows)

    call start_timer("precon")

    call fHYPRE_ClearAllErrors
    this%niter = 0

    call apply_system(this%solver, this%Ah, x(1,:))
    ASSERT(all(ieee_is_finite(x(1,:))))

    call apply_system(this%solver, this%Ah, x(2,:))
    ASSERT(all(ieee_is_finite(x(2,:))))

    call stop_timer("precon")

  contains

    subroutine apply_system(solver, A, x)

      use,intrinsic :: ieee_arithmetic, only: ieee_is_finite

      type(hypre_obj), intent(in) :: solver, A
      real(r8), intent(inout) :: x(:)

      real(r8) :: rel_resid_norm
      integer :: num_iterations

      ASSERT(all(ieee_is_finite(x)))

      !! Initialize the Hypre RHS vector.
      call fHYPRE_IJVectorInitialize(this%bh, ierr)
      call fHYPRE_IJVectorSetValues(this%bh, this%nrows, this%rows, x, ierr)
      call fHYPRE_IJVectorAssemble(this%bh, ierr)
      INSIST(ierr == 0)

      !! Initialize the Hypre initial guess vector.
      x = 0.0_r8  !TODO: is there another way to set the hypre vector to 0?
      call fHYPRE_IJVectorInitialize(this%xh, ierr)
      call fHYPRE_IJVectorSetValues(this%xh, this%nrows, this%rows, x, ierr)
      call fHYPRE_IJVectorAssemble(this%xh, ierr)
      INSIST(ierr == 0)

      !! Solve the system.
      call start_timer("AMSSolve")
      call fHYPRE_AMSSolve(solver, A, this%bh, this%xh, ierr)
      if (ierr /= 0) then
        INSIST(ierr == HYPRE_ERROR_CONV)
        stat = 1
        call fHYPRE_ClearAllErrors
      else
        stat = 0
      end if
      call stop_timer("AMSSolve")

      !! Retrieve the solution vector from HYPRE.
      !call fHYPRE_AMSProjectOutGradients(this%solver, this%xh, ierr) SEGFAULTS
      !INSIST(ierr == 0)
      call fHYPRE_IJVectorGetValues(this%xh, this%nrows, this%rows, x, ierr)
      ASSERT(all(ieee_is_finite(x)))
      INSIST(ierr == 0)

      !! TODO use hypre routine to project out gradient. Might want to save this for the final
      !! solution?

      call fHYPRE_AMSGetNumIterations(solver, num_iterations, ierr)
      call fHYPRE_AMSGetFinalRelativeResidualNorm(solver, rel_resid_norm, ierr)
      INSIST(ierr == 0)
      this%niter = this%niter + num_iterations

      !! DEBUGGING
      ASSERT(ieee_is_finite(rel_resid_norm)) !! WARN: this should never be tripping, let alone sporadically
      !print '("    AMS precon iter, rel_norm: ",i6,es13.3)', num_iterations, rel_resid_norm

    end subroutine apply_system

  end subroutine apply

  !! See https://hypre.readthedocs.io/en/latest/solvers-ams.html#
  !! for the description of the Poisson matrices that AMS can use.

  subroutine create_poisson_matrices(mesh, ebc_nodes, alpha, beta, A_alpha_h, A_beta_h)

    use simplex_geometry, only: tet_face_normal

    type(simpl_mesh), intent(in), target :: mesh
    integer, intent(in) :: ebc_nodes(:)
    real(r8), intent(in) :: alpha(:), beta(:)
    type(hypre_obj), intent(inout) :: A_alpha_h, A_beta_h

    integer :: i, j, k, l, n
    real(r8) :: p(3,4), tmp1(10), tmp2(10)
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix) :: A_alpha, A_beta

    ASSERT(size(alpha) == mesh%ncell)
    ASSERT(size(beta) == mesh%ncell)

    allocate(g)
    call g%init(mesh%node_imap)
    call g%add_clique(mesh%cnode)
    call g%add_complete
    call A_alpha%init(g, take_graph=.true.)
    call A_beta%init(mold=A_alpha)

    do i = 1, mesh%ncell
      p = tet_face_normal(mesh%x(:,mesh%cnode(:,i)))
      l = 0
      do k = 1, 4
        do j = 1, k
          l = l + 1
          tmp1(l) = (alpha(i)/(9*abs(mesh%volume(i)))) * dot_product(p(:,j), p(:,k))
          tmp2(l) =  (beta(i)/(9*abs(mesh%volume(i)))) * dot_product(p(:,j), p(:,k))
        end do
      end do
      call A_beta%add_to(mesh%cnode(:,i), tmp2)

      l = 0
      do k = 1, 4
        do j = 1, k-1
          l = l + 1
          tmp1(l) = tmp1(l) + (beta(i)*abs(mesh%volume(i))/20.0_r8)
        end do
        l = l + 1
        tmp1(l) = tmp1(l) + (beta(i)*abs(mesh%volume(i))/10.0_r8)
      end do
      call A_alpha%add_to(mesh%cnode(:,i), tmp1)
    end do

    !! Modifications due to essential BC
    do i = 1, size(ebc_nodes)
      n = ebc_nodes(i)
      call A_alpha%project_out(n)
      call A_alpha%set(n, n, 1.0_r8)
      call A_beta%project_out(n)
      call A_beta%set(n, n, 1.0_r8)
    end do

    call A_alpha%copy_to_ijmatrix(A_alpha_h)
    call A_beta%copy_to_ijmatrix(A_beta_h)

  end subroutine create_poisson_matrices

end module fdme_ams_precon_type
