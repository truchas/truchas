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

    real(r8), allocatable :: interior_nodes(:) ! list of nodes in the zero-conductivity region
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

    integer :: ipar, i, ierr
    real(r8) :: rpar

    this%model => model
    this%mesh => model%mesh
    call this%A%init(mold=model%A(1,1))
    this%nrows = this%mesh%edge_imap%onp_size
    this%ilower = this%mesh%edge_imap%first_gid
    this%iupper = this%mesh%edge_imap%last_gid
    this%rows = [ (i, i = this%ilower, this%iupper) ] ! global row indices for this process
    allocate(this%interior_nodes(this%mesh%nnode))

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%lh, ierr)
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

      call fHYPRE_IJVectorCreate(ilower, iupper, this%xnh, ierr)
      call fHYPRE_IJVectorCreate(ilower, iupper, this%ynh, ierr)
      call fHYPRE_IJVectorCreate(ilower, iupper, this%znh, ierr)
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
    real(r8) :: mtr2

    !ASSERT(all(ieee_is_finite(this%A%values)))

    call start_timer("precon")

    !! Preconditioner matrix
    this%A%values(:) = this%model%A(1,1)%values - this%model%A(1,2)%values
    !this%A%diag = this%model%A(1,1)%diag - this%model%A(1,2)%diag
    !this%A%nonz = this%model%A(1,1)%nonz - this%model%A(1,2)%nonz
    if (allocated(this%model%ebc)) then
      do j = 1, size(this%model%ebc%index)
        n = this%model%ebc%index(j)
        call this%A%set(n, n, 1.0_r8)
      end do
    end if

!NNC: I don't think this is necessary, as the AMS system is non-singular,
!     but if needed, it is a node-based vector not edge-based -- FIXME
!     But if it is needed, it needs to be fixed
!    ! List interior nodes for the AMS preconditioner. All nodes inside the
!    ! 0-conductivity region are marked 1.0.
!    associate (omega => this%model%omega, epsr => this%model%epsr, epsi => this%model%epsi, sigma => this%model%sigma)
!    this%interior_nodes = 1
!    do j = 1, this%mesh%ncell
!      ! WARN: Which of the following is correct? I think the first.
!      mtr2 = omega * (epsr(j) + epsi(j)) - sigma(j)
!      !mtr2 = omega * (epsi(j) - epsr(j)) - sigma(j)
!      if (mtr2 /= 0) this%interior_nodes(this%mesh%cnode(:,j)) = 0
!    end do
!    end associate
!
!    call fHYPRE_ClearAllErrors
!
!    ! Provide list of nodes inside the 0-conductivity region
!    call fHYPRE_IJVectorInitialize(this%lh, ierr)
!    call fHYPRE_IJVectorSetValues(this%lh, this%nrows, this%rows, this%interior_nodes, ierr)
!    call fHYPRE_IJVectorAssemble(this%lh, ierr)
!    INSIST(ierr == 0)
!    call fHYPRE_AMSSetInteriorNodes(this%solver, this%lh, ierr)
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
    real(r8) :: bs(this%nrows), xs(this%nrows)

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

end module fdme_ams_precon_type
