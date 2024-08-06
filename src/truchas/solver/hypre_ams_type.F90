!!
!! HYPRE_AMS_TYPE
!!
!! This module defines an auxiliary-space Maxwell solver class that is built
!! on the AMS solver from the Hypre library.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module hypre_ams_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fhypre
  use simpl_mesh_type
  use pcsr_matrix_type
  use parallel_communication, only: global_maxval
  implicit none
  private

  type, public :: hypre_ams
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    integer :: nrows, ilower, iupper
    integer, allocatable :: rows(:)
    type(hypre_obj) :: solver = hypre_null_obj ! HYPRE_Solver object handle
    ! HYPRE_IJVector object handles
    type(hypre_obj) :: b = hypre_null_obj
    type(hypre_obj) :: x = hypre_null_obj
    type(hypre_obj) :: xn = hypre_null_obj
    type(hypre_obj) :: yn = hypre_null_obj
    type(hypre_obj) :: zn = hypre_null_obj
    type(hypre_obj) :: interior_nodes = hypre_null_obj
    ! HYPRE_IJMatrix object handles
    type(hypre_obj) :: A = hypre_null_obj
    type(hypre_obj) :: grad = hypre_null_obj
    real(r8), public :: rel_rnorm
    integer, public :: num_iter
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    final :: hypre_ams_final
  end type

contains

  subroutine hypre_ams_final(this)
    type(hypre_ams), intent(inout) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%b)) call fHYPRE_IJVectorDestroy(this%b, ierr)
    if (hypre_associated(this%x)) call fHYPRE_IJVectorDestroy(this%x, ierr)
    if (hypre_associated(this%xn)) call fHYPRE_IJVectorDestroy(this%xn, ierr)
    if (hypre_associated(this%yn)) call fHYPRE_IJVectorDestroy(this%yn, ierr)
    if (hypre_associated(this%zn)) call fHYPRE_IJVectorDestroy(this%zn, ierr)
    if (hypre_associated(this%interior_nodes)) call fHYPRE_IJVectorDestroy(this%interior_nodes, ierr)
    if (hypre_associated(this%A)) call fHYPRE_IJMatrixDestroy(this%A, ierr)
    if (hypre_associated(this%grad)) call fHYPRE_IJMatrixDestroy(this%grad, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_AMSDestroy(this%solver, ierr)
    INSIST(ierr == 0)
    call fHYPRE_ClearAllErrors
  end subroutine

  subroutine init(this, mesh, params, stat, errmsg)

    use parameter_list_type

    class(hypre_ams), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, ierr, ipar
    real(r8) :: rpar

    this%mesh => mesh

    this%nrows = this%mesh%edge_imap%onp_size
    this%ilower = this%mesh%edge_imap%first_gid
    this%iupper = this%mesh%edge_imap%last_gid
    this%rows = [(i, i=this%ilower, this%iupper)] ! global row indices for this process

    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%b, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%b, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate(this%ilower, this%iupper, this%x, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(this%x, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_AMSCreate(this%solver, ierr)
    INSIST(ierr == 0)

    call create_grad_matrix(this%mesh, this%grad)
    call fHYPRE_AMSSetDiscreteGradient(this%solver, this%grad, ierr)
    INSIST(ierr == 0)

    call create_coord_vectors(this%mesh, this%xn, this%yn, this%zn)
    call fHYPRE_AMSSetCoordinateVectors(this%solver, this%xn, this%yn, this%zn, ierr)
    INSIST(ierr == 0)

    !! Solver parameters
    call params%get('max-iter', ipar, stat, errmsg, default=100)
    if (stat /= 0) return
    call fHYPRE_AMSSetMaxIter(this%solver, ipar, ierr)
    INSIST(ierr == 0)

    call params%get('rel-tol', rpar, stat, errmsg, default=1.0e-8_r8)
    if (stat /= 0) return
    call fHYPRE_AMSSetTol(this%solver, rpar, ierr)
    INSIST(ierr == 0)

    call params%get('print-level', ipar, stat, errmsg, default=0)
    if (stat /= 0) return
    call fHYPRE_AMSSetPrintLevel(this%solver, ipar, ierr)
    INSIST(ierr == 0)

    call params%get('ams-cycle-type', ipar, stat, errmsg, default=1)
    if (stat /= 0) return
    call fHYPRE_AMSSetCycleType(this%solver, ipar, ierr)
    INSIST(ierr == 0)

    stat = 0

  end subroutine init


  subroutine setup(this, A, alpha, beta)

    class(hypre_ams), intent(inout) :: this
    type(pcsr_matrix), intent(in) :: A
    real(r8), intent(in) :: alpha(:), beta(:)

    integer :: ierr

    ASSERT(size(alpha) == this%mesh%ncell)
    ASSERT(size(beta)  == this%mesh%ncell)

    call A%copy_to_ijmatrix(this%A)

    if (global_maxval(beta) == 0.0_r8) then
      call fHYPRE_AMSSetBetaPoissonMatrix(this%solver, hypre_null_obj, ierr)
      INSIST(ierr == 0)
    else
      block
        integer :: nvalues, jlower, jupper, j
        integer, allocatable :: indices(:)
        real(r8), allocatable :: interior_nodes(:)

        allocate(interior_nodes(this%mesh%nnode), source=1.0_r8)
        do j = 1, this%mesh%ncell
          if (beta(j) > 0) interior_nodes(this%mesh%cnode(:,j)) = 0
        end do

        nvalues = this%mesh%nnode_onP
        jlower = this%mesh%node_imap%first_gid
        jupper = this%mesh%node_imap%last_gid
        indices = [(j, j=jlower, jupper)]

        call fHYPRE_IJVectorCreate(jlower, jupper, this%interior_nodes, ierr)
        call fHYPRE_IJVectorSetMaxOffProcElmts(this%interior_nodes, 0, ierr)
        call fHYPRE_IJVectorInitialize(this%interior_nodes, ierr)
        call fHYPRE_IJVectorSetValues(this%interior_nodes, nvalues, indices, interior_nodes, ierr)
        call fHYPRE_IJVectorAssemble(this%interior_nodes, ierr)
        INSIST(ierr == 0)
        call fHYPRE_AMSSetInteriorNodes(this%solver, this%interior_nodes, ierr)
        INSIST(ierr == 0)
      end block
    end if

    call fHYPRE_AMSSetup(this%solver, this%A, this%b, this%x, ierr)
    INSIST(ierr == 0)

  end subroutine setup


  subroutine solve(this, b, x, stat)

    class(hypre_ams), intent(inout) :: this
    real(r8), intent(in) :: b(:)
    real(r8), intent(inout) :: x(:)
    integer, intent(out) :: stat

    integer :: ierr

    ASSERT(size(b) >= this%mesh%nedge_onP)
    ASSERT(size(x) >= this%mesh%nedge_onP)

    call fHYPRE_ClearAllErrors

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize(this%b, ierr)
    call fHYPRE_IJVectorSetValues(this%b, this%nrows, this%rows, b, ierr)
    call fHYPRE_IJVectorAssemble(this%b, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    call fHYPRE_IJVectorInitialize(this%x, ierr)
    call fHYPRE_IJVectorSetValues(this%x, this%nrows, this%rows, x, ierr)
    call fHYPRE_IJVectorAssemble(this%x, ierr)
    INSIST(ierr == 0)

    !! Solve the system.
    call fHYPRE_AMSSolve(this%solver, this%A, this%b, this%x, ierr)
    if (ierr /= 0) then
      INSIST(ierr == HYPRE_ERROR_CONV)
      stat = 1
      call fHYPRE_ClearAllErrors
    else
      stat = 0
    end if

    !! Retrieve the solution vector from Hypre.
    call fHYPRE_IJVectorGetValues(this%x, this%nrows, this%rows, x, ierr)
    INSIST(ierr == 0)

    call fHYPRE_AMSGetNumIterations(this%solver, this%num_iter, ierr)
    INSIST(ierr == 0)
    call fHYPRE_AMSGetFinalRelativeResidualNorm(this%solver, this%rel_rnorm, ierr)
    INSIST(ierr == 0)

  end subroutine solve


  subroutine create_grad_matrix(mesh, grad)

    type(simpl_mesh), intent(in), target :: mesh
    type(hypre_obj) :: grad

    integer :: j
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix) :: matrix

    allocate(g)
    call g%init(mesh%edge_imap, mesh%node_imap)
    do j = 1, mesh%nedge_onP
      call g%add_edge(j, mesh%enode(:,j))
    end do
    call g%add_complete

    call matrix%init(g, take_graph=.true.)
    do j = 1, mesh%nedge_onP
      call matrix%set(j, mesh%enode(1,j), -1.0_r8)
      call matrix%set(j, mesh%enode(2,j), +1.0_r8)
    end do
    call matrix%copy_to_ijmatrix(grad)  !TODO: VERIFY PROCEDURE WORKS FOR NON-SQUARE MATRICES

  end subroutine create_grad_matrix


  subroutine create_coord_vectors(mesh, x, y, z)

    type(simpl_mesh), intent(in) :: mesh
    type(hypre_obj) :: x, y, z

    integer :: nvalues, jlower, jupper, j, ierr
    integer, allocatable :: indices(:)

    nvalues = mesh%nnode_onP
    jlower = mesh%node_imap%first_gid
    jupper = mesh%node_imap%last_gid
    indices = [(j, j=jlower, jupper)] ! global indices for this process

    call fHYPRE_IJVectorCreate(jlower, jupper, x, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(x, 0, ierr)
    call fHYPRE_IJVectorInitialize(x, ierr)
    call fHYPRE_IJVectorSetValues(x, nvalues, indices, mesh%x(1,:mesh%nnode_onP), ierr)
    call fHYPRE_IJVectorAssemble(x, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate(jlower, jupper, y, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(y, 0, ierr)
    call fHYPRE_IJVectorInitialize(y, ierr)
    call fHYPRE_IJVectorSetValues(y, nvalues, indices, mesh%x(2,:mesh%nnode_onP), ierr)
    call fHYPRE_IJVectorAssemble(x, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate(jlower, jupper, z, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts(z, 0, ierr)
    call fHYPRE_IJVectorInitialize(z, ierr)
    call fHYPRE_IJVectorSetValues(z, nvalues, indices, mesh%x(3,:mesh%nnode_onP), ierr)
    call fHYPRE_IJVectorAssemble(x, ierr)
    INSIST(ierr == 0)

  end subroutine create_coord_vectors

end module hypre_ams_type
