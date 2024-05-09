!!
!! TDME_CG_SOLVER_TYPE
!!
!! This module provides a derived type that defines a preconditioned CG solver
!! for the linear time step system of the discrete time-domain Maxwell equations.
!! The preconditioner uses Hiptmair smoothing on an auxiliary space.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Refactored February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module tdme_cg_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cg_solver_class
  use simpl_mesh_type
  use tdme_model_type
  use msr_matrix_type
  use mimetic_discretization
  implicit none
  private

  type, extends(cg_solver), public :: tdme_cg_solver
    type(simpl_mesh), pointer :: mesh => null()   ! unowned reference
    type(tdme_model), pointer :: model => null()  ! unowned reference
    !! Hiptmair preconditioner data
    logical, allocatable :: w0mask(:)
    type(msr_matrix) :: a0, a1
    real(r8), allocatable :: un(:), rn(:), r(:) ! persistent local workspace for pc
  contains
    procedure :: init
    procedure :: ax
    procedure :: pc
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    use upper_packed_matrix_procs, only: upm_cong_prod

    class(tdme_cg_solver), intent(out) :: this
    type(tdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    !! Local gradient matrix
    real(r8), parameter :: grad(6,4) = reshape([-1, -1, -1,  0,  0,  0, &
                                                 1,  0,  0, -1, -1,  0, &
                                                 0,  1,  0,  1,  0, -1, &
                                                 0,  0,  1,  0,  1,  1], shape(grad))

    integer :: j
    real(r8) :: a0(10)
    integer, allocatable :: ebedge(:)

    this%model => model
    this%mesh => model%mesh

    this%minitr = 1
    call params%get('cg-max-iter', this%maxitr, stat, errmsg, default=500)
    if (stat /= 0) return
    if (this%maxitr < 1) then
      stat = 1
      errmsg = 'cg-max-iter must be > 0'
      return
    end if
    this%tol = 0.0_r8 ! we do not use an absolute residual tolerance
    call params%get('cg-tol', this%red, stat, errmsg, default=1.0e-8_r8)
    if (stat /= 0) return
    if (this%red <= 0.0_r8 .or. this%red >= 0.1_r8) then
      stat = 1
      errmsg = 'cg-tol must be > 0.0 and < 0.1'
      return
    end if
    call params%get('output-level', this%output_level, stat, errmsg, default=1)
    if (stat /= 0) return

    !! Create the node space mask array: true values correspond to nodes not
    !! belonging to any E-field Dirichlet edge. TODO: replace by node list
    if (allocated(this%model%ebc)) then
      allocate(this%w0mask(this%mesh%nnode))
      this%w0mask = .true.
      do j = 1, size(this%model%ebc%index)
        this%w0mask(this%mesh%enode(:,this%model%ebc%index(j))) = .false.
      end do
      call this%mesh%node_imap%gather_offp(this%w0mask)
    end if

    !! Assemble the edge-based coefficient matrix A1
    block
      type(msr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nedge)
      call g%add_clique(this%mesh%cedge)
      call g%add_complete
      call this%a1%init(g, take_graph=.true.)
    end block

    do j = 1, this%mesh%ncell
      call this%a1%add_to(this%mesh%cedge(:,j), this%model%mtr1(:,j))
    end do

    !! Project out the rows and columns corresponding to Dirichlet edges.
    if (allocated(this%model%ebc)) then
      associate (index => this%model%ebc%index)
        do j = 1, size(index)
          call this%a1%project_out(index(j))
          call this%a1%set(index(j), index(j), 1.0_r8)
        end do
      end associate
    end if

    !! Projected System. We form the projected system defined on the nullspace
    !! of the curl operator. This corresponds to the range of the gradient
    !! operator, and thus the system is representable as node-based system

    !! Compute and assemble the node-based projected matrix A0.
    block
      type(msr_graph), pointer :: g
      allocate(g)
      call g%init(this%mesh%nnode)
      call g%add_clique(this%mesh%cnode)
      call g%add_complete
      call this%a0%init(g, take_graph=.true.)
    end block

    associate (eps => this%model%eps, sigma => this%model%sigma, dt => this%model%dt)
      do j = 1, this%mesh%ncell
        a0 = (eps(j) + 0.5_r8*dt*sigma(j)) * upm_cong_prod(6, 4, W1_matrix_WE(this%mesh, j), grad)
        call this%a0%add_to(this%mesh%cnode(:,j), a0)
      end do
    end associate

    !! Project out the rows and columns corresponding to Dirichlet edges.
    if (allocated(this%w0mask)) then
      allocate(ebedge(count(.not.this%w0mask)))
      ebedge = pack(array=[(j,j=1,this%mesh%nnode)], mask=.not.this%w0mask)
      do j = 1, size(ebedge)
        call this%a0%project_out(ebedge(j))
        call this%a0%set(ebedge(j), ebedge(j), 1.0_r8)
      end do
      deallocate(ebedge)
    end if
    !ASSERT(this%a0%is_symmetric()) ! will generally fail due to order-of-operation differences

    allocate(this%un(this%mesh%nnode), this%rn(this%mesh%nnode), this%r(this%mesh%nedge))

  end subroutine init

  subroutine ax(this, x, y)
    class(tdme_cg_solver), intent(inout) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    !ASSERT(all(x(this%model%ebc%index) == 0.0_r8))
    y = this%a1%matvec(x)
    !ASSERT(all(y(this%model%ebc%index) == 0.0_r8))
    call this%mesh%edge_imap%gather_offp(y)
  end subroutine

  subroutine pc(this, x, y)
    class(tdme_cg_solver), intent(inout) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    call hiptmair(this, x, y)
  end subroutine

  subroutine hiptmair(this, f, u)

    use mimetic_discretization

    type(tdme_cg_solver), intent(inout) :: this ! intent(in) except for workspace
    real(r8), intent(in)  :: f(:)
    real(r8), intent(out) :: u(:)

    integer :: nedge_onP, nnode_onP

    nedge_onP = this%mesh%nedge_onP
    nnode_onP = this%mesh%nnode_onP

    !ASSERT(all(f(this%model%ebc%index) == 0.0_r8))

    !! Forward Gauss-Seidel relaxation on the on-process edge system.
    u = 0.0_r8
    call gs_relaxation(this%a1, f(:nedge_onP), u, pattern='f')

    !! Update the local residual and project it to the nodes.
    this%r(:) = f - this%a1%matvec(u)
    call grad_t(this%mesh, this%r, this%rn)
    if (allocated(this%w0mask)) where (.not.this%w0mask) this%rn = 0.0_r8

    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
    this%un = 0.0_r8
    call gs_relaxation(this%a0, this%rn(:nnode_onP), this%un, pattern='fb')

    !! Update the the solution with the node-based correction.
    call grad(this%mesh, this%un, u, increment=.true.)

    !! Backward Gauss-Seidel relaxation on the on-process edge system.
    call gs_relaxation(this%a1, f(:nedge_onP), u, pattern='b')

    call this%mesh%edge_imap%scatter_offp_sum(u)
    call this%mesh%edge_imap%gather_offp(u)

    !ASSERT(all(u(this%model%ebc%index) == 0.0_r8))

  end subroutine hiptmair

  !! This debugging procedure is retained for reference. When Joule heating was
  !! first introduced, non-symmetry of the Hiptmair preconditioning in parallel
  !! resulted in CG solver failures. This procedure was key in confirming the
  !! cause of the failures and ultimately fixing the parallel symmetry issue.

  subroutine hiptmair_is_symmetric(this)

    use parallel_communication

    type(tdme_cg_solver), intent(inout), target :: this

    integer :: j, jproc, k, kproc, kstart, lun, nedge_onP, n, m, bsize_vector(nPE)
    real(r8), dimension(this%mesh%nedge) :: ej, ek, tmp
    real(r8) :: pe(this%mesh%nedge,this%mesh%edge_imap%global_size)
    character(len=1) :: edge_type(this%mesh%nedge), etype(this%mesh%edge_imap%global_size)
    real(r8) :: ej_dot_pek, ek_dot_pej, ej_dot_pej

    nedge_onP = this%mesh%nedge_onP

    call gather(this%mesh%edge_imap%onp_size, bsize_vector)

    !! Mark the on-PE edges as partition boundary, partition interior, or ignored
    tmp(1:nedge_onP)  = 0.0_r8
    tmp(nedge_onP+1:) = 1.0_r8
    call this%mesh%edge_imap%scatter_offp_max(tmp)
    where (tmp > 0.0_r8)
      edge_type = 'B'
    elsewhere
      edge_type = 'I'
    end where
    if (allocated(this%model%ebc)) then
      do j = 1, size(this%model%ebc%index)
        edge_type(this%model%ebc%index(j)) = '*'
      end do
    end if
    call gather(edge_type(:nedge_onP), etype)

    if (is_IOP) then
      open(newunit=lun,file='hiptmair.dat',status='replace',action='write')
      write(unit=lun,fmt=*) nPE
      write(unit=lun,fmt=*) bsize_vector
    end if

    if (is_IOP) print *, 'hiptmair_is_symmetric: computing matvecs...'
    n = 0
    do jproc = 1, nPE
      do j = 1, bsize_vector(jproc)
        n = n + 1

        !! Initialize the (jproc, j) unit vector
        ej = 0.0_r8
        if (this_PE == jproc) ej(j) = 1.0_r8
        call this%mesh%edge_imap%gather_offp(ej)

        call hiptmair(this, ej, pe(:,n))
        ej_dot_pej = global_dot_product(ej(:nedge_onP), pe(:nedge_onP,n))
        if (is_IOP) then
          write(unit=lun,fmt='(i2,1x,i4,1x,a1,1x,e24.16)') jproc, j, etype(n), ej_dot_pej
        end if

      end do
    end do
    if (is_IOP) print *, 'hiptmair_is_symmetric: matvecs done; computing dot products...'
    return
    n = 0
    do jproc = 1, nPE
      do j = 1, bsize_vector(jproc)
        n = n + 1

        !! Initialize the (jproc, j) unit vector
        ej = 0.0_r8
        if (this_PE == jproc) ej(j) = 1.0_r8
        call this%mesh%edge_imap%gather_offp(ej)

        m = n
        do kproc = jproc, nPE
          kstart = 1
          if (kproc == jproc) kstart = j + 1
          do k = kstart, bsize_vector(kproc)
            m = m + 1

            !! Initialize the (kproc, k) unit vector
            ek = 0.0_r8
            if (this_PE == kproc) ek(k) = 1.0_r8
            call this%mesh%edge_imap%gather_offp(ek)

            ej_dot_pek = global_dot_product(ej(:nedge_onP), pe(:nedge_onP,m))
            ek_dot_pej = global_dot_product(ek(:nedge_onP), pe(:nedge_onP,n))

            if (is_IOP) then
              write(unit=lun,fmt='(2(i2,1x,i4,1x,a1,1x),2e24.16)') &
                jproc, j, etype(n), kproc, k, etype(m), ej_dot_pek, ek_dot_pej
            end if

          end do
        end do
      end do
    end do

    close(lun)

  contains

    subroutine broadcast_from_PE(value, pe)
      character(*), intent(inout) :: value
      integer, intent(in) :: pe
      character(len(value)) :: tmp(nPE)
      call gather ([value], tmp)
      if (is_IOP) value = tmp(pe)
      call broadcast(value)
    end subroutine

  end subroutine hiptmair_is_symmetric

end module tdme_cg_solver_type
