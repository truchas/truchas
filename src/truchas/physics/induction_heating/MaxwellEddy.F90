!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

#ifdef PGI_COMPILER_WORKAROUND
#define PGI_WORKAROUND
#endif

module MaxwellEddy

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use simpl_mesh_type
  use mimetic_discretization
  use state_history_type
  use sparse_matrix
  use CGSolver
  use truchas_logging_services

use debug_EM

  implicit none
  private

  type, public :: system
    type(simpl_mesh), pointer :: mesh => null() ! spatial discretization

    !! Solution history
    real(r8) :: t, dt
    type(state_history) :: ehist, bhist
    real(r8), allocatable :: g0(:)

    !! CG control descriptor
    type(cg_desc) :: cg = cg_desc(0, 1000, 1, 0.0_r8, 1.0e-6_r8)

    real(r8), allocatable :: eps(:)
    real(r8), allocatable :: mu(:)
    real(r8), allocatable :: sigma(:)
    real(r8) :: etasq, delta

    integer, allocatable :: emask(:)

    real(r8), allocatable :: mtr1(:,:)
    real(r8), allocatable :: mtr2(:,:)
    real(r8), allocatable :: mtr3(:,:,:)

    !! Hiptmair relaxation data
    logical, allocatable :: w0mask(:)
    logical, allocatable :: w1mask(:)
    type(msr_matrix) :: a0, a1
  contains
    procedure :: init => initialize_system
    procedure :: set_initial_state
    procedure :: set_param
    procedure :: step
    !procedure :: interpolate_fields ! currently unused
    procedure :: joule_heat
  end type system

  !! Look-aside data for CG-called procedures
  type(system), pointer :: cg_sys => null()

contains

  subroutine set_initial_state(this, t, e, b)
    class(system), intent(inout) :: this
    real(r8), intent(in) :: t, e(:), b(:)
    ASSERT( size(e) == this%mesh%nedge )
    ASSERT( size(b) == this%mesh%nface )
    this%t = t
    call this%ehist%init(3, t, e)
    call this%bhist%init(3, t, b)  ! really just mvec=1 if no interpolated output.
  end subroutine

  subroutine set_param(this, minitr, maxitr, tol, red, output_level)
    class(system), intent(inout) :: this
    real(r8), intent(in), optional :: tol, red
    integer, intent(in), optional :: minitr, maxitr, output_level
    if (present(minitr)) this%cg%minitr = minitr
    if (present(maxitr)) this%cg%maxitr = maxitr
    if (present(tol))    this%cg%tol = tol
    if (present(red))    this%cg%red = red
    if (present(output_level)) this%cg%output_level = output_level
  end subroutine

  subroutine interpolate_fields(this, t, e, b)
    class(system), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: e(:), b(:)
    call this%ehist%interp_state(t, e)
    call this%bhist%interp_state(t, b)
  end subroutine


  subroutine step(this, t, e, b, status, set_bv, bndry_src)

    class(system), intent(inout), target :: this
    real(r8), intent(out) :: t, e(:), b(:)
    integer, intent(out) :: status
    optional :: bndry_src

    interface
      subroutine set_bv (t, e)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in) :: t
        real(r8), intent(inout) :: e(:)
      end subroutine set_bv
      subroutine bndry_src (t, s)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in) :: t
        real(r8), intent(out) :: s(:)
      end subroutine bndry_src
    end interface

    real(r8) :: r(this%mesh%nedge), de(this%mesh%nedge), g(this%mesh%nedge)
    real(r8), pointer :: e0(:), b0(:)
    character(len=80) :: message

    ASSERT( size(e) == this%mesh%nedge )
    ASSERT( size(b) == this%mesh%nface )

    cg_sys => this
    t = this%t + this%dt

    call this%ehist%get_last_state_view(e0)
    call this%bhist%get_last_state_view(b0)

    if (present(bndry_src)) then        ! We have a boundary source term:
      if (.not.allocated(this%g0)) then
        allocate(this%g0(this%mesh%nedge))   !   allocate the component vector
        call bndry_src (this%t, this%g0)  !   compute the previous boundary source vector
      end if
      call bndry_src (t, g)
    else if (allocated(this%g0)) then    ! deallocate the component vector, so we know that
      deallocate(this%g0)                 ! the 'previous' value will need to be computed
    end if

    call this%ehist%interp_state(t, e)  ! Initial guess: quadratic extrapolation
    call set_bv (t, e)                  ! Set boundary values in the initial guess
    call this%mesh%edge_imap%gather_offp(e)
    call tr_res (this, e0, b0, e, r)     ! Compute the residual

    if (present(bndry_src)) r = r + (0.5_r8*this%dt)*(this%g0 + g)
    call this%mesh%edge_imap%gather_offp(r)
    !call verify_vector (r, this%mesh%edge_gs, 'R:')

    de = 0.0_r8 ! initial guess for the correction
    call SolveCG (this%cg, cg_ax, cg_pc, r, de, this%mesh%nedge_onP, status)!, this%mesh%edge_gs)
    call this%mesh%edge_imap%gather_offp(de)
    !call verify_vector (de, this%mesh%edge_gs, 'DE:')
    !stop
    if (status /= 0) return
    e = e + de
    !call verify_vector (e, this%mesh%edge_gs, 'E:')

    if (this%cg%output_level >= 3) then
      write(message,fmt='(t6,a,es10.3)') 'step |de|_max=', global_maxval(abs(de))
      call TLS_info (trim(message))
    end if

    !! Advance the B-field
    b = b0 - (0.5_r8 * this%dt * this%delta**2) * curl(this%mesh, e0 + e)
    !call verify_vector (b, this%mesh%face_gs, 'B:')

    !! Update the derived data type
    this%t = t
    call this%ehist%record_state(t, e)
    call this%bhist%record_state(t, b)
    this%g0 = g

  end subroutine step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  SPECIAL PURPOSE PROCEDURES FOR TRAPEZOID-RULE TIME INTEGRATION
!!

  subroutine initialize_system (this, mesh, eps, mu, sigma, etasq, delta, dt, emask)

    class(system), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: eps(:), mu(:), sigma(:)
    real(r8), intent(in) :: dt, etasq, delta
    integer, intent(in) :: emask(:)

    integer :: i, j, k, l, m, n
    real(r8) :: c(4,6), ctm2(6,4), ctm2c(6,6), g(6,4), a0(4,4), m1(21), m2(10)
    integer, allocatable :: ebedge(:)

    ASSERT( size(eps) == mesh%ncell )
    ASSERT( size(mu) == mesh%ncell )
    ASSERT( size(sigma) == mesh%ncell )
    ASSERT( size(emask) == mesh%nedge )

    this%mesh => mesh
    this%dt = dt
    this%etasq = etasq
    this%delta = delta

    allocate(this%eps(mesh%ncell), this%mu(mesh%ncell), this%sigma(mesh%ncell), this%emask(mesh%nedge))
    this%eps = eps
    this%mu = mu
    this%sigma = sigma
    this%emask = emask

    allocate(this%mtr1(21,mesh%ncell), this%mtr2(21,mesh%ncell), this%mtr3(6,4,mesh%ncell))

    c = reshape(source=(/ 0.0_r8,  0.0_r8,  1.0_r8,  1.0_r8, &
                          0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
                          0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
                          1.0_r8,  0.0_r8,  0.0_r8,  1.0_r8, &
                         -1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, &
                          1.0_r8,  1.0_r8,  0.0_r8,  0.0_r8 /), shape=shape(c))

    do j = 1, mesh%ncell
      !! Transpose(Curl).M2
      !ctms = matmul(Transpose(c), m2(:,:))
      m2 = W2_matrix_WE(mesh, j)
      do k = 1, 6
        ctm2(k,1) = m2(1)*c(1,k) + m2(2)*c(2,k) + m2(4)*c(3,k) + m2(7)*c(4,k)
        ctm2(k,2) = m2(2)*c(1,k) + m2(3)*c(2,k) + m2(5)*c(3,k) + m2(8)*c(4,k)
        ctm2(k,3) = m2(4)*c(1,k) + m2(5)*c(2,k) + m2(6)*c(3,k) + m2(9)*c(4,k)
        ctm2(k,4) = m2(7)*c(1,k) + m2(8)*c(2,k) + m2(9)*c(3,k) + m2(10)*c(4,k)
      end do

      !! h Transpose(Curl).M2
      this%mtr3(:,:,j) = (dt/mu(j)) * ctm2

      !! (h/2)^2 Transpose(Curl).MS.Curl
      ctm2c = ((0.5_r8*dt*delta)**2/mu(j)) * matmul(ctm2, c)

      m1 = W1_matrix_WE(mesh, j)
      l = 0
      do k = 1, 6     ! NB: we traverse the elements of the upper triangular part
        do i = 1, k   ! column by column.
          l = l + 1
          this%mtr1(l,j) = (etasq*eps(j) + 0.5_r8*dt*sigma(j)) * m1(l) + ctm2c(i,k)
          this%mtr2(l,j) = (etasq*eps(j) - 0.5_r8*dt*sigma(j)) * m1(l) - ctm2c(i,k)
        end do
      end do

    end do

    !! Create the mask arrays for the spaces: any edge not tagged as a
    !! dirichlet (boundary) edge is masked (w1mask); any node on a unmasked
    !! edge is unmasked (w0mask).  'Masked' means a true value and marks
    !! a DOF that belongs to the space.

    allocate(this%w0mask(mesh%nnode), this%w1mask(mesh%nedge))
    this%w1mask = (emask == 0)
    this%w0mask = .true.
    do j = 1, mesh%nedge
      if (.not.this%w1mask(j)) this%w0mask(mesh%enode(:,j)) = .false.
    end do
    call this%mesh%node_imap%gather_offp(this%w0mask)

    !! Assemble the edge-based coefficient matrix A1
    call this%a1%init(mesh%nedge, mesh%cedge)
    do j = 1, mesh%ncell
      l = 0
      do n = 1, 6     ! NB: we traverse the elements (m,n) of the upper triangular
        do m = 1, n   ! part column by column
          l = l + 1
          call this%a1%add_to(mesh%cedge(m,j), mesh%cedge(n,j), this%mtr1(l,j))
          if (m == n) cycle
          call this%a1%add_to(mesh%cedge(n,j), mesh%cedge(m,j), this%mtr1(l,j))
        end do
      end do
    end do

    !! Project out the rows and columns corresponding to Dirichlet edges.
    allocate(ebedge(count(emask/=0)))
    ebedge = pack(array=(/(j,j=1,mesh%nedge)/), mask=(emask/=0))
    call remove_dof (this%a1, ebedge)
    deallocate(ebedge)
    !ASSERT( is_symmetric(this%a1) )

    !! Projected System.  We form the projected system defined on the nullspace
    !! of the curl operator.  This corresponds to the range of the gradient
    !! operator, and thus the system is representable as node-based system

    !! Local gradient matrix
    g = reshape(source=(/ -1.0_r8, -1.0_r8, -1.0_r8,  0.0_r8,  0.0_r8,  0.0_r8, &
                           1.0_r8,  0.0_r8,  0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
                           0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
                           0.0_r8,  0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  1.0_r8/), shape=shape(g))

    !! Compute and assemble the node-based projected matrix A0.
    call this%a0%init(mesh%nnode, mesh%cnode)
    do j = 1, this%mesh%ncell
      !if (this%fs_cell(j)) cycle
      a0 = (etasq*eps(j) + 0.5_r8*dt*sigma(j)) * matmul(transpose(g), sym_matmul(W1_matrix_WE(mesh, j), g))
      do k = 1, 4
        do i = 1, 4
          call this%a0%add_to(mesh%cnode(i,j), mesh%cnode(k,j), a0(i,k))
        end do
      end do
    end do

    !! Project out the rows and columns corresponding to Dirichlet edges.
    allocate(ebedge(count(.not.this%w0mask)))
    ebedge = pack(array=(/(j,j=1,mesh%nnode)/), mask=.not.this%w0mask)
    call remove_dof (this%a0, ebedge)
    deallocate(ebedge)
    !ASSERT( is_symmetric(this%a0) )

  end subroutine initialize_system

  subroutine remove_dof (this, dof)

    type(msr_matrix), intent(inout) :: this
    integer, intent(in) :: dof(:)

    integer :: j, m, n, lmn, lnm

    ASSERT( this%nrow == this%ncol )

    do j = 1, size(dof)
      m = dof(j)
      call this%project_out(m)
      this%diag(m) = 1.0_r8
    end do

  end subroutine remove_dof

  !!
  !! Residual for the time-discretized system (trapezoid-rule).
  !!

  subroutine tr_res (sys, e, b, e0, r)

    type(system),  intent(in)  :: sys
    real(r8), intent(in)  :: e(:), e0(:), b(:)
    real(r8), intent(out) :: r(:)

    integer :: i, j, k, l1, l2
    real(r8) :: value, el(6), e0l(6), bl(4)

    ASSERT( size(e) == sys%mesh%nedge )
    ASSERT( size(e0) == sys%mesh%nedge )
    ASSERT( size(r) == sys%mesh%nedge )
    ASSERT( size(b) == sys%mesh%nface )

    r = 0.0_r8
    do j = 1, sys%mesh%ncell
      do k = 1, 6   ! Gather local copies of the vectors
        i = sys%mesh%cedge(k,j)
        el(k)  = e(i)
        e0l(k) = e0(i)
      end do
      do k = 1, 4
        bl(k) = b(sys%mesh%cface(k,j))
      end do
      l1 = 1
      do k = 1, 6
        value = r(sys%mesh%cedge(k,j))
        l2 = l1
        do i = 1, k-1
          value = value + sys%mtr2(l2,j)*el(i) - sys%mtr1(l2,j)*e0l(i)
          l2 = l2 + 1
        end do
        do i = k, 6
          value = value + sys%mtr2(l2,j)*el(i) - sys%mtr1(l2,j)*e0l(i)
          l2 = l2 + i
        end do
        do i = 1, 4
          value = value + sys%mtr3(k,i,j) * bl(i)
        end do
        r(sys%mesh%cedge(k,j)) = value
        l1 = l1 + k
      end do
    end do

    where (sys%emask /= 0) r = 0.0_r8

  end subroutine tr_res

  !!
  !! Hiptmair relaxation procedure
  !!

  subroutine hiptmair (sys, f, u)

    type(system),  intent(in)  :: sys
    real(r8), intent(in)  :: f(:)
    real(r8), intent(out) :: u(:)

    integer :: nedge_onP, nnode_onP
    real(r8) :: un(sys%mesh%nnode), rn(sys%mesh%nnode), r(sys%mesh%nedge)

    nedge_onP = sys%mesh%nedge_onP
    nnode_onP = sys%mesh%nnode_onP

    ASSERT( all(sys%w1mask .or. (f == 0.0_r8)) )

    !! Forward Gauss-Seidel relaxation on the on-process edge system.
    u = 0.0_r8
    call gs_relaxation (sys%a1, f(:nedge_onP), u, pattern='f')

    !! Update the local residual and project it to the nodes.
    r = f - sys%a1%matvec(u)
    rn = grad_t(sys%mesh, r)
    where (.not.sys%w0mask) rn = 0.0_r8

    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
    un = 0.0_r8
    call gs_relaxation (sys%a0, rn(:nnode_onP), un, pattern='fb')

    !! Update the the solution with the node-based correction.
    u = u + grad(sys%mesh, un)

    !! Backward Gauss-Seidel relaxation on the on-process edge system.
    call gs_relaxation (sys%a1, f(:nedge_onP), u, pattern='b')

    call sys%mesh%edge_imap%scatter_offp_sum(u)
    call sys%mesh%edge_imap%gather_offp(u)

    ASSERT( all(sys%w1mask .or. (u == 0.0_r8)) )

  end subroutine hiptmair

  !!
  !! These procedures get passed to the CG solver to compute the matrix-vector
  !! product and preconditioning.  Note that they reference the (previously
  !! initialized) module variable CG_SYS which points to the additional data
  !! required by the computations.
  !!

  subroutine cg_ax (x, y)

    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)

    ASSERT( all(cg_sys%w1mask .or. (x == 0.0_r8)) )

    y = cg_sys%a1%matvec(x)
    ASSERT( all(cg_sys%w1mask .or. (y == 0.0_r8)) )

    call cg_sys%mesh%edge_imap%gather_offp(y)

  end subroutine cg_ax


  subroutine cg_pc (x, y)
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    call hiptmair (cg_sys, x, y)
  end subroutine cg_pc


  function joule_heat(this, u) result(q)

    class(system), intent(in) :: this
    real(r8), intent(in) :: u(:)
    real(r8) :: q(this%mesh%ncell)

    integer :: i, j, k, l1, l2
    real(r8) :: value, ul(6), m1(21)

    ASSERT( size(u) == this%mesh%nedge )

    q = 0.0_r8
    do j = 1, this%mesh%ncell
      if (this%sigma(j) == 0.0_r8) cycle
      do k = 1, 6  ! gather local copy of the vector
        ul(k) = u(this%mesh%cedge(k,j))
      end do
      m1 = W1_matrix_WE(this%mesh, j)
      l1 = 1
      do k = 1, 6
        value = 0.0_r8
        l2 = l1
        do i = 1, k-1
          value = value + m1(l2) * ul(i)
          l2 = l2 + 1
        end do
        do i = k, 6
          value = value + m1(l2) * ul(i)
          l2 = l2 + i
        end do
        q(j) = q(j) + ul(k) * value
        l1 = l1 + k
      end do
      q(j) = this%sigma(j) * q(j) / abs(this%mesh%volume(j))  ! want a source density
    end do

  end function joule_heat

  subroutine hiptmair_is_symmetric (sys)

    type(system), intent(in), target :: sys

    integer :: j, jproc, k, kproc, kstart, lun, nedge_onP, n, m, bsize_vector(nPE)
    real(r8), dimension(sys%mesh%nedge) :: ej, ek, tmp
    real(r8) :: pe(sys%mesh%nedge,sys%mesh%edge_imap%global_size)
    character(len=1) :: edge_type(sys%mesh%nedge), etype(sys%mesh%edge_imap%global_size)
    real(r8) :: ej_dot_pek, ek_dot_pej, ej_dot_pej

    nedge_onP = sys%mesh%nedge_onP
    cg_sys => sys

    call gather (sys%mesh%edge_imap%onp_size, bsize_vector)

    !! Mark the on-PE edges as partition boundary, partition interior, or ignored
    tmp(1:nedge_onP)  = 0.0_r8
    tmp(nedge_onP+1:) = 1.0_r8
    call sys%mesh%edge_imap%scatter_offp_max(tmp)
    where (tmp > 0.0_r8)
      edge_type = 'B'
    elsewhere
      edge_type = 'I'
    end where
    where (sys%emask /= 0) edge_type = '*'
    call gather (edge_type(:nedge_onP), etype)

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
        call sys%mesh%edge_imap%gather_offp(ej)

        call hiptmair (sys, ej, pe(:,n))
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
        call sys%mesh%edge_imap%gather_offp(ej)

        m = n
        do kproc = jproc, nPE
          kstart = 1
          if (kproc == jproc) kstart = j + 1
          do k = kstart, bsize_vector(kproc)
            m = m + 1

            !! Initialize the (kproc, k) unit vector
            ek = 0.0_r8
            if (this_PE == kproc) ek(k) = 1.0_r8
            call sys%mesh%edge_imap%gather_offp(ek)

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

    close(unit=lun)

  contains

    subroutine broadcast_from_PE (value, pe)
      character(len=*), intent(inout) :: value
      integer, intent(in) :: pe
      character(len(value)) :: tmp(nPE)
      call gather ([value], tmp)
      if (is_IOP) value = tmp(pe)
      call broadcast (value)
    end subroutine broadcast_from_PE

  end subroutine hiptmair_is_symmetric

end module MaxwellEddy
