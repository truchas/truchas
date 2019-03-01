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

  use kinds
  use parallel_communication
  use simpl_mesh_type
  use index_partitioning
  use mimetic_discretization
  use solution_history
  use sparse_matrix
  use CGSolver
  use truchas_logging_services
  
use debug_EM

  implicit none
  private
  
  public :: initialize_system, set_initial_values, set_param, step, interpolate_fields
  public :: joule_heat, destroy_system
  
  type, public :: system
    type(simpl_mesh), pointer :: mesh => null() ! spatial discretization
    
    !! Solution history
    real(kind=r8) :: t, dt
    type(history) :: ehist, bhist
    real(kind=r8), allocatable :: g0(:)
    
    !! CG control descriptor
    type(cg_desc) :: cg = cg_desc(0, 1000, 1, 0.0_r8, 1.0e-6_r8)
    
    real(kind=r8), allocatable :: eps(:)
    real(kind=r8), allocatable :: mu(:)
    real(kind=r8), allocatable :: sigma(:)
    real(kind=r8) :: etasq, delta
    
    integer, allocatable :: emask(:)
    
    real(kind=r8), allocatable :: mtr1(:,:)
    real(kind=r8), allocatable :: mtr2(:,:)
    real(kind=r8), allocatable :: mtr3(:,:,:)
    
    !! Hiptmair relaxation data
    logical, allocatable :: w0mask(:)
    logical, allocatable :: w1mask(:)
    type(msr_matrix) :: a0, a1
  end type system
  
  !! Look-aside data for CG-called procedures
  type(system), pointer :: cg_sys => null()
  
contains

  subroutine destroy_system (sys)
    type(system), intent(inout) :: sys
    type(system) :: default
    nullify(sys%mesh)
    if (allocated(sys%g0)) deallocate(sys%g0)
    if (allocated(sys%eps)) deallocate(sys%eps)
    if (allocated(sys%mu)) deallocate(sys%mu)
    if (allocated(sys%sigma)) deallocate(sys%sigma)
    if (allocated(sys%emask)) deallocate(sys%emask)
    if (allocated(sys%mtr1)) deallocate(sys%mtr1)
    if (allocated(sys%mtr2)) deallocate(sys%mtr2)
    if (allocated(sys%mtr3)) deallocate(sys%mtr3)
    if (allocated(sys%w0mask)) deallocate(sys%w0mask)
    if (allocated(sys%w1mask)) deallocate(sys%w1mask)
    call destroy_msr_matrix (sys%a0)
    call destroy_msr_matrix (sys%a1)
    call destroy (sys%ehist)
    call destroy (sys%bhist)
    sys = default ! set default initialization values
  end subroutine destroy_system

  subroutine set_initial_values (sys, t, e, b)
    type(system), intent(inout) :: sys
    real(kind=r8), intent(in) :: t, e(:), b(:)
    ASSERT( size(e) == sys%mesh%nedge )
    ASSERT( size(b) == sys%mesh%nface )
    sys%t = t
    call create_history (sys%ehist, 3, t, e)
    call create_history (sys%bhist, 3, t, b)  ! really just mvec=1 if no interpolated output.
  end subroutine set_initial_values
  
  subroutine set_param (sys, minitr, maxitr, tol, red, output_level)
    type(system), intent(inout) :: sys
    real(kind=r8), intent(in), optional :: tol, red
    integer, intent(in), optional :: minitr, maxitr, output_level
    if (present(minitr)) sys%cg%minitr = minitr
    if (present(maxitr)) sys%cg%maxitr = maxitr
    if (present(tol))    sys%cg%tol = tol
    if (present(red))    sys%cg%red = red
    if (present(output_level)) sys%cg%output_level = output_level
  end subroutine set_param
  
  subroutine interpolate_fields (sys, t, e, b)
    type(system),  intent(in)  :: sys
    real(kind=r8), intent(in)  :: t
    real(kind=r8), intent(out) :: e(:), b(:)
    call interpolate_solution (sys%ehist, t, e)
    call interpolate_solution (sys%bhist, t, b)
  end subroutine interpolate_fields
  

  subroutine step (sys, t, e, b, status, set_bv, bndry_src)
  
    type(system), intent(inout), target :: sys
    real(kind=r8), intent(out) :: t, e(:), b(:)
    integer, intent(out) :: status
    optional :: bndry_src

    interface
      subroutine set_bv (t, e)
        use kinds, only: r8
        real(kind=r8), intent(in) :: t
        real(kind=r8), intent(inout) :: e(:)
      end subroutine set_bv
      subroutine bndry_src (t, s)
        use kinds, only: r8
        real(kind=r8), intent(in) :: t
        real(kind=r8), intent(out) :: s(:)
      end subroutine bndry_src
    end interface
    
    real(kind=r8) :: r(sys%mesh%nedge), de(sys%mesh%nedge), g(sys%mesh%nedge)
    real(kind=r8), pointer :: e0(:), b0(:)
    character(len=80) :: message

    ASSERT( size(e) == sys%mesh%nedge )
    ASSERT( size(b) == sys%mesh%nface )
    
    cg_sys => sys
    t = sys%t + sys%dt
    
    e0 => most_recent_solution(sys%ehist)
    b0 => most_recent_solution(sys%bhist)
    
    if (present(bndry_src)) then        ! We have a boundary source term:
      if (.not.allocated(sys%g0)) then
        allocate(sys%g0(sys%mesh%nedge))   !   allocate the component vector
        call bndry_src (sys%t, sys%g0)  !   compute the previous boundary source vector
      end if
      call bndry_src (t, g)
    else if (allocated(sys%g0)) then    ! deallocate the component vector, so we know that
      deallocate(sys%g0)                 ! the 'previous' value will need to be computed
    end if
    
    call interpolate_solution (sys%ehist, t, e)  ! Initial guess: quadratic extrapolation
    call set_bv (t, e)                  ! Set boundary values in the initial guess
    call gather_boundary (sys%mesh%edge_ip, e)
    call tr_res (sys, e0, b0, e, r)     ! Compute the residual
    
    if (present(bndry_src)) r = r + (0.5_r8*sys%dt)*(sys%g0 + g)
    call gather_boundary (sys%mesh%edge_ip, r)
    !call verify_vector (r, sys%mesh%edge_gs, 'R:')

    de = 0.0_r8 ! initial guess for the correction
    call SolveCG (sys%cg, cg_ax, cg_pc, r, de, sys%mesh%nedge_onP, status)!, sys%mesh%edge_gs)
    call gather_boundary (sys%mesh%edge_ip, de)
    !call verify_vector (de, sys%mesh%edge_gs, 'DE:')
    !stop
    if (status /= 0) return
    e = e + de
    !call verify_vector (e, sys%mesh%edge_gs, 'E:')
    
    if (sys%cg%output_level >= 3) then
      write(message,fmt='(t6,a,es10.3)') 'step |de|_max=', global_maxval(abs(de))
      call TLS_info (trim(message))
    end if

    !! Advance the B-field
    b = b0 - (0.5_r8 * sys%dt * sys%delta**2) * curl(sys%mesh, e0 + e)
    !call verify_vector (b, sys%mesh%face_gs, 'B:')
    
    !! Update the derived data type
    sys%t = t
    call record_solution (sys%ehist, t, e)
    call record_solution (sys%bhist, t, b)
    sys%g0 = g
    
  end subroutine step
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  SPECIAL PURPOSE PROCEDURES FOR TRAPEZOID-RULE TIME INTEGRATION
!!

  subroutine initialize_system (sys, mesh, eps, mu, sigma, etasq, delta, dt, emask)
  
    type(system), intent(out) :: sys
    type(simpl_mesh), intent(in), target :: mesh
    real(kind=r8), intent(in) :: eps(:), mu(:), sigma(:)
    real(kind=r8), intent(in) :: dt, etasq, delta
    integer, intent(in) :: emask(:)
    
    integer :: i, j, k, l, m, n
    real(kind=r8) :: c(4,6), ctm2(6,4), ctm2c(6,6), g(6,4), a0(4,4), m1(21), m2(10)
    integer, allocatable :: ebedge(:)
    
    ASSERT( size(eps) == mesh%ncell )
    ASSERT( size(mu) == mesh%ncell )
    ASSERT( size(sigma) == mesh%ncell )
    ASSERT( size(emask) == mesh%nedge )
    
    sys%mesh => mesh
    sys%dt = dt
    sys%etasq = etasq
    sys%delta = delta
    
    allocate(sys%eps(mesh%ncell), sys%mu(mesh%ncell), sys%sigma(mesh%ncell), sys%emask(mesh%nedge))
    sys%eps = eps
    sys%mu = mu
    sys%sigma = sigma
    sys%emask = emask
    
    allocate(sys%mtr1(21,mesh%ncell), sys%mtr2(21,mesh%ncell), sys%mtr3(6,4,mesh%ncell))
    
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
      sys%mtr3(:,:,j) = (dt/mu(j)) * ctm2
      
      !! (h/2)^2 Transpose(Curl).MS.Curl
      ctm2c = ((0.5_r8*dt*delta)**2/mu(j)) * matmul(ctm2, c)
      
      m1 = W1_matrix_WE(mesh, j)
      l = 0
      do k = 1, 6     ! NB: we traverse the elements of the upper triangular part
        do i = 1, k   ! column by column.
          l = l + 1
          sys%mtr1(l,j) = (etasq*eps(j) + 0.5_r8*dt*sigma(j)) * m1(l) + ctm2c(i,k)
          sys%mtr2(l,j) = (etasq*eps(j) - 0.5_r8*dt*sigma(j)) * m1(l) - ctm2c(i,k)
        end do
      end do
      
    end do
    
    !! Create the mask arrays for the spaces: any edge not tagged as a
    !! dirichlet (boundary) edge is masked (w1mask); any node on a unmasked
    !! edge is unmasked (w0mask).  'Masked' means a true value and marks
    !! a DOF that belongs to the space.
    
    allocate(sys%w0mask(mesh%nnode), sys%w1mask(mesh%nedge))
    sys%w1mask = (emask == 0)
    sys%w0mask = .true.
    do j = 1, mesh%nedge
      if (.not.sys%w1mask(j)) sys%w0mask(mesh%enode(:,j)) = .false.
    end do
    call gather_boundary (sys%mesh%node_ip, sys%w0mask)
    
    !! Assemble the edge-based coefficient matrix A1
    call create_msr_matrix (mesh%nedge, mesh%cedge, sys%a1)
    do j = 1, mesh%ncell
      l = 0
      do n = 1, 6     ! NB: we traverse the elements (m,n) of the upper triangular
        do m = 1, n   ! part column by column
          l = l + 1
          call increment_msr_matrix_value (sys%a1, mesh%cedge(m,j), mesh%cedge(n,j), sys%mtr1(l,j))
          if (m == n) cycle
          call increment_msr_matrix_value (sys%a1, mesh%cedge(n,j), mesh%cedge(m,j), sys%mtr1(l,j))
        end do
      end do
    end do
    
    !! Project out the rows and columns corresponding to Dirichlet edges.
    allocate(ebedge(count(emask/=0)))
    ebedge = pack(array=(/(j,j=1,mesh%nedge)/), mask=(emask/=0))
    call remove_dof (sys%a1, ebedge)
    deallocate(ebedge)
    !ASSERT( is_symmetric(sys%a1) )

    !! Projected System.  We form the projected system defined on the nullspace
    !! of the curl operator.  This corresponds to the range of the gradient
    !! operator, and thus the system is representable as node-based system
    
    !! Local gradient matrix
    g = reshape(source=(/ -1.0_r8, -1.0_r8, -1.0_r8,  0.0_r8,  0.0_r8,  0.0_r8, &
                           1.0_r8,  0.0_r8,  0.0_r8, -1.0_r8, -1.0_r8,  0.0_r8, &
                           0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  0.0_r8, -1.0_r8, &
                           0.0_r8,  0.0_r8,  1.0_r8,  0.0_r8,  1.0_r8,  1.0_r8/), shape=shape(g))
                             
    !! Compute and assemble the node-based projected matrix A0.
    call create_msr_matrix (mesh%nnode, mesh%cnode, sys%a0)
    do j = 1, sys%mesh%ncell
      !if (sys%fs_cell(j)) cycle
      a0 = (etasq*eps(j) + 0.5_r8*dt*sigma(j)) * matmul(transpose(g), sym_matmul(W1_matrix_WE(mesh, j), g))
      do k = 1, 4
        do i = 1, 4
          call increment_msr_matrix_value (sys%a0, mesh%cnode(i,j), mesh%cnode(k,j), a0(i,k))
        end do
      end do
    end do
    
    !! Project out the rows and columns corresponding to Dirichlet edges.
    allocate(ebedge(count(.not.sys%w0mask)))
    ebedge = pack(array=(/(j,j=1,mesh%nnode)/), mask=.not.sys%w0mask)
    call remove_dof (sys%a0, ebedge)
    deallocate(ebedge)
    !ASSERT( is_symmetric(sys%a0) )

  end subroutine initialize_system
  
  subroutine remove_dof (this, dof)

    type(msr_matrix), intent(inout) :: this
    integer, intent(in) :: dof(:)

    integer :: j, m, n, lmn, lnm
    
    ASSERT( this%nrow == this%ncol )

    do j = 1, size(dof)
      m = dof(j)
      if (m < 1 .or. m > this%nrow) cycle
      this%diag(m) = 1.0_r8
      do lmn = this%xadj(m), this%xadj(m+1)-1
        this%nonz(lmn) = 0.0_r8
        n = this%adjncy(lmn)
        do lnm = this%xadj(n), this%xadj(n+1)-1
          if (this%adjncy(lnm) == m) exit
        end do
        INSIST( lnm < this%xadj(n+1) )
        this%nonz(lnm) = 0.0_r8
      end do
    end do

  end subroutine remove_dof

  !!
  !! Residual for the time-discretized system (trapezoid-rule).
  !!
  
  subroutine tr_res (sys, e, b, e0, r)
  
    type(system),  intent(in)  :: sys
    real(kind=r8), intent(in)  :: e(:), e0(:), b(:)
    real(kind=r8), intent(out) :: r(:)
    
    integer :: i, j, k, l1, l2
    real(kind=r8) :: value, el(6), e0l(6), bl(4)
    
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
    real(kind=r8), intent(in)  :: f(:)
    real(kind=r8), intent(out) :: u(:)
    
    integer :: nedge_onP, nnode_onP
    real(kind=r8) :: un(sys%mesh%nnode), rn(sys%mesh%nnode), r(sys%mesh%nedge)
    
    nedge_onP = sys%mesh%nedge_onP
    nnode_onP = sys%mesh%nnode_onP
    
    ASSERT( all(sys%w1mask .or. (f == 0.0_r8)) )
    
    !! Forward Gauss-Seidel relaxation on the on-process edge system.
    u = 0.0_r8
    call gs_relaxation (sys%a1, f(:nedge_onP), u, pattern='f')
    
    !! Update the local residual and project it to the nodes.
    r = f - matmul_msr(sys%a1, u)
    rn = grad_t(sys%mesh, r)
    where (.not.sys%w0mask) rn = 0.0_r8
    
    !! Symmetric Gauss-Seidel relaxation on the projected on-process node system.
    un = 0.0_r8
    call gs_relaxation (sys%a0, rn(:nnode_onP), un, pattern='fb')
        
    !! Update the the solution with the node-based correction.
    u = u + grad(sys%mesh, un)
    
    !! Backward Gauss-Seidel relaxation on the on-process edge system.
    call gs_relaxation (sys%a1, f(:nedge_onP), u, pattern='b')

    call scatter_boundary_sum (sys%mesh%edge_ip, u)
    call gather_boundary (sys%mesh%edge_ip, u)
    
    ASSERT( all(sys%w1mask .or. (u == 0.0_r8)) )
    
  end subroutine hiptmair
  
  !!
  !! These procedures get passed to the CG solver to compute the matrix-vector
  !! product and preconditioning.  Note that they reference the (previously
  !! initialized) module variable CG_SYS which points to the additional data
  !! required by the computations.
  !!
  
  subroutine cg_ax (x, y)
  
    real(kind=r8), intent(in)  :: x(:)
    real(kind=r8), intent(out) :: y(:)
    
    ASSERT( all(cg_sys%w1mask .or. (x == 0.0_r8)) )

    y = matmul_msr(cg_sys%a1, x)
    ASSERT( all(cg_sys%w1mask .or. (y == 0.0_r8)) )
    
    call gather_boundary (cg_sys%mesh%edge_ip, y)
    
  end subroutine cg_ax
  
  
  subroutine cg_pc (x, y)
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    call hiptmair (cg_sys, x, y)
  end subroutine cg_pc
  
  
  function joule_heat (sys, u) result (q)
  
    type(system),  intent(in) :: sys
    real(kind=r8), intent(in) :: u(:)
    real(kind=r8) :: q(sys%mesh%ncell)
    
    integer :: i, j, k, l1, l2
    real(kind=r8) :: value, ul(6), m1(21)
    
    ASSERT( size(u) == sys%mesh%nedge )
    
    q = 0.0_r8
    do j = 1, sys%mesh%ncell
      if (sys%sigma(j) == 0.0_r8) cycle
      do k = 1, 6  ! gather local copy of the vector
        ul(k) = u(sys%mesh%cedge(k,j))
      end do
      m1 = W1_matrix_WE(sys%mesh, j)
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
      q(j) = sys%sigma(j) * q(j) / abs(sys%mesh%volume(j))  ! want a source density
    end do
    
  end function joule_heat
  
  subroutine hiptmair_is_symmetric (sys)
    
    type(system), intent(in), target :: sys
    
    integer :: j, jproc, k, kproc, kstart, lun, nedge_onP, n, m, bsize_vector(nPE)
    real(kind=r8), dimension(sys%mesh%nedge) :: ej, ek, tmp
    real(kind=r8) :: pe(sys%mesh%nedge,sys%mesh%edge_ip%global_size())
    character(len=1) :: edge_type(sys%mesh%nedge), etype(sys%mesh%edge_ip%global_size())
    real(kind=r8) :: ej_dot_pek, ek_dot_pej, ej_dot_pej
    
    nedge_onP = sys%mesh%nedge_onP
    cg_sys => sys
    
    call collate (bsize_vector, sys%mesh%edge_ip%onP_size())
    
    !! Mark the on-PE edges as partition boundary, partition interior, or ignored
    tmp(1:nedge_onP)  = 0.0_r8
    tmp(nedge_onP+1:) = 1.0_r8
    call scatter_boundary_max (sys%mesh%edge_ip, tmp)
    where (tmp > 0.0_r8)
      edge_type = 'B'
    elsewhere
      edge_type = 'I'
    end where
    where (sys%emask /= 0) edge_type = '*'
    call collate (etype, edge_type(:nedge_onP))
    
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
        call gather_boundary (sys%mesh%edge_ip, ej)
        
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
        call gather_boundary (sys%mesh%edge_ip, ej)
        
        m = n
        do kproc = jproc, nPE
          kstart = 1
          if (kproc == jproc) kstart = j + 1
          do k = kstart, bsize_vector(kproc)
            m = m + 1

            !! Initialize the (kproc, k) unit vector
            ek = 0.0_r8
            if (this_PE == kproc) ek(k) = 1.0_r8
            call gather_boundary (sys%mesh%edge_ip, ek)

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
      call collate (tmp, (/ value /))
      if (is_IOP) value = tmp(pe)
      call broadcast (value)
    end subroutine broadcast_from_PE
  
  end subroutine hiptmair_is_symmetric
  
end module MaxwellEddy
