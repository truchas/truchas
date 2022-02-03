!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_hypre_hybrid_type

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds
  use pgslib_module
  use parallel_communication
  use index_partitioning
  use pcsr_matrix_type
  use hypre_hybrid_type
  use parameter_list_type
  implicit none
  
  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()
  real(r8), parameter :: PI = 3.1415926535897931_r8
  real(r8), parameter :: TWOPI = 6.2831853071795862_r8
  integer :: status = 0
  
  real(r8) :: a ! set by the tests
  integer :: nx, ny, nz ! set by the tests
  
  call init_parallel_communication(argv)

  call cg_test_1
  call cg_test_2
  call gmres_test_1
  call gmres_test_2

  call pgslib_finalize
  
  call exit (status)

contains

  subroutine cg_test_1
  
    type(pcsr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    if (is_IOP) write(*,'(/,a)') 'Running CG_TEST_1'

    a = 1.0_r8
    nx = 13; ny = 17; nz = 19

    allocate(params)
    call params%set ('krylov-method', 'cg')
    call params%set ('rel-tol', 1.0e-8_r8)
    call params%set ('conv-rate-tol', 0.8_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('cg-use-two-norm', .true.)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 0)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%graph%row_ip%onP_size()
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; if (is_IOP) x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      if (is_IOP) write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = global_maxval(abs(u-x))
    l2err = sqrt(global_sum((u-x)**2))
    if (is_IOP) write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    if (is_IOP) write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 26) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected 26 diagonally-scaled CG iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr /= 0) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected 0 AMG preconditioned CG iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 2.0e-8) then
      if (is_IOP) write(*,'(a,es9.2)') 'error: expected maxerr < 2.0e-8; got', maxerr
      status = 1
    end if
  
  end subroutine


  subroutine cg_test_2
  
    type(pcsr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    if (is_IOP) write(*,'(/,a)') 'Running CG_TEST_2'

    a = 1.0e-4_r8
    nx = 13; ny = 17; nz = 19

    allocate(params)
    call params%set ('krylov-method', 'cg')
    call params%set ('rel-tol', 1.0e-8_r8)
    call params%set ('conv-rate-tol', 0.6_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    !call params%set ('cg-use-two-norm', .true.)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 0)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%graph%row_ip%onP_size()
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; if (is_IOP) x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      if (is_IOP) write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = global_maxval(abs(u-x))
    l2err = sqrt(global_sum((u-x)**2))
    if (is_IOP) write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    if (is_IOP) write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 18) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected 18 diagonally-scaled CG iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr > 5) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected no more than 5 AMG preconditioned CG iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 2.0e-8) then
      if (is_IOP) write(*,'(a,es9.2)') 'error: expected maxerr < 2.0e-8; got', maxerr
      status = 1
    end if
  
  end subroutine


  subroutine gmres_test_1
  
    type(pcsr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    if (is_IOP) write(*,'(/,a)') 'Running GMRES_TEST_1'

    a = 1.0_r8
    nx = 13; ny = 17; nz = 19

    allocate(params)
    call params%set ('krylov-method', 'gmres')
    call params%set ('rel-tol', 1.0e-6_r8)
    call params%set ('conv-rate-tol', 0.8_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('gmres-krylov-dim', 5)
    call params%set ('amg-smoothing-sweeps', 1)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 1)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%graph%row_ip%onP_size()
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; if (is_IOP) x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      if (is_IOP) write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = global_maxval(abs(u-x))
    l2err = sqrt(global_sum((u-x)**2))
    if (is_IOP) write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    if (is_IOP) write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 19) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected 19 diagonally-scaled GMRES iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr /= 0) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected 0 AMG preconditioned GMRES iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 2.0e-6) then
      if (is_IOP) write(*,'(a,es9.2)') 'error: expected maxerr < 2.0e-6; got', maxerr
      status = 1
    end if
    
  end subroutine


  subroutine gmres_test_2
  
    type(pcsr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    if (is_IOP) write(*,'(/,a)') 'Running GMRES_TEST_2'

    a = 1.0e-4_r8
    nx = 13; ny = 17; nz = 19

    allocate(params)
    call params%set ('krylov-method', 'gmres')
    call params%set ('rel-tol', 1.0e-6_r8)
    call params%set ('conv-rate-tol', 0.6_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('gmres-krylov-dim', 5)
    call params%set ('amg-smoothing-sweeps', 1)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 1)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%graph%row_ip%onP_size()
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; if (is_IOP) x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      if (is_IOP) write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = global_maxval(abs(u-x))
    l2err = sqrt(global_sum((u-x)**2))
    if (is_IOP) write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    if (is_IOP) write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 2) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected 2 diagonally-scaled GMRES iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr > 7) then
      if (is_IOP) write(*,'(a,i0)') 'error: expected no more than 7 AMG preconditioned GMRES iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 2.0e-6) then
      if (is_IOP) write(*,'(a,es9.2)') 'error: expected maxerr < 2.0e-6; got', maxerr
      status = 1
    end if
    
  end subroutine


  subroutine create_matrix (matrix)

    type(pcsr_matrix), intent(out) :: matrix
    
    integer :: ix, iy, iz, n, j, k, ntot, nloc
    integer, allocatable :: nnbr_g(:,:), nnbr(:,:), offP_index(:)
    type(ip_desc), pointer :: row_ip
    type(pcsr_graph), pointer :: graph
    
    !! Stencil neighbors of each grid point (GLOBAL).
    if (is_IOP) then
      allocate(nnbr_g(6,NX*NY*NZ))
      do iz = 0, NZ-1
        do iy = 0, NY-1
          do ix = 0, NX-1
            n = linear_index(ix,iy,iz)
            nnbr_g(1,n) = linear_index(ix,iy,iz-1)
            nnbr_g(2,n) = linear_index(ix,iy-1,iz)
            nnbr_g(3,n) = linear_index(ix-1,iy,iz)
            nnbr_g(4,n) = linear_index(ix+1,iy,iz)
            nnbr_g(5,n) = linear_index(ix,iy+1,iz)
            nnbr_g(6,n) = linear_index(ix,iy,iz+1)
            !print *, n, ':', nnbr_g(:,n)
          end do
        end do
      end do
    else
      allocate(nnbr_g(6,0))
    end if
    
    !! Partition the unknowns.
    ntot = NX*NY*NZ
    nloc = ntot / nPE
    if (this_PE <= modulo(ntot,nPE)) nloc = nloc + 1
    !print *, 'RANK=', this_PE, ', NLOC=', nloc
    
    !! Setup the index partition for the grid points.
    allocate(row_ip)
    call row_ip%init (nloc)
    call localize_index_array (nnbr_g, row_ip, row_ip, nnbr, offP_index)
    call row_ip%add_offP_index (offP_index)
    deallocate(offP_index, nnbr_g)
    
    !! Create the parallel CSR matrix.
    allocate(graph)
    call graph%init (row_ip)
    do j = 1, size(nnbr,2)
      call graph%add_edge (j, j)
      call graph%add_edge (j, nnbr(:,j))
    end do
    call graph%add_complete
    call matrix%init (graph, take_graph=.true.)
    do j = 1, size(nnbr,2)
      call matrix%set (j, j, a + 6.0_r8)
      do k = 1, size(nnbr,1)
        call matrix%set (j, nnbr(k,j), -1.0_r8)
      end do
    end do
    deallocate(nnbr)
    
  end subroutine create_matrix
  
  integer function linear_index (ix, iy, iz)
    integer, intent(in) :: ix, iy, iz
    linear_index = 1 + modulo(ix,NX) + NX*(modulo(iy,NY) + NY*modulo(iz,NZ))
  end function linear_index

  function eigenvalue (kx, ky, kz) result (lambda)
    integer, intent(in) :: kx, ky, kz
    real(r8) :: lambda
    lambda = a + 4*(sin(PI*kx/NX)**2 + sin(PI*ky/NY)**2 + sin(PI*kz/NZ)**2)
  end function
  
  subroutine eigenvector (kx, ky, kz, u)
  
    integer, intent(in) :: kx, ky, kz
    real(r8), intent(out) :: u(:)
    
    integer :: ix, iy, iz, n
    real(r8) :: dx, dy, dz
    real(r8), allocatable :: u_g(:)
    
    dx = real(kx,kind=r8) / NX
    dy = real(ky,kind=r8) / NY
    dz = real(kz,kind=r8) / NZ

    if (is_IOP) then
      allocate(u_g(NX*NY*NZ))
      do iz = 0, NZ-1
        do iy = 0, NY-1
          do ix = 0, NX-1
            n = linear_index(ix,iy,iz)
            u_g(n) = cos(TWOPI*(dx*ix + dy*iy + dz*iz))
            !print *, n, ':', u_g(n)
          end do
        end do
      end do
    else
      allocate(u_g(0))
    end if
    
    call distribute (u, u_g)
    
  end subroutine eigenvector

end program test_hypre_hybrid_type
