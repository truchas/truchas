!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_hypre_pcg_type

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds
  use parallel_util_module, only: parallel_init
  use pgslib_module
  use parallel_communication
  use index_partitioning
  use pcsr_matrix_type
  use hypre_pcg_type
  use parameter_list_type
  implicit none
  
  !integer, parameter :: NX = 128, NY = 128, NZ = 128
  integer, parameter :: NX = 19, NY = 13, NZ = 17
  !integer, parameter :: NX = 7, NY = 5, NZ = 3
  
  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()
  type(pcsr_matrix), target :: matrix
  type(hypre_pcg) :: solver
  type(parameter_list), pointer :: params
  real(r8), allocatable :: x(:), b(:), u(:)
  integer :: nrow, kx, ky, kz, num_itr, stat, maxitr
  real(r8) :: a, maxerr, l2err, rtol
  real(r8), parameter :: PI = 3.1415926535897931_r8
  
  call parallel_init (argv)
  call init_parallel_communication
  
  a = 1.0e-6_r8
  call create_matrix (a, matrix)
  allocate(params)
  rtol = 1.0e-12_r8
  call params%set ('rel-tol', rtol)
  call params%set ('max-iter', 100)
  call params%set ('print-level', 2)
  call params%set ('amg-num-cycles', 1)
  call solver%init (matrix, params)
  call solver%setup ()
  
  nrow = matrix%graph%row_ip%onP_size()
  allocate(x(nrow), b(nrow), u(nrow))
  
  kx = 1; ky = 1; kz = 1
  call solution (kx, ky, kz, u)
  b = (a + 4*(sin(PI*kx/NX)**2 + sin(PI*ky/NY)**2 + sin(PI*kz/NZ)**2)) * u
  x = 0.0_r8
  call solver%solve (b, x, stat)
  if (stat /= 0) then
    if (is_IOP) print *, 'PCG failed to converge'
    call exit (1)
  end if
  call solver%get_metrics (num_itr)

  maxerr = global_maxval(abs(u-x))
  l2err = sqrt(global_sum((u-x)**2) / matrix%graph%row_ip%global_size())
  if (is_IOP) print '(a,i2,2(a,es9.2),/)', 'itr=', num_itr, ', maxerr=', maxerr, ', l2err=', l2err

  stat = 0
  maxitr = 11
  if (num_itr <= maxitr) then
    if (is_IOP) print '(a,i0)', 'PASS: expected number of iterations <= ', maxitr
  else
    stat = 1
    if (is_IOP) print '(a,i0)', 'FAIL: expected number of iterations <= ', maxitr
  end if

  maxerr = 20*rtol
  if (l2err <= maxerr) then
    if (is_IOP) print '(a,es10.3)', 'PASS: expected l2 error <=', maxerr
  else
    stat = 1
    if (is_IOP) print '(a,es10.3)', 'FAIL: expected l2 error <=', maxerr
  end if

  call pgslib_finalize

  call exit(stat)

contains

  subroutine create_matrix (a, matrix)

    real(r8), intent(in) :: a
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
    print *, 'RANK=', this_PE, ', NLOC=', nloc
    
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
    
  end subroutine create_matrix
  
  integer function linear_index (ix, iy, iz)
    integer, intent(in) :: ix, iy, iz
    linear_index = 1 + modulo(ix,NX) + NX*(modulo(iy,NY) + NY*modulo(iz,NZ))
  end function linear_index
  
  subroutine solution (kx, ky, kz, u)
  
    integer, intent(in) :: kx, ky, kz
    real(r8), intent(out) :: u(:)
    
    integer :: ix, iy, iz, n
    real(r8) :: dx, dy, dz
    real(r8), allocatable :: u_g(:)
    real(r8), parameter :: TWOPI = 6.2831853071795862_r8
    
    dx = real(kx,kind=r8) / NX
    dy = real(ky,kind=r8) / NY
    dz = real(kz,kind=r8) / NZ
    
    if (is_IOP) then
      allocate(u_g(NX*NY*NZ))
      do iz = 0, NZ-1
        do iy = 0, NY-1
          do ix = 0, NX-1
            n = linear_index(ix,iy,iz)
            u_g(n) = sin(TWOPI*(dx*ix + dy*iy + dz*iz))
            !print *, n, ':', u_g(n)
          end do
        end do
      end do
    else
      allocate(u_g(0))
    end if
    
    call distribute (u, u_g)
    
  end subroutine solution

end program test_hypre_pcg_type
