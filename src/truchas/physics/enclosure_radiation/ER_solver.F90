#include "f90_assert.fpp"

module ER_solver

  use kinds, only: r8
  use ER_system
  use ER_dist_encl
  use ER_encl_func
  use scalar_functions
  use parallel_communication
  use truchas_logging_services
  implicit none
  private

  public :: ERS_solver_create, ERS_solver_destroy
  public :: ERS_set_ambient, ERS_set_emissivity, ERS_set_absolute_zero, ERS_set_stefan_boltzmann
  public :: ERS_set_solver_controls
  public :: ERS_compute_residual, ERS_compute_heat_flux, ERS_solve_radiosity
  public :: ERS_precon_cheby, ERS_precon_jacobi, ERS_precon_matvec1
  public :: ERS_rhs_deriv, ERS_rhs
  public :: ERS_set_cheby_param

  type, public :: solver
    integer :: nface
    type(dist_encl), pointer :: encl => null()  ! radiation enclosure surface
    type(scafun),    pointer :: tamb => null()  ! ambient temperature function
    type(encl_func), pointer :: eps  => null()  ! emissivity enclosure function
    type(system) :: sys
    !! Chebyshev solver parameters
    real(r8) :: c, d = -1.0_r8  ! iteration parameters
    real(r8) :: tol
    integer  :: maxitr
  end type solver

contains

  subroutine ERS_solver_create (this, file, csf, color)

    use ER_file
    use index_partitioning

    type(solver), intent(out) :: this
    character(len=*), intent(in) :: file
    real(r8), intent(in) :: csf
    integer, intent(in) :: color(:)

    integer :: j, ncid, offset, bsize(nPE)
    logical :: renumbered

    if (is_IOP) call ERF_open_ro (file, ncid)

    !! Load the distributed enclosure surface.
    allocate(this%encl)
    call create_dist_encl (this%encl, ncid, color)

    if (csf /= 1.0_r8) this%encl%coord = csf * this%encl%coord

    this%nface = this%encl%nface_onP

    !! The current implementation of LOAD_VIEW_FACTORS doesn't allow for
    !! renumbering the enclosure faces, so check that no renumbering has
    !! occurred.  This should be the case if COLOR is a blocked coloring.
    renumbered = .false.
    offset = first_index(this%encl%face_ip) - 1
    do j = 1, this%encl%nface_onP
      if (this%encl%face_map(j) /= j+offset) renumbered = .true.
    end do
    INSIST(.not.global_any(renumbered))

    !! Load the distributed view factor matrix.
    call collate (bsize, this%nface)
    call broadcast (bsize)
    call load_view_factors (this%sys, ncid, bsize)

    if (is_IOP) call ERF_close (ncid)

  end subroutine ERS_solver_create

  subroutine load_view_factors (this, ncid, bsize)

    use ER_file

    type(system), intent(out) :: this
    integer, intent(in) :: ncid
    integer, intent(in) :: bsize(:)

    integer :: j, n, nface, nface_tot, start
    integer :: vf_bsize(nPE), lengths(nPE), idum0(0)
    real :: rdum0(0)
    integer, pointer :: ibuf(:) => null()
    real,    pointer :: rbuf(:) => null()

    if (is_IOP) call ERF_get_vf_dims (ncid, nface_tot, n)
    call broadcast (nface_tot)
    nface = bsize(this_PE)

    this%nface = nface
    this%offset = sum(bsize(:this_PE-1))
    this%nface_tot = nface_tot

    ASSERT( nface_tot == sum(bsize) )

    !! Read and distribute the ambient view factors.
    call allocate_collated_array (rbuf, this%nface_tot)
    if (is_IOP) call ERF_get_ambient (ncid, rbuf)
    allocate(this%amb_vf(this%nface))
    call distribute (this%amb_vf, rbuf)
    deallocate(rbuf)

    !! Read and distribute the VF matrix nonzero row counts.
    call allocate_collated_array (ibuf, this%nface_tot)
    if (is_IOP) call ERF_get_vf_rowcount (ncid, ibuf)
    allocate(this%ia(this%nface+1))
    call distribute (this%ia(2:), ibuf)
    deallocate(ibuf)

    !! Convert the row counts into the local IA indexing array.
    this%ia(1) = 1
    do j = 2, size(this%ia)
      this%ia(j) = this%ia(j) + this%ia(j-1)
    end do

    !! Determine the sizes of the distributed VF matrix.
    n = this%ia(this%nface+1) - this%ia(1)
    allocate(this%vf(n), this%ja(n))
    call collate (vf_bsize, n)
    !call broadcast (vf_bsize)

    !! Read the VF matrix in process-sized blocks, sending them to the owning
    !! processes as we go.  PGSLib only provides the distribute collective to
    !! do this, and so we distribute with 0-sized destination arrays for all
    !! processes but the receiving one.  We also use the optional LENGTHS
    !! argument to distribute (only referenced on the IO process) that gives
    !! the number of items to distribute to each process, resulting in the
    !! actual sizes of the array arguments being ignored except to check that
    !! they are sufficently large.  This allows us to simplify the calls.
    !! The code is structured for future use of MPI_Isend/Irecv, abandoning
    !! the usual 'single code path' pattern.

    if (is_IOP) then
      call ERF_get_vf_rows (ncid, this%vf, this%ja, start=1)
      if (nPE > 1) then
        n = maxval(vf_bsize(2:))
        allocate(ibuf(n), rbuf(n))
        lengths = 0
        start = 1 + vf_bsize(1)
        do n = 2, nPE
          call ERF_get_vf_rows (ncid, rbuf(:vf_bsize(n)), ibuf(:vf_bsize(n)), start)
          lengths(n) = vf_bsize(n)
          call distribute (rdum0, rbuf, lengths)
          call distribute (idum0, ibuf, lengths)
          lengths(n) = 0
          start = start + vf_bsize(n)
        end do
        deallocate(ibuf, rbuf)
      end if
    else  ! everybody else just participates in the distribute calls.
      do n = 2, nPE
        if (n == this_PE) then
          call distribute (this%vf, rdum0, lengths)
          call distribute (this%ja, idum0, lengths)
        else
          call distribute (rdum0, rdum0, lengths)
          call distribute (idum0, idum0, lengths)
        end if
      end do
    end if

  end subroutine load_view_factors

  subroutine ERS_solver_destroy (this)
    type(solver), intent(inout) :: this
    if (associated(this%encl)) then
      call destroy_dist_encl (this%encl)
      deallocate(this%encl)
    end if
    if (associated(this%tamb)) then
      call destroy (this%tamb)
      deallocate(this%tamb)
    end if
    if (associated(this%eps)) then
      call EF_destroy (this%eps)
      deallocate(this%eps)
    end if
    call system_destroy (this%sys)
  end subroutine ERS_solver_destroy

  subroutine ERS_set_ambient (this, tamb)
    type(solver), intent(inout) :: this
    type(scafun), pointer :: tamb
    this%tamb => tamb
  end subroutine ERS_set_ambient

  subroutine ERS_set_emissivity (this, eps)
    type(solver), intent(inout) :: this
    type(encl_func), pointer :: eps
    this%eps => eps
  end subroutine ERS_set_emissivity

  subroutine ERS_set_absolute_zero (this, t0)
    type(solver), intent(inout) :: this
    real(r8), intent(in) :: t0
    this%sys%t0 = t0
  end subroutine ERS_set_absolute_zero

  subroutine ERS_set_stefan_boltzmann (this, sigma)
    type(solver), intent(inout) :: this
    real(r8), intent(in) :: sigma
    this%sys%sigma = sigma
  end subroutine ERS_set_stefan_boltzmann

  subroutine ERS_set_solver_controls (this, tol, maxitr)
    type(solver), intent(inout) :: this
    real(r8), intent(in) :: tol
    integer, intent(in) :: maxitr
    this%tol = tol
    this%maxitr = maxitr
  end subroutine ERS_set_solver_controls

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERS_SOLVE_RADIOSITY
 !!
 !! This subroutine solves the radiosity system for the radiosity from the
 !! enclosure surface given the temperature on the surface.  The surface
 !! emissivities and ambient temperature parameters of the system may be
 !! time-dependent.
 !!
 !! IMPLEMENTATION NOTES
 !!
 !! (1) The Chebyshev iteration parameters c and d must be chosen so that the
 !! spectrum of the radiosity system matrix is enclosed by an ellipse in the
 !! right half plane from the family of ellipses with foci d±c.  Parameters
 !! producing the smallest such ellipse give optimal convergence.  Since
 !! the spectrum is contained in the interval (0,2), the degenerate ellipse
 !! with the min and max eigenvalues as foci is optimal; that is,
 !!
 !!         c = (lmax - lmin)/2,   d = (lmax + lmin)/2.
 !!
 !! While lmin is easily computed in practice using the power method, lmax
 !! is not.  Thus we instead choose the foci lmin and 2 - lmin, which ensures
 !! that (0, 2) is covered with a family of right hand plane ellipses.  The
 !! resulting convergence rate is diminished but not impractically so.
 !!
 !! (2)
 !! (2) The spectrum of the radiosity system is independent of its size,
 !! depending only on the enclosure geometry and the emissivity of the
 !! surface.  Thus we


  subroutine ERS_solve_radiosity (this, time, temp, qrad)

    type(solver), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: temp(:)
    real(r8), intent(inout) :: qrad(:)

    integer  :: n
    real(r8) :: err, tamb

    ASSERT(size(temp) == this%nface)
    ASSERT(size(qrad) == this%nface)

    tamb = eval(this%tamb, (/time/))
    call EF_eval(this%eps, time)

    call system_cheby_solve (this%sys, temp, tamb, this%eps%values, &
                             this%c, this%d, this%tol, this%maxitr, qrad, n, err)
    if (is_IOP) write(*,'(2x,a,i0,a,es9.3)') 'rhs-relative error(', n, ')=', err
    INSIST( n <= this%maxitr )

  end subroutine ERS_solve_radiosity

  subroutine ERS_precon_cheby (this, time, numitr, z)

    type(solver), intent(inout) :: this
    real(r8), intent(in) :: time
    integer, intent(in) :: numitr
    real(r8), intent(inout) :: z(:)

    ASSERT(size(z) == this%nface)

    call EF_eval(this%eps, time)
    call system_cheby_precon (this%sys, this%eps%values, this%c, this%d, numitr, z)

  end subroutine ERS_precon_cheby

  subroutine ERS_precon_jacobi (this, time, numitr, z)

    type(solver), intent(inout) :: this
    real(r8), intent(in) :: time
    integer, intent(in) :: numitr
    real(r8), intent(inout) :: z(:)

    ASSERT(size(z) == this%nface)

    call EF_eval(this%eps, time)
    call system_jacobi_precon (this%sys, this%eps%values, numitr, z)

  end subroutine ERS_precon_jacobi

  subroutine ERS_precon_matvec1 (this, time, z)

    type(solver), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)

    ASSERT(size(z) == this%nface)

    call EF_eval(this%eps, time)
    call system_matvec1 (this%sys, this%eps%values, z)

  end subroutine ERS_precon_matvec1

  subroutine ERS_set_cheby_param (this, time)

    type(solver), intent(inout) :: this
    real(r8), intent(in) :: time

    integer  :: n, maxitr
    real(r8) :: lmin, tol, err
    real(r8) :: z(this%nface)
    integer, allocatable :: prn_seed(:)
    character(len=255) :: msg

    call TLS_info ('    Calculating Chebyshev iteration parameters ...')

    call EF_eval(this%eps, time)

    !! Set the PRN generator seed.  This gives us reproducibility on the same
    !! platform using the same number of processes/data distribution.  But we
    !! will get different results on different platforms.
    call random_seed (size=n)
    allocate(prn_seed(n))
    prn_seed = -354891804
    call random_seed (put=prn_seed)
    deallocate(prn_seed)
    call random_number (z)  ! contains a component of the min eigenvector

    !! Chebyshev iteration parameters C and D (see Notes 1, 2).
    tol = 1.0d-4; maxitr = 10000 !TODO! MAKE THESE PARAMETERS
    call system_lambda_min (this%sys, this%eps%values, tol, maxitr, z, lmin, n, err)
    write(msg,'(6x,a,i0,a,f8.6,a,es9.3)') 'eigenvalue calculation: lmin(', n, ')=', lmin,  ', error=', err
    call TLS_info (trim(msg))
    
    !! It is not expected that the preceding eigenvalue calculation should ever
    !! not converge.  But if it doesn't we'll plunge ahead, using the most
    !! pessimistic value, and write a warning message to that effect.
    !INSIST(n <= maxitr)
    if (n > maxitr) then
      call TLS_warn ('Min eigenvalue calculation failed to converge; setting lmin=0 and continuing.')
      lmin = 0.0_r8
    end if

    this%d = 1.0_r8
    this%c = 1.0_r8 - lmin

    write(msg,'(6x,2(a,f8.6))') 'setting d=', this%d, ', c=', this%c
    call TLS_info (trim(msg))

  end subroutine ERS_set_cheby_param

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERS_COMPUTE_HEAT_FLUX
 !!
 !! Calculates the net heat flux through the enclosure surface given the
 !! radiosity from the enclosure surface.  The flux may depend on a possibly
 !! time-dependent ambient temperature.
 !!

  subroutine ERS_compute_heat_flux (this, time, qrad, temp, flux)

    type(solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: qrad(:)  ! local radiosity vector
    real(r8), intent(in)  :: temp(:)  ! local temperature vector
    real(r8), intent(out) :: flux(:)  ! local heat flux vector

    real(r8) :: tamb

    ASSERT(size(qrad) == this%nface)
    ASSERT(size(temp) == this%nface)
    ASSERT(size(flux) == this%nface)

    tamb = eval(this%tamb, (/time/))
    call EF_eval (this%eps, time)
    call system_flux (this%sys, qrad, temp, tamb, this%eps%values, flux)

  end subroutine ERS_compute_heat_flux

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERS_COMPUTE_RESIDUAL
 !!
 !! This subroutine calculates the residual of the radiosity system given the
 !! temperature and radiosity on the enclosure surface.  The emissivities and
 !! ambient temperature parameters of the system may be time-dependent.
 !!

  subroutine ERS_compute_residual (this, time, qrad, temp, res)

    type(solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: qrad(:)  ! radiosity vector
    real(r8), intent(in)  :: temp(:)  ! temperature vector
    real(r8), intent(out) :: res(:)   ! residual vector

    real(r8) :: tamb

    ASSERT(size(qrad) == this%nface)
    ASSERT(size(temp) == this%nface)
    ASSERT(size(res)  == this%nface)

    tamb = eval(this%tamb, (/time/))
    call EF_eval (this%eps, time)
    call system_residual (this%sys, qrad, temp, tamb, this%eps%values, res)

  end subroutine ERS_compute_residual

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERS_RHS_DERIV
 !!
 !! The subroutine calculates the Jacobian of the radiosity system RHS with
 !! respect to the enclosure surface temperatures.  This Jacobian is diagonal
 !! and the diagonal is returned in drhs.
 !!

  subroutine ERS_rhs_deriv (this, time, temp, drhs)

    type(solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: temp(:)  ! temperature vector
    real(r8), intent(out) :: drhs(:)  ! RHS derivative (diagonal)

    ASSERT(size(temp) == this%nface)
    ASSERT(size(drhs) == this%nface)

    call EF_eval (this%eps, time)
    call system_rhs_deriv (this%sys, temp, this%eps%values, drhs)

  end subroutine ERS_rhs_deriv

  subroutine ERS_rhs (this, time, temp, rhs)

    type(solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: temp(:)  ! temperature vector
    real(r8), intent(out) :: rhs(:)

    real(r8) :: tamb

    ASSERT(size(temp) == this%nface)
    ASSERT(size(rhs) == this%nface)

    tamb = eval(this%tamb, (/time/))
    call EF_eval (this%eps, time)
    call system_rhs (this%sys, temp, tamb, this%eps%values, rhs)

  end subroutine ERS_rhs

end module ER_solver
