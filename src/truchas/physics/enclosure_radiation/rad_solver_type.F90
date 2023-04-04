!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module rad_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use rad_system_type
  use rad_encl_func_type
  use scalar_func_class
  use parallel_communication
  use truchas_logging_services
  implicit none
  private

  type, public :: rad_solver
    integer :: nface
    class(scalar_func), allocatable :: tamb  ! ambient temperature function
    type(rad_encl_func), allocatable :: eps  ! emissivity enclosure function
    type(rad_system) :: sys
    !! Chebyshev solver parameters
    real(r8) :: c, d = -1.0_r8  ! iteration parameters
    real(r8) :: tol
    integer  :: maxitr
  contains
    procedure :: init
    procedure :: set_ambient
    procedure :: set_emissivity
    procedure :: set_absolute_zero
    procedure :: set_stefan_boltzmann
    procedure :: set_solver_controls
    procedure :: set_cheby_param
    procedure :: residual
    procedure :: heat_flux
    procedure :: rhs
    procedure :: rhs_deriv
    procedure :: precon_cheby
    procedure :: precon_jacobi
    procedure :: precon_matvec1
    procedure :: solve_radiosity
  end type rad_solver

contains

  subroutine init(this, vf)
    use encl_vf_class
    class(rad_solver), intent(out), target :: this
    class(encl_vf), intent(in), target :: vf
    this%nface = vf%nface
    call this%sys%init(vf)
  end subroutine init

  subroutine set_ambient (this, tamb)
    class(rad_solver), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: tamb
    call move_alloc (tamb, this%tamb)
  end subroutine set_ambient

  subroutine set_emissivity (this, eps)
    class(rad_solver), intent(inout) :: this
    type(rad_encl_func), allocatable, intent(inout) :: eps
    call move_alloc(eps, this%eps)
  end subroutine set_emissivity

  subroutine set_absolute_zero (this, t0)
    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: t0
    this%sys%t0 = t0
  end subroutine set_absolute_zero

  subroutine set_stefan_boltzmann (this, sigma)
    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: sigma
    this%sys%sigma = sigma
  end subroutine set_stefan_boltzmann

  subroutine set_solver_controls (this, tol, maxitr)
    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: tol
    integer, intent(in) :: maxitr
    this%tol = tol
    this%maxitr = maxitr
  end subroutine set_solver_controls

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOLVE_RADIOSITY
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


  subroutine solve_radiosity (this, time, temp, qrad, stat, numitr, error)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: temp(:)
    real(r8), intent(inout) :: qrad(:)
    integer,  intent(out) :: stat
    integer,  intent(out) :: numitr
    real(r8), intent(out) :: error

    real(r8) :: tamb

    ASSERT(size(temp) == this%nface)
    ASSERT(size(qrad) == this%nface)

    tamb = this%tamb%eval ([time])
    call this%eps%eval (time, temp)

    call this%sys%cheby_solve (temp, tamb, this%eps%values, &
                               this%c, this%d, this%tol, this%maxitr, qrad, numitr, error)
    stat = merge(0, 1, numitr <= this%maxitr)

  end subroutine solve_radiosity

  subroutine precon_cheby (this, time, temp, numitr, z)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: time, temp(:)
    integer, intent(in) :: numitr
    real(r8), intent(inout) :: z(:)

    ASSERT(size(temp) == this%nface)
    ASSERT(size(z) == this%nface)

    call this%eps%eval (time, temp)
    call this%sys%cheby_precon (this%eps%values, this%c, this%d, numitr, z)

  end subroutine precon_cheby

  subroutine precon_jacobi (this, time, temp, numitr, z)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: time, temp(:)
    integer, intent(in) :: numitr
    real(r8), intent(inout) :: z(:)

    ASSERT(size(temp) == this%nface)
    ASSERT(size(z) == this%nface)

    call this%eps%eval (time, temp)
    call this%sys%jacobi_precon (this%eps%values, numitr, z)

  end subroutine precon_jacobi

  subroutine precon_matvec1 (this, time, temp, z)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: time, temp(:)
    real(r8), intent(inout) :: z(:)

    ASSERT(size(temp) == this%nface)
    ASSERT(size(z) == this%nface)

    call this%eps%eval (time, temp)
    call this%sys%matvec1 (this%eps%values, z)

  end subroutine precon_matvec1

  subroutine set_cheby_param (this, time, temp)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in) :: time, temp(:)

    integer, parameter ::  maxitr = 10000
    real(r8), parameter :: tol = 1.0d-4
    integer  :: n
    real(r8) :: lmin, err
    real(r8) :: z(this%nface)
    integer, allocatable :: prn_seed(:)
    character(len=255) :: msg

    call TLS_info ('    calculating Chebyshev iteration parameters ...')

    call this%eps%eval (time, temp)

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
    call this%sys%lambda_min (this%eps%values, tol, maxitr, z, lmin, n, err)
    write(msg,'(6x,a,i0,a,f8.6,a,es10.3)') 'eigenvalue calculation: lmin(', n, ')=', lmin,  ', error=', err
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

  end subroutine set_cheby_param

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! HEAT_FLUX
 !!
 !! Calculates the net heat flux through the enclosure surface given the
 !! radiosity from the enclosure surface.  The flux may depend on a possibly
 !! time-dependent ambient temperature.
 !!

  subroutine heat_flux (this, time, qrad, temp, flux)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: qrad(:)  ! local radiosity vector
    real(r8), intent(in)  :: temp(:)  ! local temperature vector
    real(r8), intent(out) :: flux(:)  ! local heat flux vector

    real(r8) :: tamb

    ASSERT(size(qrad) == this%nface)
    ASSERT(size(temp) == this%nface)
    ASSERT(size(flux) == this%nface)

    tamb = this%tamb%eval([time])
    call this%eps%eval (time, temp)
    call this%sys%flux (qrad, temp, tamb, this%eps%values, flux)

  end subroutine heat_flux

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RESIDUAL
 !!
 !! This subroutine calculates the residual of the radiosity system given the
 !! temperature and radiosity on the enclosure surface.  The emissivities and
 !! ambient temperature parameters of the system may be time-dependent.
 !!

  subroutine residual (this, time, qrad, temp, res)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: qrad(:)  ! radiosity vector
    real(r8), intent(in)  :: temp(:)  ! temperature vector
    real(r8), intent(out) :: res(:)   ! residual vector

    real(r8) :: tamb

    ASSERT(size(qrad) == this%nface)
    ASSERT(size(temp) == this%nface)
    ASSERT(size(res)  == this%nface)

    tamb = this%tamb%eval([time])
    call this%eps%eval (time, temp)
    call this%sys%residual (qrad, temp, tamb, this%eps%values, res)

  end subroutine residual

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERS_RHS_DERIV
 !!
 !! The subroutine calculates the Jacobian of the radiosity system RHS with
 !! respect to the enclosure surface temperatures.  This Jacobian is diagonal
 !! and the diagonal is returned in drhs.
 !!

  subroutine rhs_deriv (this, time, temp, drhs)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: temp(:)  ! temperature vector
    real(r8), intent(out) :: drhs(:)  ! RHS derivative (diagonal)

    real(r8) :: tamb
    real(r8), allocatable :: fdinc(:), deps(:)

    ASSERT(size(temp) == this%nface)
    ASSERT(size(drhs) == this%nface)

    fdinc = max(1.0_r8, abs(temp)) * sqrt(epsilon(1.0_r8))
    call this%eps%eval (time, temp + fdinc)
    deps = this%eps%values
    call this%eps%eval (time, temp - fdinc)
    deps = (deps - this%eps%values) / (2*fdinc)

    tamb = this%tamb%eval([time])
    call this%eps%eval (time, temp)
    call this%sys%rhs_deriv (temp, tamb, this%eps%values, deps, drhs)

  end subroutine rhs_deriv

  subroutine rhs (this, time, temp, b)

    class(rad_solver), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: temp(:)  ! temperature vector
    real(r8), intent(out) :: b(:)

    real(r8) :: tamb

    ASSERT(size(temp) == this%nface)
    ASSERT(size(b) == this%nface)

    tamb = this%tamb%eval([time])
    call this%eps%eval (time, temp)
    call this%sys%rhs (temp, tamb, this%eps%values, b)

  end subroutine rhs

end module rad_solver_type
