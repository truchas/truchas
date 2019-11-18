!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module rad_system_type

  use kinds
  use parallel_communication
  implicit none
  private

  type, public :: rad_system
    real(r8) :: sigma ! Stefan-Boltzmann constant
    real(r8) :: t0    ! absolute-zero temperature
    integer :: nface      ! number of faces (number of rows) on this process
    integer :: offset     ! difference between global index and local row index
    integer :: nface_tot  ! total number of faces (number of columns)
    integer, allocatable :: ia(:), ja(:)
    real, allocatable :: vf(:), amb_vf(:)
  contains
    procedure :: flux
    procedure :: residual
    procedure :: matvec1
    procedure :: rhs
    procedure :: rhs_deriv
    procedure :: cheby_solve
    procedure :: cheby_precon
    procedure :: jacobi_precon
    procedure :: lambda_min
    procedure :: lambda_max
  end type rad_system

contains

  !!
  !! Heat flux computational kernel.
  !!

  subroutine flux (this, q, t, tamb, eps, f)

    class(rad_system), intent(in) :: this
    real(r8),  intent(in)  :: q(:)  ! local radiosity vector
    real(r8),  intent(in)  :: t(:)
    real(r8),  intent(in)  :: tamb  ! ambient temperature
    real(r8),  intent(in)  :: eps(:)
    real(r8),  intent(out) :: f(:)  ! local heat flux vector

    integer :: i, j
    real(r8) :: s, global_q(this%nface_tot)

    ASSERT(size(q) == this%nface)
    ASSERT(size(t) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(size(f) == this%nface)

    call collate (global_q, q)
    call broadcast (global_q)

    do j = 1, this%nface
      s = 0.0_r8
      do i = this%ia(j), this%ia(j+1)-1
        s = s + this%vf(i) * global_q(this%ja(i))
      end do
      f(j) = eps(j)*this%sigma*(t(j)-this%t0)**4 - eps(j)*s
    end do

    if (allocated(this%amb_vf)) f = f - (this%sigma*(tamb-this%t0)**4)*eps*this%amb_vf

!    do j = 1, this%nface
!      s = q(j)
!      do i = this%ia(j), this%ia(j+1)-1
!        s = s - this%vf(i) * global_q(this%ja(i))
!      end do
!      f(j) = s
!    end do
!
!    if (allocated(this%amb_vf)) f = f - (this%sigma*(tamb-this%t0)**4) * this%amb_vf

  end subroutine flux

  !!
  !!  Radiosity system residual -- core computational kernel
  !!

  subroutine residual (this, q, t, tamb, eps, r)

    class(rad_system), intent(in) :: this
    real(r8), intent(in)  :: q(:)
    real(r8), intent(in)  :: t(:)
    real(r8), intent(in)  :: tamb
    real(r8), intent(in)  :: eps(:)
    real(r8), intent(out) :: r(:)

    integer :: i, j
    real(r8) :: s, global_q(this%nface_tot)

    ASSERT(size(q) == this%nface)
    ASSERT(size(t) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(size(r) == this%nface)

    call collate (global_q, q)
    call broadcast (global_q)

    do j = 1, size(r)
      s = 0.0_r8
      do i = this%ia(j), this%ia(j+1)-1
        s = s + this%vf(i) * global_q(this%ja(i))
      end do
      r(j) = eps(j)*this%sigma*(t(j)-this%t0)**4 - q(j) + (1.0_r8-eps(j))*s
    end do

    if (allocated(this%amb_vf)) r = r + (this%sigma*(tamb-this%t0)**4)*(1.0_r8-eps)*this%amb_vf

  end subroutine residual

  subroutine matvec1 (this, eps, z)

    class(rad_system), intent(in) :: this
    real(r8), intent(in) :: eps(:)
    real(r8), intent(inout) :: z(:)

    integer :: i, j
    real(r8) :: s, global_z(this%nface_tot)

    ASSERT(size(eps) == this%nface)
    ASSERT(size(z) == this%nface)

    call collate (global_z, z)
    call broadcast (global_z)

    do j = 1, this%nface
      s = 0.0_r8
      do i = this%ia(j), this%ia(j+1)-1
        s = s + this%vf(i) * global_z(this%ja(i))
      end do
      z(j) = eps(j)*s
    end do

  end subroutine matvec1

  !!
  !! Jacobian of the radiosity system RHS with respect to the temperatures.
  !! This is diagonal, and DRHS returns the diagonal.
  !!

  subroutine rhs_deriv (this, t, eps, drhs)

    class(rad_system), intent(in) :: this
    real(r8), intent(in) :: t(:)
    real(r8), intent(in) :: eps(:)
    real(r8), intent(out) :: drhs(:)

    ASSERT(size(t) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(size(drhs) == this%nface)

    drhs = 4*this%sigma*eps*(t-this%t0)**3

  end subroutine rhs_deriv

  subroutine rhs (this, t, tamb, eps, b)

    class(rad_system), intent(in)  :: this
    real(r8), intent(in)  :: tamb, t(:), eps(:)
    real(r8), intent(out) :: b(:)

    ASSERT(size(t) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(size(b) == this%nface)

    b = eps * this%sigma * (t-this%t0)**4
    if (allocated(this%amb_vf)) b = b + (this%sigma*(tamb-this%t0)**4)*(1.0_r8-eps)*this%amb_vf

  end subroutine rhs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CHEBY_SOLVE
 !!
 !! This subroutine solves the radiosity system using the Chebyshev iterative
 !! method.
 !!
 !! The Chebyshev iteration parameters c and d should define a family of
 !! ellipses in the right-hand plane with foci d±c that contain the spectrum
 !! of the system coefficient matrix.
 !!
 !! The radiosity system matrix has the form I - B.  The matrix B has real
 !! eigenvalues (a consequence of view factor reciprocity) contained in the
 !! interval (-1, 1) (B >= 0 with row sums < 1).  Thus in this case then the
 !! spectrum of I - B is contained in (0,2) and for optimal convergence the
 !! ellipse foci are set to the minimum and maximum eigenvalues:
 !!
 !!               d = (lmax + lmin)/2,   c = (lmax - lmin)/2.
 !!
 !! The routines lambda_min and lambda_max can be used to compute estimates
 !! of these extreme eigenvalues.  The minimum eigenvalue is the most important
 !! to obtain and fortunately this seems to be easily calculated.  The maximum
 !! eigenvalue is difficult to find with the implemented power method due to
 !! poor separation from other eigenvalues, but it is not important to find it.
 !! Typically it seems to be close to 1, and using this value for lmax seems
 !! to be acceptable even when the actual maximum eigenvalue is larger.
 !! Using the known bound of 2 is safest but cuts the convergence rate by
 !! up to half.
 !!
 !! IMPLEMENTATION NOTES
 !!
 !! 1) In "Adaptive procedure for estimating parameters for the nonsymmetric
 !!     Tchebychev iteration", Numer. Math. 31, 183-208 (1978), Manteuffel
 !!     describes procedures for dynamically determining the c and d
 !!     parameters.  It would be great to implement one of these procedures
 !!     and eliminate the need to specify c and d.
 !!

  subroutine cheby_solve (this, t, tamb, eps, c, d, tol, maxitr, q, numitr, error)

    class(rad_system), intent(in) :: this
    real(r8), intent(in) :: t(:)    ! face temperatures
    real(r8), intent(in) :: tamb    ! ambient temperature
    real(r8), intent(in) :: eps(:)  ! face emissivities
    real(r8), intent(in) :: c, d    ! chebyshev iteration parameters
    real(r8), intent(in) :: tol     ! convergence tolerance
    integer,  intent(in) :: maxitr
    real(r8), intent(inout) :: q(:) ! face radiosities
    integer,  intent(out)   :: numitr
    real(r8), intent(out)   :: error

    integer :: i, j
    real(r8) :: rhs_norm, res_norm, s, psi, omega
    real(r8), dimension(this%nface) :: rhs, r, p
    real(r8) :: global_q(this%nface_tot)

    ASSERT(size(q) == this%nface)
    ASSERT(size(t) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(maxitr > 0)
    ASSERT(tol > 0.0_r8)
    ASSERT(d > 0.0_r8)
    ASSERT(d - abs(c) > 0.0_r8)

    rhs = eps*this%sigma*(t-this%t0)**4 + &
          this%sigma*(tamb-this%t0)**4*(1.0_r8-eps)*this%amb_vf

    rhs_norm = global_l2norm(rhs)

    !! Handle the pathological case of a zero rhs.
    if (rhs_norm == 0.0_r8) then
      q = 0.0_r8
      numitr = 0
      error  = 0.0_r8
      return
    end if

    numitr = 0
    do

      !! Compute the residual of the current iterate.
      call collate (global_q, q)
      call broadcast (global_q)
      do j = 1, size(r)
        s = 0.0_r8
        do i = this%ia(j), this%ia(j+1)-1
          s = s + this%vf(i)*global_q(this%ja(i))
        end do
        r(j) = rhs(j) - q(j) + (1.0_r8-eps(j))*s
      end do

      if (numitr > 0) then

        !! Test residual of the last iterate.
        res_norm = global_l2norm(r)
        if (res_norm <= tol*rhs_norm) exit

        numitr = numitr + 1
        if (numitr > maxitr) exit

        if (numitr == 2) then ! starting values for the update parameters
          psi = 0.5_r8*(c/d)**2
          omega = 2*d/(2*d**2-c**2)
        else  ! recursion for update parameters
          psi = (c*omega/2.0_r8)**2
          omega = 1.0_r8/(d - omega*(c/2)**2)
        end if

        !! Next iterate.
        p = r + psi*p
        q = q + omega*p

      else  ! first iterate is a special case

        numitr = numitr + 1

        p = r
        q = q + (1.0_r8/d)*p

      end if

    end do

    error = res_norm/rhs_norm

  end subroutine cheby_solve

  subroutine cheby_precon (this, eps, c, d, numitr, q)

    class(rad_system), intent(in) :: this
    real(r8), intent(in) :: eps(:)  ! face emissivities
    real(r8), intent(in) :: c, d    ! chebyshev iteration parameters
    integer,  intent(in) :: numitr  ! number of iterations
    real(r8), intent(inout) :: q(:)  ! RHS on entry, approx solution on return

    integer :: i, j, n
    real(r8) :: s, psi, omega
    real(r8), dimension(this%nface) :: r, p, rhs
    real(r8) :: global_q(this%nface_tot)

    ASSERT(size(q) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(numitr > 0)
    ASSERT(d > 0.0_r8)
    ASSERT(d - abs(c) > 0.0_r8)

    rhs = q
    p = q
    q = (1.0_r8/d)*p

    do n = 2, numitr

      !! Compute the residual of the current iterate.
      call collate (global_q, q)
      call broadcast (global_q)
      do j = 1, size(r)
        s = 0.0_r8
        do i = this%ia(j), this%ia(j+1)-1
          s = s + this%vf(i)*global_q(this%ja(i))
        end do
        r(j) = rhs(j) - q(j) + (1.0_r8-eps(j))*s
      end do

      if (n == 2) then ! starting values for the update parameters
        psi = 0.5_r8*(c/d)**2
        omega = 2*d/(2*d**2-c**2)
      else  ! recursion for update parameters
        psi = (c*omega/2.0_r8)**2
        omega = 1.0_r8/(d - omega*(c/2)**2)
      end if

      !! Next iterate.
      p = r + psi*p
      q = q + omega*p

    end do

  end subroutine cheby_precon

  subroutine jacobi_precon (this, eps, numitr, z)

    class(rad_system), intent(in) :: this
    real(r8), intent(in) :: eps(:)  ! face emissivities
    integer,  intent(in) :: numitr  ! number of iterations
    real(r8), intent(inout) :: z(:)  ! RHS on entry, approx solution on return

    integer :: i, j, n
    real(r8) :: s, rhs(this%nface), global_z(this%nface_tot)

    ASSERT(size(z) == this%nface)
    ASSERT(size(eps) == this%nface)
    ASSERT(numitr > 0)

    if (numitr == 1) return

    rhs = z

    do n = 2, numitr
      call collate (global_z, z)
      call broadcast (global_z)
      do j = 1, this%nface
        s = 0.0_r8
        do i = this%ia(j), this%ia(j+1)-1
          s = s + this%vf(i)*global_z(this%ja(i))
        end do
        z(j) = rhs(j) + (1.0_r8-eps(j))*s
      end do
    end do

  end subroutine jacobi_precon

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! LAMBDA_MIN / LAMBDA_MAX
 !!
 !! These subroutines calculate the minimum and maximum eigenvalues of the
 !! radiosity system matrix using the power method. The matrix has the form
 !! I - B, where B >= 0 with real eigenvalues contained in the open interval
 !! (-1,1).  The minimum eigenvalue is obtained by applying the power method
 !! to the matrix I - (I - B) = B.  Under mild conditions, Perron-Frobenius
 !! guarantees a simple real eigenvalue equal to the spectral radius of B.
 !! In practice, the power method converges quickly and reliably to obtain
 !! the max eigenvalue of B.  In contrast, the power method converges very
 !! slowly to the maximum eigenvalue of I-B due to the poor separation of
 !! the smallest eigenvalues of B.
 !!
 !! The input argument EPS specifies the emissivities of the enclosure faces.
 !! An initial, nonzero vector that contains a component of desired eigenvector
 !! must be provided in Q; it need not be normalized.  A random vector is a
 !! good choice.  The calculated normalized eigenvector and corresponding
 !! eigenvalue are returned in Q and LMIN or LMAX, and satisfy
 !!
 !!     || (I-B)*Q - L*Q|| < TOL*L
 !!
 !! NB (May 2011): The above analysis for finding the minimum eigenvalue --
 !! the maximum eigenvalue of the non-negative B -- is not quite correct.
 !! Having a symmetric structure, B is effectively irreducible (it may
 !! upon permutation be block diagonal with irreducible blocks, but then
 !! the eigenvalue problem completely decouples into separate subproblems).
 !! But it may be cyclic and not a primitive matrix, having multiple
 !! eigenvalues with modulus equal to the spectral radius, and in this
 !! case the power method will not converge.  This is easily fixed, though,
 !! by applying the power method to a shifted B, beta*I + B, with beta>0
 !! since the shifted matrix is primitive.  See Wood and O'Neill,
 !! ANZIAM J. 45(E) 2007; also ANZIAM J. 48 2007.
 !! TODO: Those papers suggest that the method of Collatz would be a
 !! much more robust method of strictly bounding the maximum eigenvalue
 !! which is precisely what we want to do.
 !!

  subroutine lambda_min (this, eps, tol, maxitr, q, lmin, numitr, error)

    class(rad_system), intent(in) :: this
    real(r8), intent(in)    :: eps(:)
    real(r8), intent(in)    :: tol
    integer,  intent(in)    :: maxitr
    real(r8), intent(inout) :: q(:)
    real(r8), intent(out)   :: lmin
    integer,  intent(out)   :: numitr
    real(r8), intent(out)   :: error

    integer :: i, j, n
    real(r8) :: s, theta, shift, global_q(this%nface_tot), p(this%nface)

    ASSERT(size(eps) == this%nface)
    ASSERT(size(eps) == size(q))
    ASSERT(tol > 0.0_r8)
    ASSERT(maxitr > 0)

    p = q
    
    shift = 0.1_r8  ! anything > 0 will do (but effect the convergence rate)

    do n = 0, maxitr
      !! Replicate the global eigenvector iterate across all processes.
      s = global_l2norm(p)
      INSIST( s > 0.0_r8 )
      q = p / s
      call collate (global_q, q)
      call broadcast (global_q)

      !! Matrix-vector product.
      do j = 1, this%nface
        !s = 0.0_r8
        s = shift * q(j)
        do i = this%ia(j), this%ia(j+1)-1
          s = s + (1.0_r8-eps(j))*this%vf(i)*global_q(this%ja(i))
        end do
        p(j) = s
      end do

      !! Eigenvalue estimate.
      theta = global_dot_product(q, p)

      !! Convergence test.
      error = global_l2norm(p - theta*q)
      if (error < tol*abs(1+shift-theta)) exit
    end do

    lmin = 1 + shift - theta
    numitr = n

  end subroutine lambda_min

  subroutine lambda_max (this, eps, tol, maxitr, q, lmax, numitr, error)

    class(rad_system), intent(in) :: this
    real(r8), intent(in)    :: eps(:)
    real(r8), intent(in)    :: tol
    integer,  intent(in)    :: maxitr
    real(r8), intent(inout) :: q(:)
    real(r8), intent(out)   :: lmax
    integer,  intent(out)   :: numitr
    real(r8), intent(out)   :: error

    integer :: i, j, n
    real(r8) :: s, theta, global_q(this%nface_tot), p(this%nface)

    ASSERT(size(eps) == this%nface)
    ASSERT(size(eps) == size(q))
    ASSERT(tol > 0.0_r8)
    ASSERT(maxitr > 0)

    p = q

    do n = 0, maxitr
      !! Replicate the global eigenvector iterate across all processes.
      s = global_l2norm(p)
      INSIST( s > 0.0_r8 )
      q = p / s
      call collate (global_q, q)
      call broadcast (global_q)

      !! Matrix-vector product.
      do j = 1, this%nface
        s = q(j)
        do i = this%ia(j), this%ia(j+1)-1
          s = s - (1.0_r8-eps(j))*this%vf(i)*global_q(this%ja(i))
        end do
        p(j) = s
      end do

      !! Eigenvalue estimate.
      theta = global_dot_product(q, p)

      !! Convergence test.
      error = global_l2norm(p - theta*q)
      if (error < tol*abs(theta)) exit
    end do

    lmax = theta
    numitr = n

  end subroutine lambda_max

  function global_l2norm (x) result (l2norm)
    real(r8), intent(in) :: x(:)
    real(r8) :: l2norm, s
    integer :: j
    s = 0.0_r8
    do j = 1, size(x)
      s = s + x(j)**2
    end do
    l2norm = sqrt(global_sum(s))
  end function global_l2norm

end module rad_system_type
