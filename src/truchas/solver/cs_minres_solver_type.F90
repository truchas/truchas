!!
!! CS_MINRES_SOLVER_TYPE
!!
!! This module defines a linear solver for complex symmetric linear systems
!! that uses the CS-MINRES Krylov subspace method [1].
!!
!! This implementation is intended for consistent systems with use of a real
!! symmetric positive definite preconditioner. It does not go to the extra
!! effort of computing the QLP factorization from the QR factorization and
!! thus is not as well suited to very badly conditioned systems.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! [1] Sou-Cheng Choi, "Minimal Residual Methods for Complex Symmetric,
!!     Skew Symmetric, and Skew Hermitian Systems", arXiv:1304.6782, (2014)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cs_minres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op_class
  use parameter_list_type
  use mpi
  implicit none
  private

  type, public :: cs_minres_solver
    integer :: comm = MPI_COMM_WORLD, nproc, rank, lun
    integer :: max_iter
    real(r8) :: rel_tol
    logical :: verbose
    ! solver info
    integer :: num_iter
    real(r8) :: rel_rnorm, rnorm, Anorm, xnorm
  contains
    procedure :: init
    procedure :: solve
  end type

contains

  subroutine init(this, params, comm)

    use,intrinsic :: iso_fortran_env, only: output_unit

    class(cs_minres_solver), intent(out) :: this
    type(parameter_list) :: params
    integer, intent(in), optional :: comm

    integer :: ierr

    if (present(comm)) this%comm = comm
    call MPI_Comm_rank(this%comm, this%rank, ierr)
    call MPI_Comm_size(this%comm, this%nproc, ierr)

    call params%get('rel-tol', this%rel_tol, default=epsilon(1.0_r8))
    call params%get('max-iter', this%max_iter, default=10000)
    call params%get('verbose', this%verbose, default=.true.)
    this%verbose = this%verbose .and. (this%rank == 0)
    call params%get('output-unit', this%lun, default=output_unit)

  end subroutine


  subroutine solve(this, lin_op, b, x, stat, msg)

    use zvector_class

    class(cs_minres_solver), intent(inout) :: this
    class(complex_lin_op), intent(inout) :: lin_op
    class(zvector), intent(in) :: b
    class(zvector), intent(inout) :: x
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: msg

    integer :: k, ofreq
    real(r8) :: rel_rnorm, rnorm, Anorm, xnorm
    real(r8) :: cs, beta_1, beta_k, beta_km1, beta_kp1
    complex(r8) :: sn, alpha_k, gamma_k, delta_k, delta_kp1, dtmp, eps_k, eps_kp1, phi_k, tau_k
    class(zvector), allocatable :: q, vbar, z_kp1, z_k, z_km1
    class(zvector), allocatable :: w_k, w_km1, w_km2, tmp

    if (this%verbose) write(this%lun,'(a,es8.2,a,i0,a)') &
        'CS-MINRES(rel-tol=', this%rel_tol, ',max-iter=', this%max_iter, ')'

    !! Reset solver info
    this%num_iter = 0
    this%rel_rnorm = 0
    this%rnorm = 0
    this%Anorm = 0
    this%xnorm = 0

    ofreq = 1  ! starting output frequency when verbose

    call x%clone(q)
    call x%clone(vbar)
    call x%clone(z_k)
    call x%clone(z_kp1)

    !! Compute the initial z_1 and q_1 vectors, and beta_1 (k=0)
    call z_kp1%copy(b)
    call lin_op%precon(z_kp1, q)
    beta_1 = q%dotc(z_kp1) ! dotc value should be real
    if (beta_1 < 0.0_r8) then
      stat = -2
      msg  = 'preconditioner is not positive definite'
      if (this%verbose) write(this%lun,'(a,i0,": ",a)') 'stat=', stat, msg
      return
    end if
    beta_1 = sqrt(beta_1)

    call x%setval(0.0_r8)

    if (beta_1 == 0) then ! trivial solution
      stat = 0
      msg  = 'b=0 and exact solution is x=0'
      if (this%verbose) write(this%lun,'(a,i0,": ",a)') 'stat=', stat, msg
      return
    end if

    !! Initialize other quantities (k=0)
    beta_k = 0
    beta_kp1 = beta_1
    phi_k = beta_1

    rnorm = beta_1
    Anorm = 0.0_r8
    xnorm = 0.0_r8
    rel_rnorm = rnorm/(Anorm*xnorm + beta_1)

    if (this%verbose) write(this%lun,'(a)') '   iter  rel rnorm      rnorm      Anorm      xnorm'

    do k = 1, this%max_iter

      !! Lanczos process
      beta_km1 = beta_k; beta_k = beta_kp1 ! shift
      call move_alloc(z_k, z_km1); call move_alloc(z_kp1, z_k) ! shift
      call vbar%conjg(q) ! done with q_k, q now used as temporary
      call vbar%scale(cmplx(1/beta_k,kind=r8))
      call lin_op%matvec(vbar, q) ! this and next line forms p_k in q
      if (k > 1) call q%update(cmplx(-beta_k/beta_km1,kind=r8), z_km1) ! done with z_{k-1}
      alpha_k = vbar%dotu(q)
      call q%update(-alpha_k/beta_k, z_k) ! form z_{k+1} in q
      call move_alloc(q, z_kp1)
      call move_alloc(z_km1, q) ! recycle z_{k-1) storage for q_{k+1}
      call lin_op%precon(z_kp1, q) ! computes q_{k+1}
      beta_kp1 = z_kp1%dotc(q) ! dotc value should be real
      if (beta_kp1 < 0.0_r8) then
        stat = -2
        msg  = 'preconditioner is not positive definite'
        exit
      end if
      beta_kp1 = sqrt(beta_kp1)

      !! Apply previous Q_{k-1,k} to column k, k+1
      if (k == 1) then ! no previous Q; set initial gamma_1 and delta_2
        gamma_k = alpha_k
        delta_kp1 = beta_kp1
      else
        delta_k = delta_kp1 ! shift
        dtmp = cs*delta_k + sn*alpha_k
        gamma_k = conjg(sn)*delta_k - cs*alpha_k
        delta_k = dtmp
        if (k > 2) eps_k = eps_kp1 ! shift
        eps_kp1 = sn*beta_kp1
        delta_kp1 = -cs*beta_kp1
      end if

      !! Compute current Q_{k,k+1} and apply to column k and RHS
      call SymOrtho((gamma_k), cmplx(beta_kp1,kind=r8), cs, sn, gamma_k)
      tau_k = cs*phi_k
      phi_k = conjg(sn)*phi_k

      !! Next basis vector w_k
      select case (k)
      case (1)
        call move_alloc(vbar, w_k)
        call x%clone(vbar)
      case (2)
        call move_alloc(w_k, w_km1)
        call move_alloc(vbar, w_k)
        call w_k%update(-delta_k, w_km1)
        call x%clone(vbar)
      case (3)
        call move_alloc(w_km1, w_km2)
        call move_alloc(w_k, w_km1)
        call move_alloc(vbar, w_k)
        call w_k%update(-delta_k, w_km1, -eps_k, w_km2)
        call x%clone(vbar)
        call move_alloc(w_km2, tmp)
      case (4:)
        call move_alloc(w_km1, w_km2)
        call move_alloc(w_k, w_km1)
        call move_alloc(vbar, w_k)
        call w_k%update(-delta_k, w_km1, -eps_k, w_km2)
        call move_alloc(tmp, vbar)
        call move_alloc(w_km2, tmp)
      end select
      call w_k%scale(1/gamma_k)

      !! Next solution iterate x_k
      call x%update(tau_k, w_k)

      !! Output norms periodically. Lagged by one iteration in case of early
      !! termination of next iteration; final iteration output after exit.
      if (this%verbose .and. mod(k-1,ofreq) == 0) then
        if (any(k-1 == [100, 1000])) ofreq = 10*ofreq
        write(this%lun,'(i7,*(es11.2))') k-1, rel_rnorm, rnorm, Anorm, xnorm
      end if

      !! Update norms
      rnorm = abs(phi_k)
      if (k == 1) then
        Anorm = sqrt(beta_kp1**2 + abs(alpha_k)**2) ! Eqn B.2
      else
        Anorm = max(Anorm, sqrt(beta_kp1**2 + abs(alpha_k)**2 + beta_k**2)) ! Eqn B.2
      end if
      xnorm = x%norm2() ! no cheap estimate available without QLP
      rel_rnorm = rnorm/(Anorm*xnorm + beta_1) ! NRBE condition from Table B.1

      if (rel_rnorm <= this%rel_tol) then ! converged
        stat = 1
        msg  = 'solution converged within tolerance'
        exit
      end if

    end do

    if (k > this%max_iter) then
      stat = -1
      msg  = 'maximum number of iterations limit'
    end if

    !! Output info from final iteration
    if (this%verbose) then
      if (stat < 0) k = k - 1 ! last successful iteration when failed
      write(this%lun,'(i7,*(es11.2))') k, rel_rnorm, rnorm, Anorm, xnorm
      write(this%lun,'(a,i0,": ",a)') 'stat=', stat, msg
    end if

    !! Store solver info
    this%num_iter = k
    this%rel_rnorm = rel_rnorm
    this%rnorm = rnorm
    this%Anorm = Anorm
    this%xnorm = xnorm

  end subroutine solve


  pure subroutine SymOrtho(a, b, c, s, r)

    complex(r8), intent(in) :: a, b
    real(r8), intent(out) :: c
    complex(r8), intent(out) :: s, r

    real(r8) :: absa, absb, t

    absa = abs(a)
    absb = abs(b)

    if (absb == 0) then
      c = 1.0_r8
      s = 0.0_r8
      r = a
    else if (absa == 0) then
      c = 0.0_r8
      s = 1.0_r8
      r = b
    else if (absb >= absa) then
      t = absa/absb
      c = 1/sqrt(1+t*t) ! temporary
      s = c*conjg(b/absb)*(a/absa)
      c = c*t
      r = b/conjg(s)
    else
      t = absb/absa
      c = 1/sqrt(1+t*t)
      s = c*t*conjg(b/absb)*(a/absa)
      r = a/c
    end if

  end subroutine SymOrtho

end module cs_minres_solver_type
