! -*- Mode: f90 -*-

!   anorm - appropriate norm of coefficient matrix
!       NOTE: used only when ABS(stopping) = stopping_axb.
!
!   bnorm - appropriate norm of source vector
!       NOTE: used only when ABS(stopping) = stopping_axb or stopping_b.
!
!   eps - convergence criterion [DEFAULT = 1e-7]
!
!   iter - number of iterations completed
!       NOTE: if iter = 0 and status = 0, the initial guess was found to be the
!             solution, so no iterations were required
!
!   itmax - maximum number of iterations to allow [NO DEFAULT]
!       NOTE: itmax must be less than or equal to the
!             dimensions of rnorm, err, rnormt, and errt
!
!   luout - logical unit to which output should be written
!           [DEFAULT is 6, i.e. stdout]
!
!   luerr - logical unit to which error messages should be
!           written [DEFAULT is 7 or 0, i.e. stderr]
!
!   lutty - logical unit to which intermediate status messages
!           should be written [DEFAULT is 6, i.e. stdout]
!
!   method - indicates which method is being used
!
!            use named constants set in constants_module:
!                  method_LU
!                  method_CG
!                  method_GMRES
!                  method_FGMRES
!                  method_BCGS
!                  method_TFQMR
!                  method_SOR
!                  method_SSOR
!                  method_Jacobi
!                  method_IC
!                  method_ILU
!
!   NOTE: There's really not much point in using BCGS, except for
!         research purposes, since TFQMR is safer and not much
!         more expensive computationally.
!
!   subspace - number of old vectors to use in truncated methods
!          like GMRES (not required for non-truncated methods)
!          [DEFAULT = 20]
!
!   NOTE: Typically, the number of iterations (and the likelihood of
!         failure to converge) decreases as subspace increases.
!         However, the cost of the algorithm, both in operations and
!         memory required, increases as subspace increases.
!   
!   norm_type - determines which norm to use
!      0  ==>  infinity norm
!      1  ==>  1-norm
!      2  ==>  2-norm for vectors, Frobenius norm for matrices
!              [DEFAULT]
!
!   NOTES:
!    - For a vector x of length n:
!                                  
!        infinity norm  ==>   max   | x |
!                           1<=i<=n    i
!
!                     ---
!                     \  
!        1-norm  ==>  /   | x |
!                     ---    i
!                      n
!
!                          ---
!                          \    2          T
!        2-norm  ==>  sqrt /   x  = sqrt (x x)
!                          ---  i
!                           n
!
!    - For an n-by-m matrix A:
!                                     m
!                                    ---
!                                    \  
!        infinity norm  ==>    max   /   | a  |
!                            1<=i<=n ---    ij
!                                    j=1
!                              n 
!                             ---
!                             \  
!        1-norm  ==>    max   /   | a  |
!                     1<=j<=m ---    ij
!                             i=1
!                                   n   m
!                                  --- ---
!                                  \   \    2
!        Frobenius norm  ==>  sqrt /   /   a
!                                  --- ---  ij
!                                  i=1 j=1
!
!     - Note that the matrix 1-norm is expensive to calculate for
!       matrices stored in some formats (e.g. ELL).
!
!   omega - relaxation parameter [DEFAULT = 1.0]
!
!   outmode - output mode:
!      0  ==>  no output [DEFAULT]
!      1  ==>  write out error messages only
!      2  ==>  same as iout = 1, plus write out warning messages
!      3  ==>  same as iout = 2, plus write a one-line summary consisting of
!              the iteration number, the norms of the calculated residual and
!              error estimate, and the norms of the true residual and error
!              estimate (if calculated) each time convergence is checked
!      4  ==>  same as iout = 3, plus write out the coefficient,
!              preconditioner (if there is one), source, initial guess and
!              converged solution
!      5  ==>  same as iout = 4, plus write out the current iterate and true
!              residual (if computed) at each iteration
!      6  ==>  same as iout = 5, plus write out vectors at each stage of the
!              iteration, details about the preconditioner if available, etc.
!              (WARNING: this generates reams of output)
!
!   output_frequency - frequency with which to print a msg. as to the
!       progress of the iteration
!
!   residual - determines how often to compute and/or use true
!            residual (b-Ax) and how often to compute error estimate
!            based on true residual (errt)
!      0  ==>  never [DEFAULT]
!      1  ==>  compute b-Ax and errt each time a new iterate is
!              obtained, but never replace recursively-generated
!              residual (r) with b-Ax
!      2  ==>  replace r with b-Ax every MAX(10,SQRT(n)) iterations
!      3  ==>  compute b-Ax every iteration, but only replace r
!              with b-Ax every MAX(10,SQRT(n)) iterations
!
!   NOTES:
!    - Residual = 1 and residual = 3 are mainly for investigating
!      the behavior of solvers, since time is wasted calculating
!      but not using b-Ax.
!
!    - Use of residual > 1 can be very beneficial to CG, and is
!      recommended.  For production use and when CPU time is a
!      consideration, use residual = 2.
!
!    - Values of residual > 1 are ignored by GMRES, since b-Ax
!      is recalculated automatically each subspace iterations.
!
!    - Values of residual > 1 are meaningless for stationary
!      methods (Jacobi, SOR, SSOR), since they do not generate
!      a residual during the course of the algorithm.
!
!    - Values of residual > 1 are ignored by BCGS and TFQMR.
!      It doesn't seem to benefit these methods, and seems to harm
!      both in some situations.
!
!   save_locals - whether locally-allocated arrays should be deallocated
!                 prior to returning or saved for future use
!
!   steps - number of Jacobi or SSOR steps to take when
!           precond = 3 or 4 [DEFAULT = 1]
!
!   stopping - determines which stopping test to use
!          (if negative, then the true residual will be used)
!      0  ==>  || x - xold || / || x ||          [DEFAULT FOR STATIONARY METHODS]
!      1  ==>  || r || / ( || A ||*|| x || + || b || )
!      2  ==>  || r || / || b ||              [DEFAULT FOR NONSTATIONARY METHODS]
!      3  ==>  || r || / || x ||
!      4  ==>  || r || / || r0 ||
!      5  ==>  || r ||
!
!   where A, x, and b are the coefficient matrix, current estimate
!   of the solution vector, and source vector of the linear system
!   being solved, Ax = b.  The vector r is the residual, r = b - Ax,
!   and r0 is the initial residual (the residual computed using the
!   initial guess for x).
!
!   Recommendations: CG     -> stopping = 2
!                    GMRES  -> stopping = 2
!                    BCGS   -> stopping = -2
!                    TFQMR  -> stopping = -2
!                    Jacobi -> stopping = 0
!                    SOR -> stopping = 0
!                    SSOR -> stopping = 0
!
!   NOTES:
!    - The norm to use is set by norm.
!
!    - Stopping test 0 is fine for the stationary methods (Jacobi,
!      SOR, SSOR) since 1) the residual is not directly available,
!      and 2) these methods converge nearly monotonically anyway.
!      However, this test is completely inappropriate and dangerous
!      for the nonstationary methods (CG, GMRES, BCGS, TFQMR).
!      Therefore it is not allowed for these methods.
!
!      Note that this is a change from previous versions of UbikSolve,
!      which in fact used stopping test 0 as the default since it
!      was the only test that worked for all methods.  In retrospect
!      that was a bad idea.
!
!    - Stopping test 2 is probably the best for most situations.
!      However, it can lead to difficulties when
!                  || A ||*|| x || >> || b ||
!      In that case, use stopping test 1.
!
!    - Stopping tests 1 and 3 increase the cost of GMRES
!      significantly, since a new iterate is normally not computed
!      at each iteration.
!
!    - Note that if stopping test 1 is used, the norm of the
!      coefficient is needed.  This value is calculated unless the
!      storage format is unknown (istore = 0), in which case
!      the value must be passed in as anorm.
!
!    - Note that when the initial guess of x = 0 is used, test 4 is
!      equivalent to test 2.
!
!    - Stopping test 5 should be used carefully (if ever).
!
!    - Stopping test 2, based on the true residual, is recommended
!      for BCGS and TFQMR simply to protect against false
!      convergence.  However, it adds significant overhead to each
!      iteration.
!
!      A better (more efficient) way to protect against false
!      convergence is to use a stopping test based on the
!      recursively-generated residual (stopping > 0), and then compute
!      the true residual and corresponding error estimate after
!      'convergence'.  If the error estimate based on the true
!      residual does not meet the desired convergence criterion, the
!      the final iterate can be used as the initial guess for
!      additional iterations with the same (or a different) method.
!
!    - For the stationary iterative methods (Jacobi, SOR, SSOR), the
!      true residual is always used when stopping /= 0, since an
!      estimated residual isn't built up during the iteration.  So,
!      e.g., stopping = 1 is the same as stopping = -1.  This adds quite a
!      bit of expense to each iteration, but is useful if you want
!      to plot the residual history.
!
!   small - small floating-point value
!          [DEFAULT = EPSILON(1.0d0)]
!
!   smallest - smallest floating-point value to worry about
!          [DEFAULT = TINY(1.0)]
!
!   for reference, on the x86 architecture:
!
!        PRECISION(1.0)   =  6
!        DIGITS(1.0)      =  24
!        EPSILON(1.0)     =  1.19E-07
!        TINY(1.0)        =  1.17E-38
!        HUGE(1.0)        =  3.40E+38
!        
!        PRECISION(1.0d0) =  15
!        DIGITS(1.0d0)    =  53
!        EPSILON(1.0d0)   =  2.22E-016
!        TINY(1.0d0)      =  2.22E-308
!        HUGE(1.0d0)      =  1.80E+308
!
!   NOTES:
!    - If the initial residual norm is smaller than this value,
!      convergence is assumed.
!
!--------------------------------------------------------------------------
!
! Since preconditioning is no longer chosen by a member of control_type,
! the following info will eventually be moved to documentation.
!
!   Preconditioning recommendations:
!
!     - M-step Jacobi, with steps > 1, vectorizes well, but typically
!       is not very effective.
!
!     - M-step SSOR is often very effective, but is difficult to
!       vectorize (the implementation in UbikSolve is not vectorized).
!
!     - For the incomplete factorization preconditioners, omega
!       controls the degree to which fill-in elements are added back
!       into the diagonal:
!         omega = 0  ==>  standard incomplete factorization
!         omega = 1  ==>  modified incomplete factorization
!         0<omega<1  ==>  relaxed incomplete factorization
!
!     - The incomplete factorization preconditioners can be quite
!       effective in reducing the number of iterations, but are
!       often so expensive to compute as to negate any benefit in
!       using them.  In addition both the generation and application
!       of incomplete factorization preconditioners is difficult to
!       vectorize, and the implementations in UbikSolve are not 
!       vectorized.  One way to make them more useful when solving a
!       number of linear systems with similar coefficients (e.g. a 
!       time-dependent problem) is to amortize the cost of computing
!       them over several systems.  That is, call the solver routine
!       several times before updating the factorization.
!
!     - While incomplete Cholesky (IC) factorization preconditioners
!       are appropriate only for symmetric coefficients, incomplete
!       LU (ILU) factorization preconditioners can be used for
!       either symmetric or nonsymmetric coefficients.  Typically,
!       for a symmetric coefficient, one of the IC factorizations
!       should be preferrable to one of the ILU factorizations, since
!       the number of operations required to generate an IC
!       factorization would be about 1/2 that required for an ILU
!       factorization.  However, if the coefficient is stored in ELL
!       format, the nature of the data structure causes generation
!       of an IC factorization to be highly inefficient.  Hence, for
!       coefficients stored in ELL format, ILU factorization
!       preconditioners are preferrable to IC factorization
!       preconditioners, even if the coefficient is symmetric.
!
!--------------------------------------------------------------------------
module ubik_control_type_module
  use ubik_types_module
  implicit none
  save

  interface create
     module procedure create_control
  end interface

  interface destroy
     module procedure destroy_control
  end interface

  interface Ubik_nullify
     module procedure nullify_control
  end interface

  interface set_status_alloc_failure
     module procedure set_status_alloc_failure_cntrl
  end interface

  interface alloc_failure
     module procedure alloc_failure_control
  end interface

#ifndef dec
  interface set
     module procedure set_control
  end interface
#endif

  interface initialize
     module procedure initialize_control
  end interface

  interface cpu
     module procedure cpu_default
     module procedure cpu_iter
  end interface

  interface set_cpu
     module procedure set_cpu_default
     module procedure set_cpu_iter
  end interface

  interface rnorm
     module procedure rnorm_default
     module procedure rnorm_iter
  end interface

  interface set_rnorm
     module procedure set_rnorm_default
     module procedure set_rnorm_iter
  end interface

  interface err
     module procedure err_default
     module procedure err_iter
  end interface

  interface set_err
     module procedure set_err_default
     module procedure set_err_iter
  end interface

  interface rnormt
     module procedure rnormt_default
     module procedure rnormt_iter
  end interface

  interface set_rnormt
     module procedure set_rnormt_default
     module procedure set_rnormt_iter
  end interface

  interface errt
     module procedure errt_default
     module procedure errt_iter
  end interface

  interface set_errt
     module procedure set_errt_default
     module procedure set_errt_iter
  end interface

  ! internal named constants
  integer(int_type), parameter :: outmode_none__ = 0
  integer(int_type), parameter :: outmode_errors__ = 1
  integer(int_type), parameter :: outmode_warnings__ = 2
  integer(int_type), parameter :: outmode_summary__ = 3
  integer(int_type), parameter :: outmode_system__ = 4
  integer(int_type), parameter :: outmode_iterates__ = 5
  integer(int_type), parameter :: outmode_full__ = 6

  integer(int_type), parameter :: residual_default__ = 0
  integer(int_type), parameter :: residual_compute_only__ = 1
  integer(int_type), parameter :: residual_update__ = 2
  integer(int_type), parameter :: residual_compute_update__ = 3

  integer(int_type), parameter :: stdout__ = 6
#ifdef hpux
  integer(int_type), parameter :: stderr__ = 7
#else
  integer(int_type), parameter :: stderr__ = 0
#endif

  ! these should be internal as well, but I needed them outside
  integer(int_type), parameter :: stopping_relchg__ = 0
  integer(int_type), parameter :: stopping_axb__ = 1
  integer(int_type), parameter :: stopping_b__ = 2
  integer(int_type), parameter :: stopping_x__ = 3
  integer(int_type), parameter :: stopping_r0__ = 4
  integer(int_type), parameter :: stopping_r__ = 5
  integer(int_type), parameter :: stopping_min__ = 0
  integer(int_type), parameter :: stopping_max__ = 5

  private
  public :: &
       
       ! derived type containing params
       control_type, &
       
       ! initialization methods
#ifdef dec
       set_control, &
#else
       set, &
#endif
       set_defaults, &
       initialize, &
       create, &
       Ubik_nullify, &
       destroy

  public :: &
       
       ! routines to set params
       set_itmax, &
       set_luout, &
       set_luerr, &
       set_lutty, &
       set_method, &
       set_subspace, set_nold, &
       set_norm_type, &
       set_outmode_none, &
       set_outmode_errors, &
       set_outmode_warnings, &
       set_outmode_summary, &
       set_outmode_system, &
       set_outmode_iterates, &
       set_outmode_full, &
       set_output_frequency, &
       set_residual, &
       set_residual_default, &
       set_residual_compute_only, &
       set_residual_update, &
       set_residual_compute_upd, &
       set_scope, &
       set_scope_global, &
       set_scope_local

  public :: &
       set_steps, &
       set_stopping_relchg, &
       set_stopping_axb, &
       set_stopping_b, &
       set_stopping_x, &
       set_stopping_r0, &
       set_stopping_r, &
       set_stopping_true, &
       set_stopping_recursive, &
       set_anorm, &
       set_bnorm, set_need_bnorm, &
       set_eps, &
       set_omega, &
       set_small, &
       set_smallest, &
       set_status, &
       set_cpu, &
       set_rnormt, &
       set_errt, &
       set_rnorm, &
       set_err, &
       set_xold, set_need_xold

  public :: &
       
       ! routines to access params
       freq, &
       iter, &
       itmax, &
       luout, &
       luerr, &
       lutty, &
       matvecs, &
       method, &
       preconds, &
       subspace, nold, &
       norm_type, &
       outmode_none, &
       outmode_errors, &
       outmode_warnings, &
       outmode_summary, &
       outmode_system, &
       outmode_iterates, &
       outmode_full, &
       output_frequency, &
       residual, &
       residual_default, &
       residual_compute_only, &
       residual_update, &
       residual_compute_update, &
       scope, &
       scope_global, &
       scope_local

  public :: &
       steps, &
       stopping_true, &
       stopping_relchg, &
       stopping_axb, &
       stopping_b, &
       stopping_x, &
       stopping_r0, &
       stopping_r, &
       anorm, &
       bnorm, &
       eps, &
       omega, &
       small, &
       smallest, &
       cpu, &
       rnormt, &
       errt, &
       rnorm, &
       err, &
       xold

  public :: &
       converged, &
       invalid_input, &
       alloc_failure, &
       internal_error, &
       itmax_exceeded, &
       breakdown, &
       matmul_error, &
       precond_error, &
       next_iter, &
       set_status_converged, &
       set_status_invalid_input, &
       set_status_alloc_failure, &
       set_status_internal_error, &
       set_status_itmax_exceeded, &
       set_status_breakdown, &
       set_status_matmul_error, &
       set_status_precond_error, &
       set_status_next_iter

  ! remainder are public to rest of UbikSolve, but not made public to host codes
  public :: &
       set_iter, increment_iter, &
       set_matvecs, increment_matvecs, &
       set_preconds, increment_preconds, &
       need_bnorm, &
       need_xold, &
       update_residual, &
       check_stopping_stationary, &
       check_stopping_nonstationary, &
       normalize_cpu

  type control_type
     private

     ! dimensioned to number of unknowns
     real(real_type), dimension(:), pointer :: xold_ => NULL()

     ! dimensioned to itmax
     real(real_type), dimension(:), pointer :: cpu_ => NULL()
     real(real_type), dimension(:), pointer :: rnormt_ => NULL()
     real(real_type), dimension(:), pointer :: errt_ => NULL()
     real(real_type), dimension(:), pointer :: rnorm_ => NULL()
     real(real_type), dimension(:), pointer :: err_ => NULL()

     real(real_type) :: anorm_
     real(real_type) :: bnorm_
     real(real_type) :: eps_
     real(real_type) :: omega_
     real(real_type) :: small_
     real(real_type) :: smallest_

     integer(int_type) :: status_
     integer(int_type) :: freq_
     integer(int_type) :: iter_
     integer(int_type) :: itmax_
     integer(int_type) :: luout_
     integer(int_type) :: luerr_
     integer(int_type) :: lutty_
     integer(int_type) :: matvecs_
     integer(int_type) :: method_
     integer(int_type) :: norm_type_
     integer(int_type) :: outmode_
     integer(int_type) :: output_frequency_
     integer(int_type) :: preconds_
     integer(int_type) :: subspace_

     integer(int_type) :: residual_
     integer(int_type) :: scope_
     integer(int_type) :: steps_
     integer(int_type) :: stopping_

     logical(logical_type) :: outmode_none_
     logical(logical_type) :: outmode_errors_
     logical(logical_type) :: outmode_warnings_
     logical(logical_type) :: outmode_summary_
     logical(logical_type) :: outmode_system_
     logical(logical_type) :: outmode_iterates_
     logical(logical_type) :: outmode_full_

     logical(logical_type) :: stopping_true_
     logical(logical_type) :: stopping_relchg_
     logical(logical_type) :: stopping_axb_
     logical(logical_type) :: stopping_b_
     logical(logical_type) :: stopping_x_
     logical(logical_type) :: stopping_r0_
     logical(logical_type) :: stopping_r_

     logical(logical_type) :: residual_default_
     logical(logical_type) :: residual_compute_only_
     logical(logical_type) :: residual_update_
     logical(logical_type) :: residual_compute_update_

     logical(logical_type) :: save_locals_
     logical(logical_type) :: need_bnorm_
     logical(logical_type) :: need_xold_

  end type control_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Private methods
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_outmode_bools (control)
    type(control_type), intent(inout) :: control

    control%outmode_none_     = control%outmode_ == outmode_none__
    control%outmode_errors_   = control%outmode_ >= outmode_errors__
    control%outmode_warnings_ = control%outmode_ >= outmode_warnings__
    control%outmode_iterates_ = control%outmode_ >= outmode_iterates__
    control%outmode_system_   = control%outmode_ >= outmode_system__
    control%outmode_summary_  = control%outmode_ >= outmode_summary__
    control%outmode_full_     = control%outmode_ >= outmode_full__

    return
  end subroutine set_outmode_bools

  subroutine set_stopping_bools (control)
    type(control_type), intent(inout) :: control

    control%stopping_relchg_ = control%stopping_ == stopping_relchg__
    control%stopping_axb_    = control%stopping_ == stopping_axb__
    control%stopping_b_      = control%stopping_ == stopping_b__
    control%stopping_r0_     = control%stopping_ == stopping_r0__
    control%stopping_r_      = control%stopping_ == stopping_r__
    control%stopping_x_      = control%stopping_ == stopping_x__

    return
  end subroutine set_stopping_bools

  subroutine set_residual_bools (control)
    type(control_type), intent(inout) :: control

    control%residual_default_        = control%residual_ == residual_default__
    control%residual_compute_only_   = control%residual_ == residual_compute_only__
    control%residual_update_         = control%residual_ == residual_update__
    control%residual_compute_update_ = control%residual_ == residual_compute_update__

    return
  end subroutine set_residual_bools

  function destroy_partial(control) result(status)
    type(control_type), intent(inout) :: control
    integer(int_type) :: status

    status = 0

    if (associated(control%errt_)) then
       deallocate( control%errt_, stat=status )
    end if
    if (associated(control%rnormt_)) then
       deallocate( control%rnormt_, stat=status )
    end if
    if (associated(control%err_)) then
       deallocate( control%err_, stat=status )
    end if
    if (associated(control%rnorm_)) then
       deallocate( control%rnorm_, stat=status )
    end if
    if (associated(control%cpu_)) then
       deallocate( control%cpu_, stat=status )
    end if

    return
  end function destroy_partial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! General public methods
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine create_control (control, m)
    use ubik_constants_module

    integer(int_type), intent(in), optional :: m
    type(control_type), intent(inout) :: control
    integer(int_type) :: status, n

    if (PRESENT(m)) then
       n = m
    else
       n = itmax(control)
    end if

    if (n <= 0) then
       control%status_ = invalid_input__
       print *, 'badness in control::create!!'
       return
    end if

    if (associated(control%cpu_)) then

       ! components are already allocated
       if (SIZE(control%cpu_) >= n+1) then

          ! sizes are adequate
          return
       else

          ! components are not sized correctly, and hence need to be dealloated
          ! and then reallocated
          status = destroy_partial(control)

          if (status /= 0) then
             call set_status_alloc_failure(control)
             return
          end if
       end if
    end if

    allocate( control%cpu_(0:n), stat=status )
    allocate( control%rnorm_(0:n), stat=status )
    allocate( control%err_(0:n), stat=status )
    allocate( control%rnormt_(0:n), stat=status )
    allocate( control%errt_(0:n), stat=status )

    if (status /= 0) then
       call set_status_alloc_failure(control)
    end if

    return
  end subroutine create_control

  subroutine nullify_control (control)
    type(control_type), intent(inout) :: control

    nullify( control%errt_ )
    nullify( control%rnormt_ )
    nullify( control%err_ )
    nullify( control%rnorm_ )
    nullify( control%cpu_ )
    nullify( control%xold_ )

    return
  end subroutine nullify_control

  subroutine destroy_control (control)
    type(control_type), intent(inout) :: control
    integer(int_type) :: status

    status = destroy_partial(control)
    if (status /= 0) then
       call set_status_alloc_failure(control)
       return
    end if

    if (associated(control%xold_)) then
       deallocate( control%xold_, stat=status )
       if (status /= 0) then
          call set_status_alloc_failure(control)
       end if
    end if

    return
  end subroutine destroy_control

  !*********************************************************************
  ! Sets reasonable default values for control parameters.
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Out:
  !   control - structure containing control parameters
  !*********************************************************************
  subroutine set_defaults (control)
    use Ubik_ParallelInfo_module, only: initialize_parallel, &
                                   scope_global__, scope_local__

    type(control_type), intent(inout) :: control

    ! don't save local arrays
    control%save_locals_ = .false.

    ! initialize parallel info
    call initialize_parallel ()

    ! use recursively-generated residual rather than b-Ax by default
    call set_stopping_recursive (control)

    call set_itmax (control, 500)
    call set_iter (control, 0)
    call set_freq (control, 0)
    call set_luout (control, stdout__)
    call set_luerr (control, stderr__)
    call set_lutty (control, stdout__)
    call set_matvecs (control, 0)
    call set_method (control, 99)
    call set_norm_type (control, 2)
    call set_output_frequency (control, 0)
    call set_preconds (control, 0)
    call set_residual (control, residual_default__)
    call set_scope (control, scope_global__)
    call set_steps (control, 1)
    call set_subspace (control, 20)
    call set_anorm (control, 0.0d0)
    call set_bnorm (control, 0.0d0)
    call set_eps (control, 1.0d-7)
    call set_omega (control, 1.0d0)
    call set_small (control, EPSILON(1.0d0))
    call set_smallest (control, DBLE(TINY(1.0)))

    call set_outmode_warnings (control)
    call set_stopping_b (control)

    return
  end subroutine set_defaults

  subroutine set_control (control, freq, iter, itmax, luout, luerr, lutty, &
       matvecs, method, norm_type, output_frequency, preconds, residual, &
       scope, steps, subspace, anorm, bnorm, eps, omega, small, smallest)

    integer(int_type), intent(in), optional :: freq, iter, itmax, luout, &
         luerr, lutty, matvecs, method, norm_type, output_frequency, preconds, &
         residual, scope, steps, subspace
    real(real_type), intent(in), optional :: anorm, bnorm, eps, omega, small, &
         smallest

    type(control_type), intent(inout) :: control

    if (PRESENT(itmax)) then
       call set_itmax (control, itmax)
    end if
    if (PRESENT(iter)) then
       call set_iter (control, iter)
    end if
    if (PRESENT(freq)) then
       call set_freq (control, freq)
    end if
    if (PRESENT(luout)) then
       call set_luout (control, luout)
    end if
    if (PRESENT(luerr)) then
       call set_luerr (control, luerr)
    end if
    if (PRESENT(lutty)) then
       call set_lutty (control, lutty)
    end if
    if (PRESENT(matvecs)) then
       call set_matvecs (control, matvecs)
    end if
    if (PRESENT(method)) then
       call set_method (control, method)
    end if
    if (PRESENT(norm_type)) then
       call set_norm_type (control, norm_type)
    end if
    if (PRESENT(output_frequency)) then
       call set_output_frequency (control, output_frequency)
    end if
    if (PRESENT(preconds)) then
       call set_preconds (control, preconds)
    end if
    if (PRESENT(residual)) then
       call set_residual (control, residual)
    end if
    if (PRESENT(scope)) then
       call set_scope (control, scope)
    end if
    if (PRESENT(steps)) then
       call set_steps (control, steps)
    end if
    if (PRESENT(subspace)) then
       call set_subspace (control, subspace)
    end if
    if (PRESENT(anorm)) then
       call set_anorm (control, anorm)
    end if
    if (PRESENT(bnorm)) then
       call set_bnorm (control, bnorm)
    end if
    if (PRESENT(eps)) then
       call set_eps (control, eps)
    end if
    if (PRESENT(omega)) then
       call set_omega (control, omega)
    end if
    if (PRESENT(small)) then
       call set_small (control, small)
    end if
    if (PRESENT(smallest)) then
       call set_smallest (control, smallest)
    end if

    return
  end subroutine set_control

  !*********************************************************************
  ! Initialize control parameters (for use at the start of a solve).
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Out:
  !   control - structure containing control parameters
  !*********************************************************************
  subroutine initialize_control (control)
    use Ubik_ParallelInfo_module, only: initialize_scope

    type(control_type), intent(inout) :: control

    ! obtain CPU time at start of solution
    call set_cpu (control, 0)

    ! local or global scope (this *must* be done before calling
    ! set_freq_default, since it performs a vector sum)
    call initialize_scope (control%scope_)

    ! set logical vars - shouldn't be necessary
    call set_outmode_bools (control)
    call set_stopping_bools (control)

    ! cpu time and residual and error norm arrays
    control%cpu_ = 0.0
    control%err_ = 0.0
    control%errt_ = 0.0
    control%rnorm_ = 0.0
    control%rnormt_ = 0.0

    ! freq at which to use true residual - minimum of every 10 iterations
    if (freq(control) == 0) then
       call set_freq_default (control)
    end if

    ! Initialize iteration counter.
    call set_iter (control, 0)

    ! determine if || b || will be needed
    if ((stopping_axb(control) .or. stopping_b(control)) .and. bnorm(control) <= 0.0) then
       call set_need_bnorm (control)
    else
       control%need_bnorm_ = .FALSE.
    end if

    ! determine if || b || will be needed
    if ((stopping_relchg(control))) then
       call set_need_xold (control)
    else
       control%need_xold_ = .FALSE.
    end if

    return
  end subroutine initialize_control

  ! update residual?
  subroutine set_freq_default (control)
    use Ubik_Vector_module

    type(control_type), intent(inout) :: control

    ! freq at which to use true residual - set minimum of every 10 iterations
    call set_freq (control, MAX(10, INT(SQRT(Ubik_Sum((/REAL(control%itmax_, real_type)/))))))

    return
  end subroutine set_freq_default

  function update_residual (control)
    type(control_type), intent(inout) :: control
    logical(logical_type) :: update_residual

    if (freq(control) == 0) then
       call set_freq_default (control)
    end if

    update_residual = &
         iter(control) == 0 .or. &
         ( &
           (MOD(iter(control),freq(control)) == 0) .and. &
           (residual_update(control) .or. residual_compute_update(control)) &
         )

    return
  end function update_residual

  function need_bnorm(control)
    logical(logical_type) :: need_bnorm
    type(control_type), intent(in) :: control

    need_bnorm = control%need_bnorm_
    return
  end function need_bnorm

  subroutine set_need_bnorm(control)
    type(control_type), intent(inout) :: control

    control%need_bnorm_ = .TRUE.
    return
  end subroutine set_need_bnorm

  function need_xold(control)
    logical(logical_type) :: need_xold
    type(control_type), intent(in) :: control

    need_xold = control%need_xold_
    return
  end function need_xold

  subroutine set_need_xold(control)
    type(control_type), intent(inout) :: control

    control%need_xold_ = .TRUE.
    return
  end subroutine set_need_xold

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Specific public methods related to each component of control_type
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! true residual update frequency
  function freq(control)
    integer(int_type) :: freq
    type(control_type), intent(in) :: control

    freq = control%freq_
    return
  end function freq

  subroutine set_freq (control, frq)
    integer(int_type), intent(in) :: frq
    type(control_type), intent(inout) :: control

    control%freq_ = frq
    return
  end subroutine set_freq

  ! iteration counter
  function iter(control)
    integer(int_type) :: iter
    type(control_type), intent(in) :: control

    iter = control%iter_
    return
  end function iter

  subroutine set_iter (control, it)
    integer(int_type), intent(in) :: it
    type(control_type), intent(inout) :: control

    control%iter_ = it
    return
  end subroutine set_iter

  subroutine increment_iter (control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%iter_ = control%iter_ + 1
    if (control%iter_ > control%itmax_) then
       control%iter_ = control%iter_ - 1
       control%status_ = itmax_exceeded__
    end if
    return
  end subroutine increment_iter

  ! max iters to allow
  function itmax(control)
    integer(int_type) :: itmax
    type(control_type), intent(in) :: control

    itmax = control%itmax_
    return
  end function itmax

  subroutine set_itmax (control, itmax)
    integer(int_type), intent(in) :: itmax
    type(control_type), intent(inout) :: control

    control%itmax_ = itmax
    return
  end subroutine set_itmax

  ! logical unit for normal output
  function luout(control)
    integer(int_type) :: luout
    type(control_type), intent(in) :: control

    luout = control%luout_
    return
  end function luout

  subroutine set_luout (control, luout)
    integer(int_type), intent(in) :: luout
    type(control_type), intent(inout) :: control

    control%luout_ = luout
    return
  end subroutine set_luout

  ! logical unit for error output
  function luerr(control)
    integer(int_type) :: luerr
    type(control_type), intent(in) :: control

    luerr = control%luerr_
    return
  end function luerr

  subroutine set_luerr (control, luerr)
    integer(int_type), intent(in) :: luerr
    type(control_type), intent(inout) :: control

    control%luerr_ = luerr
    return
  end subroutine set_luerr

  ! logical unit for tty output
  function lutty(control)
    integer(int_type) :: lutty
    type(control_type), intent(in) :: control

    lutty = control%lutty_
    return
  end function lutty

  subroutine set_lutty (control, lutty)
    integer(int_type), intent(in) :: lutty
    type(control_type), intent(inout) :: control

    control%lutty_ = lutty
    return
  end subroutine set_lutty

  ! matrix-vector multiplications
  function matvecs(control)
    integer(int_type) :: matvecs
    type(control_type), intent(in) :: control

    matvecs = control%matvecs_
    return
  end function matvecs

  subroutine set_matvecs (control, matvecs)
    integer(int_type), intent(in) :: matvecs
    type(control_type), intent(inout) :: control

    control%matvecs_ = matvecs
    return
  end subroutine set_matvecs

  subroutine increment_matvecs (control)
    type(control_type), intent(inout) :: control

    control%matvecs_ = control%matvecs_ + 1
    return
  end subroutine increment_matvecs

  ! method
  function method(control)
    integer(int_type) :: method
    type(control_type), intent(in) :: control

    method = control%method_
    return
  end function method

  subroutine set_method (control, method)
    integer(int_type), intent(in) :: method
    type(control_type), intent(inout) :: control

    control%method_ = method
    return
  end subroutine set_method

  ! Krylov subspace size for, e.g. GMRES
  function nold(control)
    integer(int_type) :: nold
    type(control_type), intent(in) :: control

    nold = control%subspace_
    return
  end function nold

  function subspace(control)
    integer(int_type) :: subspace
    type(control_type), intent(in) :: control

    subspace = control%subspace_
    return
  end function subspace

  subroutine set_nold (control, nold)
    integer(int_type), intent(in) :: nold
    type(control_type), intent(inout) :: control

    control%subspace_ = nold
    return
  end subroutine set_nold

  subroutine set_subspace (control, subspace)
    integer(int_type), intent(in) :: subspace
    type(control_type), intent(inout) :: control

    control%subspace_ = subspace
    return
  end subroutine set_subspace

  ! which norm type to use
  function norm_type(control)
    integer(int_type) :: norm_type
    type(control_type), intent(in) :: control

    norm_type = control%norm_type_
    return
  end function norm_type

  subroutine set_norm_type (control, norm_type)
    integer(int_type), intent(in) :: norm_type
    type(control_type), intent(inout) :: control

    if (norm_type >= 0 .and. norm_type <= 2) then
       control%norm_type_ = norm_type
    else
       ! warning?, die?
    end if

    return
  end subroutine set_norm_type

  ! output mode
  function outmode(control)
    integer(int_type) :: outmode
    type(control_type), intent(in) :: control

    outmode = control%outmode_
    return
  end function outmode

  function outmode_none(control)
    logical(logical_type) :: outmode_none
    type(control_type), intent(in) :: control

    outmode_none = control%outmode_none_
    return
  end function outmode_none

  function outmode_errors(control)
    logical(logical_type) :: outmode_errors
    type(control_type), intent(in) :: control

    outmode_errors = control%outmode_errors_
    return
  end function outmode_errors

  subroutine set_outmode_none (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_none__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_none

  subroutine set_outmode_errors (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_errors__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_errors

  function outmode_warnings(control)
    logical(logical_type) :: outmode_warnings
    type(control_type), intent(in) :: control

    outmode_warnings = control%outmode_warnings_
    return
  end function outmode_warnings

  subroutine set_outmode_warnings (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_warnings__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_warnings

  function outmode_iterates(control)
    logical(logical_type) :: outmode_iterates
    type(control_type), intent(in) :: control

    outmode_iterates = control%outmode_iterates_
    return
  end function outmode_iterates

  subroutine set_outmode_iterates (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_iterates__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_iterates

  function outmode_summary(control)
    logical(logical_type) :: outmode_summary
    type(control_type), intent(in) :: control

    outmode_summary = control%outmode_summary_
    return
  end function outmode_summary

  subroutine set_outmode_summary (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_summary__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_summary

  function outmode_system(control)
    logical(logical_type) :: outmode_system
    type(control_type), intent(in) :: control

    outmode_system = control%outmode_system_
    return
  end function outmode_system

  subroutine set_outmode_system (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_system__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_system

  function outmode_full(control)
    logical(logical_type) :: outmode_full
    type(control_type), intent(in) :: control

    outmode_full = control%outmode_full_
    return
  end function outmode_full

  subroutine set_outmode_full (control)
    type(control_type), intent(inout) :: control

    control%outmode_ = outmode_full__
    call set_outmode_bools (control)
    return
  end subroutine set_outmode_full

  ! output frequency
  function output_frequency(control)
    integer(int_type) :: output_frequency
    type(control_type), intent(in) :: control

    output_frequency = control%output_frequency_
    return
  end function output_frequency

  subroutine set_output_frequency (control, frq)
    integer(int_type), intent(in) :: frq
    type(control_type), intent(inout) :: control

    control%output_frequency_ = frq
    return
  end subroutine set_output_frequency

  ! preconditionings
  function preconds(control)
    integer(int_type) :: preconds
    type(control_type), intent(in) :: control

    preconds = control%preconds_
    return
  end function preconds

  subroutine set_preconds (control, preconds)
    integer(int_type), intent(in) :: preconds
    type(control_type), intent(inout) :: control

    control%preconds_ = preconds
    return
  end subroutine set_preconds

  subroutine increment_preconds (control)
    type(control_type), intent(inout) :: control

    control%preconds_ = control%preconds_ + 1
    return
  end subroutine increment_preconds

  ! control residual updating
  function residual(control)
    integer(int_type) :: residual
    type(control_type), intent(in) :: control

    residual = control%residual_
    return
  end function residual

  subroutine set_residual (control, residual)
    integer(int_type), intent(in) :: residual
    type(control_type), intent(inout) :: control

    control%residual_ = residual
    call set_residual_bools (control)
    return
  end subroutine set_residual

  function residual_default(control)
    logical(logical_type) :: residual_default
    type(control_type), intent(in) :: control

    residual_default = control%residual_default_
    return
  end function residual_default

  subroutine set_residual_default (control)
    type(control_type), intent(inout) :: control

    control%residual_ = residual_default__
    call set_residual_bools (control)
    return
  end subroutine set_residual_default

  function residual_compute_only(control)
    logical(logical_type) :: residual_compute_only
    type(control_type), intent(in) :: control

    residual_compute_only = control%residual_compute_only_
    return
  end function residual_compute_only

  subroutine set_residual_compute_only (control)
    type(control_type), intent(inout) :: control

    control%residual_ = residual_compute_only__
    call set_residual_bools (control)
    return
  end subroutine set_residual_compute_only

  function residual_update(control)
    logical(logical_type) :: residual_update
    type(control_type), intent(in) :: control

    residual_update = control%residual_update_
    return
  end function residual_update

  subroutine set_residual_update (control)
    type(control_type), intent(inout) :: control

    control%residual_ = residual_update__
    call set_residual_bools (control)
    return
  end subroutine set_residual_update

  function residual_compute_update(control)
    logical(logical_type) :: residual_compute_update
    type(control_type), intent(in) :: control

    residual_compute_update = control%residual_compute_update_
    return
  end function residual_compute_update

  subroutine set_residual_compute_upd (control)
    type(control_type), intent(inout) :: control

    control%residual_ = residual_compute_update__
    call set_residual_bools (control)
    return
  end subroutine set_residual_compute_upd

  ! local/global scope
  function scope(control)
    integer(int_type) :: scope
    type(control_type), intent(in) :: control

    scope = control%scope_
    return
  end function scope

  function scope_local(control)
    use Ubik_ParallelInfo_module, only: scope_global__, scope_local__
    logical(logical_type) :: scope_local
    type(control_type), intent(in) :: control

    scope_local = ( control%scope_ == scope_local__ )
    return
  end function scope_local

  function scope_global(control)
    use Ubik_ParallelInfo_module, only: scope_global__, scope_local__
    logical(logical_type) :: scope_global
    type(control_type), intent(in) :: control

    scope_global = ( control%scope_ == scope_global__ )
    return
  end function scope_global

  subroutine set_scope (control, scope)
    integer(int_type), intent(in) :: scope
    type(control_type), intent(inout) :: control

    control%scope_ = scope
    return
  end subroutine set_scope

  subroutine set_scope_local (control)
    use Ubik_ParallelInfo_module, only: scope_global__, scope_local__

    type(control_type), intent(inout) :: control

    control%scope_ = scope_local__
    return
  end subroutine set_scope_local

  subroutine set_scope_global (control)
    use Ubik_ParallelInfo_module, only: scope_global__, scope_local__

    type(control_type), intent(inout) :: control

    control%scope_ = scope_global__
    return
  end subroutine set_scope_global

  ! number of steps to take in multistep preconditioners
  function steps(control)
    integer(int_type) :: steps
    type(control_type), intent(in) :: control

    steps = control%steps_
    return
  end function steps

  subroutine set_steps (control, steps)
    integer(int_type), intent(in) :: steps
    type(control_type), intent(inout) :: control

    control%steps_ = steps
    return
  end subroutine set_steps

  ! stopping test
  function stopping_true(control)
    logical(logical_type) :: stopping_true
    type(control_type), intent(in) :: control

    stopping_true = control%stopping_true_
    return
  end function stopping_true

  function stopping_relchg(control)
    logical(logical_type) :: stopping_relchg
    type(control_type), intent(in) :: control

    stopping_relchg = control%stopping_relchg_
    return
  end function stopping_relchg

  function stopping_axb(control)
    logical(logical_type) :: stopping_axb
    type(control_type), intent(in) :: control

    stopping_axb = control%stopping_axb_
    return
  end function stopping_axb

  function stopping_b(control)
    logical(logical_type) :: stopping_b
    type(control_type), intent(in) :: control

    stopping_b = control%stopping_b_
    return
  end function stopping_b

  function stopping_x(control)
    logical(logical_type) :: stopping_x
    type(control_type), intent(in) :: control

    stopping_x = control%stopping_x_
    return
  end function stopping_x

  function stopping_r0(control)
    logical(logical_type) :: stopping_r0
    type(control_type), intent(in) :: control

    stopping_r0 = control%stopping_r0_
    return
  end function stopping_r0

  function stopping_r(control)
    logical(logical_type) :: stopping_r
    type(control_type), intent(in) :: control

    stopping_r = control%stopping_r_
    return
  end function stopping_r

  subroutine set_stopping_relchg (control)
    type(control_type), intent(inout) :: control

    control%stopping_ = stopping_relchg__
    call set_stopping_bools (control)
    return
  end subroutine set_stopping_relchg

  subroutine set_stopping_axb (control)
    type(control_type), intent(inout) :: control

    control%stopping_ = stopping_axb__
    call set_stopping_bools (control)
    return
  end subroutine set_stopping_axb

  subroutine set_stopping_b (control)
    type(control_type), intent(inout) :: control

    control%stopping_ = stopping_b__
    call set_stopping_bools (control)
    return
  end subroutine set_stopping_b

  subroutine set_stopping_x (control)
    type(control_type), intent(inout) :: control

    control%stopping_ = stopping_x__
    call set_stopping_bools (control)
    return
  end subroutine set_stopping_x

  subroutine set_stopping_r0 (control)
    type(control_type), intent(inout) :: control

    control%stopping_ = stopping_r0__
    call set_stopping_bools (control)
    return
  end subroutine set_stopping_r0

  subroutine set_stopping_r (control)
    type(control_type), intent(inout) :: control

    control%stopping_ = stopping_r__
    call set_stopping_bools (control)
    return
  end subroutine set_stopping_r

  subroutine set_stopping_true (control)
    type(control_type), intent(inout) :: control

    control%stopping_true_ = .TRUE.
    return
  end subroutine set_stopping_true

  subroutine set_stopping_recursive (control)
    type(control_type), intent(inout) :: control

    control%stopping_true_ = .FALSE.
    return
  end subroutine set_stopping_recursive

  subroutine check_stopping_common (control)
    use Ubik_WriteString_module

    type(control_type), intent(inout) :: control
    integer(int_type) :: luerr_

    luerr_ = luerr(control)

    if (stopping_axb(control) .and. anorm(control) <= smallest(control)) then
       if (outmode_warnings(control)) then
          call spew (luerr_, ' ')
          call spew (luerr_, '*************************************')
          call spew (luerr_, 'WARNING')
          call spew (luerr_, '  Attempt to use stopping test based on')
          call spew (luerr_, '  || r || / ( || A ||*|| x || + || b || ) with an')
          call spew (luerr_, '  invalid value for norm of coefficient, || A ||.')
          call spew (luerr_, ' ')
          call spew (luerr_, '  Make sure control%anorm is set to || A || prior to')
          call spew (luerr_, '  calling this routine.')
          call spew (luerr_, ' ')
          call spew (luerr_, '  Changing stopping test to || r || / || b ||...')
          call spew (luerr_, '*************************************')
       end if
       call set_stopping_b (control)
    end if

    return
  end subroutine check_stopping_common

  subroutine check_stopping_nonstationary (control)
    use Ubik_WriteString_module

    type(control_type), intent(inout) :: control
    integer(int_type) :: luerr_

    luerr_ = luerr(control)

    call check_stopping_common (control)

    ! relative change stopping test is inappropriate for nonstationary methods,
    ! so warn against it, but honor the user's wishes
    if (stopping_relchg(control)) then
       if (outmode_warnings(control)) then
          call spew (luerr_, ' ')
          call spew (luerr_, '*******************************')
          call spew (luerr_, 'WARNING')
          call spew (luerr_, '  Using a stopping test based on relative change')
          call spew (luerr_, '  in iterates, || x - xold || / || x ||.')
          call spew (luerr_, ' ')
          call spew (luerr_, '  This stopping test is not recommended with this solver,')
          call spew (luerr_, '  but if that is what you really want, ok.')
          call spew (luerr_, '*******************************')
       end if
    end if

    return
  end subroutine check_stopping_nonstationary

  subroutine check_stopping_stationary (control)
    use Ubik_WriteString_module

    type(control_type), intent(inout) :: control
    integer(int_type) :: luerr_

    luerr_ = luerr(control)

    call check_stopping_common (control)

    ! Check to make sure that if one of the residual-based stopping
    ! tests is chosen it is one based on the true residual.  If not,
    ! correct and continue.
    if (.not.stopping_relchg(control)) then
       if (outmode_warnings(control)) then
          call spew (luerr_, ' ')
          call spew (luerr_, '*******************************')
          call spew (luerr_, 'WARNING')
          call spew (luerr_, '  Using stopping test based on residual.')
          call spew (luerr_, ' ')
          call spew (luerr_, '  Since a residual estimate is not built as a natural')
          call spew (luerr_, '  by-product of stationary iterations we must calculate')
          call spew (luerr_, '  b-Ax explicitly each iteration.')
          call spew (luerr_, ' ')
          call spew (luerr_, '  Note that this is costly...')
          call spew (luerr_, '*******************************')
       end if
       if (.not.stopping_true(control)) then
          call set_stopping_true (control)
       end if
    end if

    return
  end subroutine check_stopping_stationary

  ! norm of coefficient
  function anorm(control)
    real(real_type) :: anorm
    type(control_type), intent(in) :: control

    anorm = control%anorm_
    return
  end function anorm

  subroutine set_anorm (control, anorm)
    real(real_type) :: anorm
    type(control_type), intent(inout) :: control

    control%anorm_ = anorm
    return
  end subroutine set_anorm

  ! norm of source vector
  function bnorm(control)
    real(real_type) :: bnorm
    type(control_type), intent(in) :: control

    bnorm = control%bnorm_
    return
  end function bnorm

  subroutine set_bnorm (control, bnorm)
    real(real_type) :: bnorm
    type(control_type), intent(inout) :: control

    control%bnorm_ = bnorm
    return
  end subroutine set_bnorm

  ! convergence criterion
  function eps(control)
    real(real_type) :: eps
    type(control_type), intent(in) :: control

    eps = control%eps_
    return
  end function eps

  subroutine set_eps (control, eps)
    real(real_type) :: eps
    type(control_type), intent(inout) :: control

    control%eps_ = eps
    return
  end subroutine set_eps

  ! relaxation parameter
  function omega(control)
    real(real_type) :: omega
    type(control_type), intent(in) :: control

    omega = control%omega_
    return
  end function omega

  subroutine set_omega (control, omega)
    real(real_type) :: omega
    type(control_type), intent(inout) :: control

    control%omega_ = omega
    return
  end subroutine set_omega

  ! small numbers
  function small(control)
    real(real_type) :: small
    type(control_type), intent(in) :: control

    small = control%small_
    return
  end function small

  subroutine set_small (control, small_number)
    real(real_type), intent(in) :: small_number
    type(control_type), intent(inout) :: control

    control%small_ = small_number
    return
  end subroutine set_small

  function smallest(control)
    real(real_type) :: smallest
    type(control_type), intent(in) :: control

    smallest = control%smallest_
    return
  end function smallest

  subroutine set_smallest (control, smallest_number)
    real(real_type), intent(in) :: smallest_number
    type(control_type), intent(inout) :: control

    control%smallest_ = smallest_number
    return
  end subroutine set_smallest

  ! xold
  function xold(control)
    type(control_type), intent(in) :: control
    real(real_type), dimension(SIZE(control%xold_)) :: xold

    xold = control%xold_
    return
  end function xold

  subroutine set_xold (control, x)
    real(real_type), dimension(:) :: x
    type(control_type), intent(inout) :: control
    integer(int_type) :: status

    if (.not.ASSOCIATED(control%xold_)) then
       allocate (control%xold_(SIZE(x)), stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(control)
          return
       end if
    end if

    control%xold_ = x
    return
  end subroutine set_xold

  ! status
  function status(control)
    real(real_type) :: status
    type(control_type), intent(in) :: control

    status = control%status_
    return
  end function status

  subroutine set_status (control, status)
    integer(int_type), intent(in) :: status
    type(control_type), intent(inout) :: control

    control%status_ = status
    return
  end subroutine set_status

  function converged(control)
    use ubik_constants_module

    logical(logical_type) :: converged
    type(control_type), intent(in) :: control

    converged = (control%status_ == converged__)
    return
  end function converged

  function invalid_input(control)
    use ubik_constants_module

    logical(logical_type) :: invalid_input
    type(control_type), intent(in) :: control

    invalid_input = (control%status_ == invalid_input__)
    return
  end function 

  function alloc_failure_control(control)
    use ubik_constants_module

    logical(logical_type) :: alloc_failure_control
    type(control_type), intent(in) :: control

    alloc_failure_control = (control%status_ == alloc_failure__)
    return
  end function alloc_failure_control

  function internal_error(control)
    use ubik_constants_module

    logical(logical_type) :: internal_error
    type(control_type), intent(in) :: control

    internal_error = (control%status_ == internal_error__)
    return
  end function internal_error

  function itmax_exceeded(control)
    use ubik_constants_module

    logical(logical_type) :: itmax_exceeded
    type(control_type), intent(in) :: control

    itmax_exceeded = (control%status_ == itmax_exceeded__)
    return
  end function itmax_exceeded

  function breakdown(control)
    use ubik_constants_module

    logical(logical_type) :: breakdown
    type(control_type), intent(in) :: control

    breakdown = (control%status_ == breakdown__)
    return
  end function breakdown

  function matmul_error(control)
    use ubik_constants_module

    logical(logical_type) :: matmul_error
    type(control_type), intent(in) :: control

    matmul_error = (control%status_ == matmul_error__)
    return
  end function matmul_error

  function precond_error(control)
    use ubik_constants_module

    logical(logical_type) :: precond_error
    type(control_type), intent(in) :: control

    precond_error = (control%status_ == precond_error__)
    return
  end function precond_error

  function next_iter(control)
    use ubik_constants_module

    logical(logical_type) :: next_iter
    type(control_type), intent(in) :: control

    next_iter = (control%status_ == next_iter__)
    return
  end function next_iter

  subroutine set_status_success(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = success__
    return
  end subroutine set_status_success

  subroutine set_status_converged(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = converged__
    return
  end subroutine set_status_converged

  subroutine set_status_invalid_input(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = invalid_input__
    return
  end subroutine set_status_invalid_input

  subroutine set_status_alloc_failure_cntrl(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = alloc_failure__
    return
  end subroutine set_status_alloc_failure_cntrl

  subroutine set_status_internal_error(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = internal_error__
    return
  end subroutine set_status_internal_error

  subroutine set_status_itmax_exceeded(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = itmax_exceeded__
    return
  end subroutine set_status_itmax_exceeded

  subroutine set_status_breakdown(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = breakdown__
    return
  end subroutine set_status_breakdown

  subroutine set_status_matmul_error(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = matmul_error__
    return
  end subroutine set_status_matmul_error

  subroutine set_status_precond_error(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = precond_error__
    return
  end subroutine set_status_precond_error

  subroutine set_status_next_iter(control)
    use ubik_constants_module

    type(control_type), intent(inout) :: control

    control%status_ = next_iter__
    return
  end subroutine set_status_next_iter

  ! cpu time
  function cpu_default(control)
    real(real_type) :: cpu_default
    type(control_type), intent(in) :: control

    cpu_default = control%cpu_(control%iter_)
    return
  end function cpu_default

  function cpu_iter(control, it)
    real(real_type) :: cpu_iter
    integer(int_type), intent(in) :: it
    type(control_type), intent(in) :: control

    cpu_iter = control%cpu_(it)
    return
  end function cpu_iter

  subroutine set_cpu_default (control)
    type(control_type), intent(inout) :: control

    call CPU_Time (control%cpu_(control%iter_))
    return
  end subroutine set_cpu_default

  subroutine set_cpu_iter (control, it)
    integer(int_type), intent(in) :: it
    type(control_type), intent(inout) :: control

    call CPU_Time (control%cpu_(it))
    return
  end subroutine set_cpu_iter

  subroutine normalize_cpu (control)
    integer(int_type) :: it
    type(control_type), intent(inout) :: control

    it = control%iter_
    control%cpu_(1:it) = control%cpu_(1:it) - control%cpu_(0)

    return
  end subroutine normalize_cpu

  ! residual norm
  function rnorm_default(control)
    real(real_type) :: rnorm_default
    type(control_type), intent(in) :: control

    rnorm_default = control%rnorm_(control%iter_)
    return
  end function rnorm_default

  function rnorm_iter(control, iteration)
    real(real_type) :: rnorm_iter
    integer(int_type), intent(in) :: iteration
    type(control_type), intent(in) :: control

    rnorm_iter = control%rnorm_(iteration)
    return
  end function rnorm_iter

  subroutine set_rnorm_default (control, rnorm)
    real(real_type), intent(in) :: rnorm
    type(control_type), intent(inout) :: control

    call set_rnorm (control, control%iter_, rnorm)
    return
  end subroutine set_rnorm_default

  subroutine set_rnorm_iter (control, it, rnorm)
    integer(int_type), intent(in) :: it
    real(real_type), intent(in) :: rnorm
    type(control_type), intent(inout) :: control

    control%rnorm_(it) = rnorm
    return
  end subroutine set_rnorm_iter

  ! true residual norm
  function rnormt_default(control)
    real(real_type) :: rnormt_default
    type(control_type), intent(in) :: control

    rnormt_default = control%rnormt_(control%iter_)
    return
  end function rnormt_default

  function rnormt_iter(control, it)
    real(real_type) :: rnormt_iter
    integer(int_type), intent(in) :: it
    type(control_type), intent(in) :: control

    rnormt_iter = control%rnormt_(it)
    return
  end function rnormt_iter

  subroutine set_rnormt_default (control, rnormt)
    real(real_type), intent(in) :: rnormt
    type(control_type), intent(inout) :: control

    call set_rnormt (control, control%iter_, rnormt)
    return
  end subroutine set_rnormt_default

  subroutine set_rnormt_iter (control, it, rnormt)
    integer(int_type), intent(in) :: it
    real(real_type), intent(in) :: rnormt
    type(control_type), intent(inout) :: control

    control%rnormt_(it) = rnormt
    return
  end subroutine set_rnormt_iter

  ! error
  function err_default(control)
    real(real_type) :: err_default
    type(control_type), intent(in) :: control

    err_default = control%err_(control%iter_)
    return
  end function err_default

  function err_iter(control, it)
    real(real_type) :: err_iter
    integer(int_type), intent(in) :: it
    type(control_type), intent(in) :: control

    err_iter = control%err_(it)
    return
  end function err_iter

  subroutine set_err_default (control, err)
    real(real_type), intent(in) :: err
    type(control_type), intent(inout) :: control

    control%err_(control%iter_) = err
    return
  end subroutine set_err_default

  subroutine set_err_iter (control, it, err)
    integer(int_type), intent(in) :: it
    real(real_type), intent(in) :: err
    type(control_type), intent(inout) :: control

    control%err_(it) = err
    return
  end subroutine set_err_iter

  ! true error
  function errt_default(control)
    real(real_type) :: errt_default
    type(control_type), intent(in) :: control

    errt_default = control%errt_(control%iter_)
    return
  end function errt_default

  function errt_iter(control, it)
    real(real_type) :: errt_iter
    integer(int_type), intent(in) :: it
    type(control_type), intent(in) :: control

    errt_iter = control%errt_(it)
    return
  end function errt_iter

  subroutine set_errt_default (control, errt)
    real(real_type), intent(in) :: errt
    type(control_type), intent(inout) :: control

    control%errt_(control%iter_) = errt
    return
  end subroutine set_errt_default

  subroutine set_errt_iter (control, it, errt)
    integer(int_type), intent(in) :: it
    real(real_type), intent(in) :: errt
    type(control_type), intent(inout) :: control

    control%errt_(it) = errt
    return
  end subroutine set_errt_iter

end module ubik_control_type_module
