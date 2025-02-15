!!
!! CS_MINRES_SOLVER_TYPE
!!
!! This module defines a linear solver for complex symmetric linear systems
!! that uses the CS-MINRES/CS-MINRES-QLP Krylov subspace method [1].
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! [1] Sou-Cheng Choi, "Minimal Residual Methods for Complex Symmetric,
!!     Skew Symmetric, and Skew Hermitian Systems", arXiv:1304.6782, (2014)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NB: the preliminary solve procedure is a direct translation of the matlab
!! code and needs very badly to be reorganized and cleaned up.
!!

module cs_minres_solver2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op2_class
  use parameter_list_type
  use mpi
  implicit none
  private

  type, public :: cs_minres_solver2
    integer :: maxit
    logical :: show
    real(r8) :: rtol, maxxnorm, Acondlim, TranCond
    ! solver diagnostics
    integer :: flag, iter, Miter, QLPiter
    real(r8) :: relres, relAres, Anorm, Acond, xnorm, Axnorm
    real(r8), allocatable :: resvec(:), Aresvec(:)
    integer :: comm = MPI_COMM_WORLD, nproc, rank
  contains
    procedure :: init
    procedure :: solve
  end type

  interface num2str
    procedure :: dnum2str, znum2str
  end interface

contains

  subroutine init(this, params, comm)

    class(cs_minres_solver2), intent(out) :: this
    type(parameter_list) :: params
    integer, intent(in), optional :: comm

    integer :: ierr

    if (present(comm)) this%comm = comm
    call MPI_Comm_rank(this%comm, this%rank, ierr)
    call MPI_Comm_size(this%comm, this%nproc, ierr)

    call params%get('rel-tol', this%rtol, default=epsilon(1.0_r8))
    call params%get('max-iter', this%maxit, default=10000)
    call params%get('maxxnorm', this%maxxnorm, default=1e7_r8)
    call params%get('Acondlim', this%Acondlim, default=1e15_r8)
    call params%get('TranCond', this%TranCond, default=1e7_r8)
    call params%get('show', this%show, default=.true.)
    this%show = this%show .and. (this%rank == 0)

  end subroutine

!!!function [x, out_param] = csminresqlp(varargin)
!!!%csminresqlp: min-length solution to complex symmetric (possibly singular) Ax=b or min||Ax-b||.
!!!%
!!!%   X = csminresqlp(A,B) solves the system of linear equations A*X=B
!!!%   or the least-squares problem min norm(B-A*X) if A is singular.
!!!%   The N-by-N matrix A must be real symmetric or complex symmetric, but
!!!%   need not be positive definite or nonsingular.  It may be double or single.
!!!%   The rhs vector B must have length N.  It may be real or complex,
!!!%   double or single.
!!!%
!!!%   X = csminresqlp(AFUN,B) accepts a function handle AFUN instead of
!!!%   the matrix A.  Y = AFUN(X) returns the matrix-vector product Y=A*X.
!!!%   In all of the following syntaxes, A can be replaced by AFUN.
!!!%
!!!%   X = csminresqlp(A,B,RTOL) specifies a stopping tolerance.
!!!%   If RTOL=[] or is absent, a default value is used.
!!!%   (Similarly for all later input parameters.)
!!!%   Default RTOL=1e-15.
!!!%
!!!%   X = csminresqlp(A,B,RTOL,MAXIT)
!!!%   specifies the maximum number of iterations.  Default MAXIT=4*N.
!!!%
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M) and
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M1,M2)uses a matrix M or M=M1*M2 as
!!!%   preconditioner.  M must be positive definite and symmetric or Hermitian.
!!!%   It may be a function handle MFUN such that Y=MFUN(X) returns Y=M\X.
!!!%   If M=[], a preconditioner is not applied.
!!!%
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M,SHIFT) and
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M1,M2,SHIFT)
!!!%   solves (A - SHIFT*I)X = B, or the corresponding least-squares problem
!!!%   if (A - SHIFT*I) is singular, where SHIFT is a real or complex scalar.
!!!%   Default SHIFT=0.
!!!%
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M,SHIFT,MAXXNORM,ACONDLIM,TRANCOND) and
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M1,M2,SHIFT,MAXXNORM,ACONDLIM,TRANCOND)
!!!%   specifies three parameters associated with singular or
!!!%   ill-conditioned systems (A - SHIFT*I)*X = B.
!!!%
!!!%   MAXXNORM is an upper bound on NORM(X).
!!!%   Default MAXXNORM=1e7.
!!!%
!!!%   ACONDLIM is an upper bound on ACOND, an estimate of COND(A).
!!!%   Default ACONDLIM=1e15.
!!!%
!!!%   TRANCOND is a real scalar >= 1.
!!!%   If TRANCOND>1,        a switch is made from MINRES iterations to
!!!%                         CS-MINRES iterationsd when ACOND >= TRANCOND.
!!!%   If TRANCOND=1,        all iterations will be CS-MINRES iterations.
!!!%   If TRANCOND=ACONDLIM, all iterations will be conventional MINRES
!!!%                         iterations (which are slightly cheaper).
!!!%   Default TRANCOND=1e7.
!!!%
!!!%   X = csminresqlp(A,B,RTOL,MAXIT,M,SHIFT,MAXXNORM,ACONDLIM,TRANCOND,SHOW)
!!!%   specifies the printing option.
!!!%   If SHOW=true,  an iteration log will be output.
!!!%   If SHOW=false, the log is suppressed.
!!!%   Default SHOW=true.
!!!%
!!!%
!!!%   [X,FLAG] = csminresqlp(A,B,...) returns a convergence FLAG:
!!!%   -1 (beta2=0)  B and X are eigenvectors of (A - SHIFT*I).
!!!%    0 (beta1=0)  B = 0.  The exact solution is X = 0.
!!!%    1 X solves the compatible (possibly singular) system (A - SHIFT*I)X = B
!!!%      to the desired tolerance:
!!!%         RELRES = RNORM / (ANORM*XNORM + NORM(B)) <= RTOL,
!!!%      where
!!!%              R = B - (A - SHIFT*I)X and RNORM = norm(R).
!!!%    2 X solves the incompatible (singular) system (A - SHIFT*I)X = B
!!!%      to the desired tolerance:
!!!%         RELARES = ARNORM / (ANORM * RNORM) <= RTOL,
!!!%      where
!!!%              AR = (A - SHIFT*I)R and ARNORM = NORM(AR).
!!!%    3 Same as 1 with RTOL = EPS.
!!!%    4 Same as 2 with RTOL = EPS.
!!!%    5 X converged to an eigenvector of (A - SHIFT*I).
!!!%    6 XNORM exceeded MAXXNORM.
!!!%    7 ACOND exceeded ACONDLIM.
!!!%    8 MAXIT iterations were performed before one of the previous
!!!%      conditions was satisfied.
!!!%    9 The system appears to be exactly singular.  XNORM does not
!!!%      yet exceed MAXXNORM, but would if further iterations were
!!!%      performed.
!!!%
!!!%    [X,FLAG,ITER,MITER,QLPITER] = csminresqlp(A,B,...) returns the
!!!%    number of iterations performed, with ITER = MITER + QLPITER.
!!!%    MITER   is the number of conventional MINRES iterations.
!!!%    QLPITER is the number of CS-MINRES-QLP iterations.
!!!%
!!!%    [X,FLAG,ITER,MITER,QLPITER,RELRES,RELARES] = csminresqlp(A,B,...)
!!!%    returns relative residuals for (A - SHIFT*I)X = B and the
!!!%    associated least-squares problem.  RELRES and RELARES are
!!!%    defined above in the description of FLAG.
!!!%
!!!%    [X,FLAG,ITER,MITER,QLPITER,RELRES,RELARES,ANORM,ACOND,XNORM,AXNORM] =
!!!%       csminresqlp(A,B,...) returns
!!!%       ANORM,  an estimate of the 2-norm of A-SHIFT*I.
!!!%       ACOND,  an estimate of COND(A-SHIFT*I,2).
!!!%       XNORM,  a recurred estimate of NORM(X).
!!!%       AXNORM, a recurred estimate of NORM((A-SHIFT*I)X)
!!!%
!!!%    [X,FLAG,ITER,MITER,QLPITER,RELRES,RELARES,ANORM,ACOND,XNORM,AXNORM,...
!!!%       RESVEC,ARESVEC] = csminresqlp(A,B,...) returns
!!!%       RESVEC,  a vector of estimates of NORM(R) at each iteration,
!!!%                including NORM(B) as the first entry.
!!!%       ARESVEC, a vector of estimates of NORM((A-SHIFT*I)R) at each
!!!%                iteration, including NORM((A-SHIFT*I)B) as the first entry.
!!!%       RESVEC and ARESVEC have length ITER+1.
!!!%
!!!%
!!!% EXAMPLE 1:
!!!%   n = 100;                                       e = ones(n,1);
!!!%   A = i * spdiags([-2*e 4*e -2*e],-1:1,n,n);    M = spdiags(4*e,0,n,n);
!!!%
!!!%   b = i * sum(A,2);             rtol = 1e-10;   maxit = 50;
!!!%   x = csminresqlp(A,b,rtol,maxit,M);
!!!%
!!!%   Alternatively, use this matrix-vector product function:
!!!%     function y = Afun(x,n)
!!!%       y = 4*x;
!!!%       y(2:n)   = y(2:n)   - 2*x(1:n-1);
!!!%       y(1:n-1) = y(1:n-1) - 2*x(2:n);
!!!%       y = i * y;
!!!%   as input to csminresqlp:
!!!%   n = 100;
!!!%   A = @(x)Afun(x,n);
!!!%   x = csminresqlp(A,b,rtol,maxit,M);
!!!%
!!!% EXAMPLE 2: i*A, where A is Laplacian on a 50 by 50 grid, singular and indefinite.
!!!%   n = 50;  N = n^2;  e = ones(n,1);
!!!%   B = spdiags([e e e], -1:1, n, n);
!!!%   A = i * sparse([],[],[],N,N,(3*n-2)^2);
!!!%   for i=1:n
!!!%     A((i-1)*n+1:i*n,(i-1)*n+1:i*n) = B;
!!!%     if i*n+1 < n*n,   A(i*n+1:(i+1)*n,(i-1)*n+1:i*n)     = B; end
!!!%     if (i-2)*n+1 > 0, A((i-2)*n+1:(i-1)*n,(i-1)*n+1:i*n) = B; end
!!!%   end
!!!%   b = i * sum(A,2);   rtol   = 1e-5;    shift = 0;    maxxnorm = 1e5;
!!!%   M = [];         Acondlim = [];    show  = true;     tranCond = 1e7;
!!!%   x = csminresqlp(A,b,rtol,N,M,shift,maxxnorm,Acondlim,tranCond,show);
!!!%
!!!% EXAMPLE 3: A is diagonal, singular and indefinite.
!!!%   h = 1;  a = -10; b = -a; n = 2*b/h + 1;
!!!%   A = i * spdiags((a:h:b)', 0, n, n);
!!!%   b = i * ones(n,1);  rtol   = 1e-6;    shift = 0;    maxxnorm = 1e2;
!!!%   M = [];         Acondlim = [];    show  = true;     tranCond = 1e7;
!!!%   x = csminresqlp(A,b,rtol,N,M,shift,maxxnorm,Acondlim,tranCond,show);
!!!%
!!!% See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, PCG, QMR, SYMMLQ,
!!!% TFQMR, CHOLINC, FUNCTION_HANDLE.
!!!% Also MINRES-QLP, MINRES, SYMMLQ, LSQR, CGLS downloadable from
!!!%   http://www.stanford.edu/group/SOL/software.html
!!!%
!!!% REFERENCES:
!!!% Sou-Cheng T. Choi,
!!!% Minimal Residual Methods for Complex Symmetric, Skew Symmetric,
!!!% and Skew Hermitian Systems. Report ANL/MCS-P3028-0812,
!!!% Computation Institute, University of Chicago, 2013.
!!!% http://home.uchicago.edu/sctchoi/CSMINRES20.pdf
!!!%
!!!% PLEASE CITE:
!!!% The above reference and this Matlab software:
!!!%
!!!% Sou-Cheng T. Choi, "CS-MINRES-QLP, version 1.1" [Matlab/GNU-Octave Software], 2017.
!!!% Available from http://home.uchicago.edu/sctchoi/
!!!%
!!!% CURRENT / FUTURE RELEASES of csminresqlp:
!!!%   http://home.uchicago.edu/sctchoi/
!!!%   http://mypages.iit.edu/~schoi32/
!!!%   http://www.mathworks.com/matlabcentral/fileexchange/61151-cs-minres-qlp
!!!
!!!% MODIFICATION HISTORY:
!!!%   11 Apr 2011: Created CSMINRES.m from MINRESQLPs.m
!!!%
!!!% KNOWN BUGS:
!!!%   DD Mmm YYYY: ---
!!!%
!!!% NOTES:
!!!%   DD Mmm YYYY: ---
!!!%
!!!% AUTHOR: Sou-Cheng (Terrya) Choi, Computation Institute, University of
!!!% Chicago
!!!%
!!!% COPYRIGHT NOTICE:
!!!%
!!!%   This is Copyrighted Material. The software is COPYRIGHTED by the
!!!%   original authors.
!!!%
!!!% COPYLEFT NOTICE:
!!!%
!!!%   Permission is granted to make and distribute verbatim copies of this
!!!%   file, provided that the entire file is copied **together** as a
!!!%   **unit**.
!!!%
!!!%   The purpose of this permission notice is to allow you to make copies
!!!%   of the software and distribute them to others, for free or for a fee,
!!!%   subject to the constraint that you maintain the entire file here
!!!%   as a unit.  This enables people you give the software to be aware
!!!%   of its origins, to ask questions of us by e-mail, to request
!!!%   improvements, obtain later releases, and so forth.
!!!%
!!!%   If you seek permission to copy and distribute translations of this
!!!%   software into another language, please e-mail a specific request to
!!!%   saunders@stanford.edu and scchoi@stanford.edu.
!!!%
!!!%   If you seek permission to excerpt a **part** of the software library,
!!!%   for example to appear in a scientific publication, please e-mail a
!!!%   specific request to saunders@stanford.edu and scchoi@stanford.edu.
!!!%
!!!% COMMENTS?
!!!%
!!!%   Email sctchoi@uchicago.edu
!!!%
!!!% DISCLAIMER: This software is prorvided to the public free of charge. The
!!!% author does not provide any guarantee.
!!!%


  subroutine solve(this, lin_op, b, x)

    use zvector_class

    class(cs_minres_solver2), intent(inout) :: this
    class(complex_lin_op2), intent(inout) :: lin_op
    class(zvector), intent(in) :: b
    class(zvector), intent(inout) :: x

    !integer :: n, Miter, j
    integer :: Miter
    complex(r8) :: shift
    logical :: debug
    character(:), allocatable :: first, last, msg(:), head
    real(r8), allocatable :: resvec(:), Aresvec(:)

    !complex(r8), dimension(size(x)) :: r2, r3
    class(zvector), allocatable :: r2, r3
    real(r8) :: beta1

    integer :: flag0, flag, iter, QLPiter, lines, headlines
    real(r8) :: beta, betal, betan, gmin, cs, cr1, cr2
    complex(r8) :: sn, sr1, sr2, dbar, gbar
    real(r8) :: rnorm, xnorm, xnorml, xl2norm, xnorm_tmp, Axnorm, Anorm, Acond, relres, relresl, pnorm
    real(r8) :: abs_gama, Anorml, gminl, gminl2, Acondl, rnorml, rootl, Arnorm, Arnorml, relAres, relAresl, epsx, t1, t2
    real(r8) :: direct_rnorm, direct_xnorm

    complex(r8) :: phi, tau, taul, taul2
    complex(r8) :: gama, gamal, gamal2, gamal3, gama_tmp, gamal_tmp, gama_QLP, gamal_QLP
    complex(r8) :: epln, eplnn
    complex(r8) :: vepln, veplnl, veplnl2, vepln_QLP
    complex(r8) :: dlta, dltan, dlta_tmp, dlta_QLP
    complex(r8) :: eta, etal, etal2
    complex(r8) :: u, ul, ul2, ul3, ul4, u_QLP, ul_QLP

    !complex(r8), dimension(size(x)) ::  w, wl, wl2, r1, xl2
    class(zvector), allocatable :: w, wl, wl2, r1, xl2

    !complex(r8) :: v(size(x)), alfa
    complex(r8) ::alfa
    class(zvector), allocatable :: v, vbar

    complex(r8), parameter :: zzero = 0

    call x%clone(r2)
    call x%clone(r3)
    call x%clone(w)
    call x%clone(wl)
    call x%clone(wl2)
    call x%clone(r1)
    call x%clone(xl2)
    call x%clone(v)
    call x%clone(vbar)

    shift = 0.0_r8

    debug = .false.

!!!%%  Check inputs and set default values.
!!!
!!!precon = true;
!!!if (isempty(M) && (isempty(M1) || isempty(M2))),  precon = false; end
!!!% if nargin <  2,                             error('Not enough input parameters');  end
!!!% if nargin <  3 || ~exist('rtol'    ,'var') || isempty(rtol)    , rtol     = eps; end
!!!% if nargin <  4 || ~exist('maxit'   ,'var') || isempty(maxit)   , maxit    = 4*n;   end
!!!% if nargin <  5 || ~exist('M'       ,'var') || isempty(M)       , precon   = false; end
!!!% if nargin <  6 || ~exist('shift'   ,'var') || isempty(shift)   , shift    = 0;     end
!!!% if nargin <  7 || ~exist('maxxnorm','var') || isempty(maxxnorm), maxxnorm = 1e7;  end
!!!% if nargin <  8 || ~exist('Acondlim','var') || isempty(Acondlim), Acondlim = 1e15;  end
!!!% if nargin <  9 || ~exist('TranCond','var') || isempty(TranCond), TranCond = 1e7;   end
!!!% if nargin < 10 || ~exist('show'    ,'var') || isempty(show)    , show     = true;  end
!!!%if nargin< 11 || ~exist('disable' ,'var') || isempty(disable) , disable  = false; end

    allocate(resvec(this%maxit+1), source=0.0_r8)
    allocate(Aresvec(this%maxit+1), source=0.0_r8)

!!!%isComplexSymmetric = (~isreal(A)) && (normest(A-A.') < eps * normest(A));
!!!%isComplexSymmetric = isComplexSymmetric(A);
!!!if (ischar(A))
!!!  if (length(A) > 0)
!!!    A = str2func(A);
!!!  else
!!!    error('Empty string for function name A');
!!!  end
!!!end
!!!if (precon && ischar(M))
!!!  if (length(M) > 0)
!!!    M = str2func(M);
!!!  else
!!!    error('Empty string for function name M');
!!!  end
!!!end

!!!%% Set up {beta1, p, v} for the first Lanczos vector v1.
!!!r2    = full(b);    % r2    = b
!!!r3    = r2;         % r3    = b
!!!beta1 = norm(r2);   % beta1 = norm(b)
!!!if precon
!!!  persistent R
!!!  clear R
!!!  r3    = minresxxxM(M, M1, M2, r2,inputs{:});   % M*r3  = b
!!!  beta1 = r3'*r2;             % beta1 = b'*inv(M)*b
!!!  if beta1 < 0
!!!    error('"M" appears to be indefinite.');
!!!  else
!!!    beta1 = sqrt(beta1);
!!!  end
!!!end

  !r2(1:n) = b(1:n)
  !r3(1:n) = r2(1:n)
  call r2%copy(b)
  call r3%copy(r2)
  !beta1 = this%global_nrm2(r2)
  beta1 = r2%norm2()
  call lin_op%precon(r2, r3)
  !beta1 = this%global_dotc(r3, r2)  ! dotc value should be real
  beta1 = r3%dotc(r2) ! dotc value should be real
  if (beta1 < 0.0_r8) then
    this%flag = -3  ! previously unused value
    return
  end if
  beta1 = sqrt(beta1)

!!!%% Initialize other quantities.
!!!flag0    = -2;     flag     = flag0;
!!!iter     = 0;      QLPiter  = 0;
!!!lines    = 1;      headlines= 20;
!!!beta     = 0;      tau      = 0;          taul     = 0;      phi      = beta1;
!!!betan    = beta1;  gmin     = 0;          cs       = -1;     sn       = 0;
!!!cr1      = 1;      sr1      = 0;          cr2      = -1;     sr2      = 0;
!!!dltan    = 0;      eplnn    = 0;          gama     = 0;      gamal    = 0;
!!!gamal2   = 0;      eta      = 0;          etal     = 0;      etal2    = 0;
!!!vepln    = 0;      veplnl   = 0;          veplnl2  = 0;      ul3      = 0;
!!!ul2      = 0;      ul       = 0;          u        = 0;      rnorm    = betan;
!!!xnorm    = 0;      xl2norm  = 0;          Axnorm   = 0;
!!!Anorm    = 0;      Acond    = 1;
!!!relres   = rnorm / (beta1 + 1e-50);       % Safeguard for beta1 = 0
!!!x        = zeros(n,1);
!!!w        = zeros(n,1);
!!!wl       = zeros(n,1);
!!!r1       = zeros(n,1);

    flag0    = -2;     flag     = flag0
    iter     = 0;      QLPiter  = 0
    lines    = 1;      headlines= 20
    beta     = 0;      tau      = 0;          taul     = 0;      phi      = beta1
    betan    = beta1;  gmin     = 0;          cs       = -1;     sn       = 0
    cr1      = 1;      sr1      = 0;          cr2      = -1;     sr2      = 0
    dltan    = 0;      eplnn    = 0;          gama     = 0;      gamal    = 0
    gamal2   = 0;      eta      = 0;          etal     = 0;      etal2    = 0
    vepln    = 0;      veplnl   = 0;          veplnl2  = 0;      ul3      = 0
    ul2      = 0;      ul       = 0;          u        = 0;      rnorm    = betan
    xnorm    = 0;      xl2norm  = 0;          Axnorm   = 0
    Anorm    = 0;      Acond    = 1
    relres   = rnorm / (beta1 + 1e-50_r8) ! Safeguard for beta1 = 0

    !x(1:n) = 0.0_r8
    !w = 0.0_r8
    !wl = 0.0_r8
    !r1 = 0.0_r8
    call x%setzero
    call w%setzero
    call wl%setzero
    call r1%setzero

!!!if ~isempty(resvec)
!!!  resvec(1)  = beta1;
!!!end

    resvec(1) = beta1

!!!%% print header if show
!!!first = 'Enter csminresqlp.  ';
!!!last  = 'Exit csminresqlp.  ';
!!!msg=[' beta2 = 0.  b and x are eigenvectors                   '   % -1
!!!     ' beta1 = 0.  The exact solution is  x = 0               '   %  0
!!!     ' A solution to Ax = b found, given rtol                 '   %  1
!!!     ' Min-length solution for singular LS problem, given rtol'   %  2
!!!     ' A solution to Ax = b found, given eps                  '   %  3
!!!     ' Min-length solution for singular LS problem, given eps '   %  4
!!!     ' x has converged to an eigenvector                      '   %  5
!!!     ' xnorm has exceeded maxxnorm                            '   %  6
!!!     ' Acond has exceeded Acondlim                            '   %  7
!!!     ' The iteration limit was reached                        '   %  8
!!!     ' Least-squares problem but no converged solution yet    ']; %  9

    first = 'Enter csminresqlp.  '
    last  = 'Exit csminresqlp.  '
    msg = [' beta2 = 0.  b and x are eigenvectors                   ', &  ! -1
           ' beta1 = 0.  The exact solution is  x = 0               ', &  !  0
           ' A solution to Ax = b found, given rtol                 ', &  !  1
           ' Min-length solution for singular LS problem, given rtol', &  !  2
           ' A solution to Ax = b found, given eps                  ', &  !  3
           ' Min-length solution for singular LS problem, given eps ', &  !  4
           ' x has converged to an eigenvector                      ', &  !  5
           ' xnorm has exceeded maxxnorm                            ', &  !  6
           ' Acond has exceeded Acondlim                            ', &  !  7
           ' The iteration limit was reached                        ', &  !  8
           ' Least-squares problem but no converged solution yet    ']    !  9

!!!if show
!!!  fprintf('\n%s%s', first, 'Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||')
!!!  fprintf('\nn      =%7g   ||b||    =%10.3e   shift    =%10.3e   rtol     =%10.3e',...
!!!             n,            beta1,             shift,             rtol)
!!!  fprintf('\nmaxit  =%7g   maxxnorm =%10.3e   Acondlim =%10.3e   TranCond =%10.3e',...
!!!             maxit,        maxxnorm,          Acondlim,          TranCond)
!!!  fprintf('\nprecon =%7g\n' , precon)
!!!  head = '    iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm';
!!!  fprintf('\n%s\n', head)
!!!end

  if (this%show) then
    write(*,'(2a)') first, 'Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||'
    !write(*,'("n      =",i7,"   ||b||    =",es10.3,"   shift    =",a,"   rtol     =",es10.3)') &
    !    this%n, beta1, num2str(shift), this%rtol
    write(*,'("||b||    =",es10.3,"   shift    =",a,"   rtol     =",es10.3)') &
        beta1, num2str(shift), this%rtol
    write(*,'("maxit  =",i7,"   maxxnorm =",es10.3,"   Acondlim =",es10.3,"   TranCond =",es10.3)') &
        this%maxit, this%maxxnorm, this%Acondlim, this%TranCond
    head = '    iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm'
    write(*,'(a)') head
  end if

!!!if beta1==0, flag = 0; end   % b = 0 => x = 0.  We will skip the main loop.

  if (beta1 == 0) then
    flag = 0
    !x = 0
    return
  end if


!!!%% Main iteration
!!!while flag == flag0 && iter < maxit

  do while (flag == flag0 .and. iter < this%maxit)

!!!%% Lanczos
!!!  iter  = iter + 1;
!!!  betal = beta;        beta = betan;
!!!  v = r3*(1/beta);
!!!  if (~isComplexSymmetric)
!!!      r3 = minresxxxA(A,v,inputs{:});
!!!  else
!!!      r3 = minresxxxA(A,conj(v),inputs{:});
!!!  end
!!!  if shift ~= 0
!!!      if (~isComplexSymmetric)
!!!          r3 = r3 - shift*v;
!!!      else
!!!          r3 = r3 - shift*conj(v);
!!!      end
!!!  end
!!!  if iter  >  1,       r3 = r3 - (beta/betal)*r1; end
!!!  if (~isComplexSymmetric)
!!!    alfa = real(v' * r3);  %%%%% Allow for Hermitian A.  Must get real alfa here.
!!!  else
!!!    alfa = v' * r3;        % a complex number, v' * r3 != r3' * v
!!!  end
!!!  r3   = r3 - (alfa/beta)*r2;     r1 = r2;     r2 = r3;

    iter = iter + 1
    betal = beta; beta = betan
    !v(1:n) = r3(1:n)*(1/beta)
    call v%update(cmplx(1/beta,kind=r8), r3, zzero)
    call vbar%conjg(v)
    call lin_op%matvec(vbar, r3)
    !if (shift /= 0.0_r8) r3(1:n) = r3(1:n) - shift*conjg(v(1:n))
    if (shift /= 0.0_r8) call r3%update(-shift, vbar)
    !if (iter > 1) r3(1:n) = r3(1:n) - (beta/betal)*r1(1:n)
    if (iter > 1) call r3%update(cmplx(-beta/betal,kind=r8), r1)
    !alfa = this%global_dotc(v, r3)
    alfa = v%dotc(r3)
    !r3(1:n) = r3(1:n) - (alfa/beta)*r2(1:n); r1(1:n) = r2(1:n); r2(1:n) = r3(1:n)
    call r3%update(-alfa/beta, r2)
    call r1%copy(r2)
    call r2%copy(r3)

!!!  if ~precon
!!!    betan = norm(r3);
!!!    if iter == 1       % Something special can happen
!!!      if betan == 0    % beta2 = 0
!!!        if alfa == 0   % alfa1 = 0
!!!          flag = 0;    % Ab = 0 and x = 0         ("A" = (A - shift*I))
!!!          break
!!!        else
!!!          flag = -1;   % Ab = alfa1 b,  x = b/alfa1, an eigenvector
!!!            if (~isComplexSymmetric)
!!!              x    = full(b)/alfa;
!!!            else
!!!              x    = full(conj(b))/alfa;
!!!            end
!!!          break
!!!        end
!!!      end
!!!    end
!!!  else
!!!    r3 = minresxxxM(M, M1, M2, r2,inputs{:});     betan = r2'*r3;
!!!    if betan > 0
!!!      betan = sqrt(betan);
!!!    else
!!!      error('"M" appears to be indefinite or singular.');
!!!    end
!!!  end
!!!  pnorm  = norm([betal alfa betan]);

  call lin_op%precon(r2, r3)
  !betan = this%global_dotc(r2, r3) ! dotc value should be real
  betan = r2%dotc(r3) ! dotc value should be real
  if (betan < 0.0_r8) then
    this%flag = -3  ! previously unused value
    return
  end if
  betan = sqrt(betan)
  pnorm = sqrt(betal**2 + abs(alfa)**2 + betan**2)

!!!  if debug
!!!    fprintf('\n\nLanczos iteration %d :\n', iter);
!!!    fprintf('\n  v_%d     = ', iter );  fprintf('%s ', num2str(v(1:min(n,5))') );
!!!    fprintf('\n  r1_%d    = ', iter );  fprintf('%s ', num2str(r1(1:min(n,5))') );
!!!    fprintf('\n  r2_%d    = ', iter );  fprintf('%s ', num2str(r2(1:min(n,5))') );
!!!    fprintf('\n  r3_%d    = ', iter );  fprintf('%s ', num2str(r3(1:min(n,5))') );
!!!    fprintf('\n  alpha_%d = %s, beta_%d = %s, beta_%d = %s pnorm_%d = %s ',...
!!!      iter, num2str(alfa), iter, num2str(beta), iter+1, num2str(betan), iter, num2str(pnorm) );
!!!  end

  if (debug) then
    write(*,'(//,"Lanczos iteration ",i0," :",/)') iter
    !write(*,'("v_",i0,"     = ",*(a,:,", "))') iter, (num2str(v(j)), j=1,min(this%n,5))
    !write(*,'("r1_",i0,"    = ",*(a,:,", "))') iter, (num2str(r1(j)), j=1,min(this%n,5))
    !write(*,'("r2_",i0,"    = ",*(a,:,", "))') iter, (num2str(r2(j)), j=1,min(this%n,5))
    !write(*,'("r3_",i0,"    = ",*(a,:,", "))') iter, (num2str(r3(j)), j=1,min(this%n,5))
    write(*,'("alpha_",i0," = ",a, ", beta_",i0," = ",a,", beta_",i0," = ",a,", pnorm_",i0," = ",a)') &
      iter, num2str(alfa), iter, num2str(beta), iter+1, num2str(betan), iter, num2str(pnorm)
  end if

!!!%% Apply previous left reflection Q_{k-1}
!!!  dbar  = dltan;
!!!  dlta  = cs*dbar       + sn*alfa;    epln     = eplnn;
!!!  gbar  = conj(sn)*dbar - cs*alfa;    eplnn    = sn*betan;
!!!  dltan = -cs*betan;            dlta_QLP = dlta;

    dbar  = dltan
    dlta  = cs*dbar       + sn*alfa;    epln     = eplnn
    gbar  = conjg(sn)*dbar - cs*alfa;    eplnn    = sn*betan
    dltan = -cs*betan;            dlta_QLP = dlta

!!!  if debug
!!!    fprintf('\n\nApply previous left reflection Q_{%d,%d}:\n', iter-1, iter');
!!!    fprintf('\n  c_%d     = %s, s_%d    = %s', ...
!!!      iter-1, num2str(cs), iter-1, num2str(sn) );
!!!    fprintf('\n  dlta_%d = %s, gbar_%d = %s', ...
!!!      iter, num2str(dlta), iter, num2str(gbar) );
!!!    fprintf('\n  epln_%d = %s, dbar_%d = %s', ...
!!!      iter+1, num2str(eplnn), iter+1, num2str(dltan) );
!!!  end

    if (debug) then
      write(*,'(/,"Apply previous left reflection Q_{",i0,",",i0,"}:",/)') iter-1, iter
      write(*,'("  c_",i0,"     = ",a,", s_",i0,"    = ",a)') &
        iter-1, num2str(cs), iter-1, num2str(sn)
      write(*,'("  dlta_",i0," = ",a,", gbar_",i0," = ",a)') &
        iter, num2str(dlta), iter, num2str(gbar)
      write(*,'("  epln_",i0," = ",a,", dbar_",i0," = ",a)') &
        iter+1, num2str(eplnn), iter+1, num2str(dltan)
    end if

!!!%% Compute the current left reflection Q_k
!!!  gamal3 = gamal2;     gamal2 = gamal;     gamal    = gama;
!!!  [cs,sn,gama] = SymOrtho(gbar, betan);     gama_tmp = gama;
!!!  taul2  = taul;       taul   = tau;       tau      = cs      *phi;
!!!  Axnorm = norm([Axnorm tau]);             phi      = conj(sn)*phi;

    gamal3 = gamal2; gamal2 = gamal; gamal = gama
    call SymOrtho(gbar, cmplx(betan,kind=r8), cs, sn, gama); gama_tmp = gama
    taul2 = taul; taul = tau; tau = cs*phi
    Axnorm = sqrt(Axnorm**2 + abs(tau)**2); phi = conjg(sn)*phi

!!!  if debug
!!!    fprintf('\n\nCompute the current left reflection Q_{%d,%d}:\n', iter, iter+1      );
!!!    fprintf('\n  c_%d     = %s, s_%d    = %s ', iter, num2str(cs), iter, num2str(sn)  );
!!!    fprintf('\n  tau_%d   = %s, phi_%d  = %s ', iter, num2str(tau), iter, num2str(phi));
!!!    fprintf('\n  gama_%d = %s ', iter, num2str(gama)  );
!!!  end

    if (debug) then
      write(*,'(/,"Compute the current left reflection Q_{",i0,",",i0,"}:",/)') iter, iter+1
      write(*,'("  c_",i0,"     = ",a,", s_",i0,"    = ",a)') iter, num2str(cs), iter, num2str(sn)
      write(*,'("  tau_",i0,"   = ",a,", phi_",i0,"  = ",a)') iter, num2str(tau), iter, num2str(phi)
      write(*,'("  gama_",i0," = ",a)') iter, num2str(gama)
    end if

!!!%% Apply the previous right reflection P{k-2,k}
!!!  if iter > 2
!!!    veplnl2  = veplnl;     etal2 = etal;     etal = eta;
!!!    dlta_tmp = sr2*vepln - cr2*dlta;
!!!    veplnl   = cr2*vepln + conj(sr2)*dlta;
!!!    dlta     = dlta_tmp;   eta = conj(sr2)*gama;   gama = -cr2*gama;
!!!  end

    if (iter > 2) then
      veplnl2 = veplnl; etal2 = etal; etal = eta
      dlta_tmp = sr2*vepln - cr2*dlta
      veplnl   = cr2*vepln + conjg(sr2)*dlta
      dlta = dlta_tmp; eta = conjg(sr2)*gama; gama = -cr2*gama
    end if

!!!  if debug
!!!    fprintf('\n\nApply the previous right reflections P_{%d,%d}:\n', iter-2, iter)
!!!    fprintf('\n  cr2_%d   = %s, sr2_%d    = %s', ...
!!!      iter, num2str(cr2), iter, num2str(sr2) );
!!!    fprintf('\n  gama_%d = %s, gama_%d  = %s, gama_%d = %s', ...
!!!      iter-2, num2str(gamal2), iter-1, num2str(gamal), iter, num2str(gama));
!!!    fprintf('\n  dlta_%d = %s, vepln_%d = %s, eta_%d   = %s', ...
!!!      iter, num2str(dlta), iter-1, num2str(veplnl), iter, num2str(eta) );
!!!  end

    if (debug) then
      write(*,'(/,"Apply the previous right reflections P_{",i0,",",i0,"}:",/)') iter-2, iter
      write(*,'("  cr2_",i0,"   = ",a,", sr2_",i0,"    = ",a)') &
        iter, num2str(cr2), iter, num2str(sr2)
      write(*,'("  gama_",i0," = ",a,", gama_",i0,"  = ",a,", gama_",i0," = ",a)') &
        iter-2, num2str(gamal2), iter-1, num2str(gamal), iter, num2str(gama)
      write(*,'("  dlta_",i0," = ",a,", vepln_",i0," = ",a,", eta_",i0,"   = ",a)') &
        iter, num2str(dlta), iter-1, num2str(veplnl), iter, num2str(eta)
    end if

!!!%% Compute the current right reflection P{k-1,k}, P_12, P_23,...
!!!  if iter > 1
!!!    [cr1, sr1, gamal] = SymOrtho(conj(gamal), conj(dlta)); gamal = conj(gamal);
!!!    vepln =   conj(sr1)*gama;
!!!    gama  = - cr1*gama;
!!!    if debug
!!!       fprintf('\n\nCompute the second current right reflections P_{%d,%d}:\n', iter-1, iter');
!!!       fprintf('\n  cr1_%d   = %s, sr1_%d   = %s', ...
!!!                    iter, num2str(cr1), iter, num2str(sr1) );
!!!       fprintf('\n  gama_%d = %s, gama_%d = %s, vepln_%d = %s',...
!!!                    iter-1, num2str(gamal), iter, num2str(gama), iter, num2str(vepln) );
!!!      end
!!!  end

    if (iter > 1) then
      call SymOrtho(conjg(gamal), conjg(dlta), cr1, sr1, gamal); gamal = conjg(gamal)
      vepln = conjg(sr1)*gama
      gama  = - cr1*gama
      if (debug) then
         write(*,'(/,"Compute the second current right reflections P_{",i0,",",i0,"}:",/)') iter-1, iter
         write(*,'("  cr1_",i0,"   = ",a,", sr1_",i0,"   = ",a)') &
                    iter, num2str(cr1), iter, num2str(sr1)
         write(*,'("  gama_",i0," = ",a,", gama_",i0," = ",a,", vepln_",i0," = ",a)') &
                      iter-1, num2str(gamal), iter, num2str(gama), iter, num2str(vepln)
      end if
    end if

!!!%% Update xnorm
!!!  xnorml = xnorm;     ul4 = ul3;     ul3   = ul2;
!!!  if iter > 2
!!!    ul2 = (taul2 - etal2*ul4 - veplnl2*ul3) / gamal2;
!!!  end
!!!  if iter > 1
!!!    ul = ( taul  - etal *ul3 - veplnl *ul2) / gamal;
!!!  end
!!!  xnorm_tmp = norm([xl2norm ul2 ul]);
!!!  if abs(gama) > pnorm * realmin && xnorm_tmp < maxxnorm
!!!    u = (tau - eta*ul2 - vepln*ul) / gama;
!!!    if norm([xnorm_tmp u]) > maxxnorm
!!!      u = 0;      flag = 6;
!!!    end
!!!  else
!!!    u = 0;     flag = 9;
!!!  end
!!!  xl2norm = norm([xl2norm ul2]);
!!!  xnorm   = norm([xl2norm ul u]);


  xnorml = xnorm; ul4 = ul3; ul3 = ul2
  if (iter > 2) ul2 = (taul2 - etal2*ul4 - veplnl2*ul3) / gamal2
  if (iter > 1) ul = ( taul  - etal *ul3 - veplnl *ul2) / gamal
  xnorm_tmp = sqrt(xl2norm**2 + abs(ul2)**2 + abs(ul)**2)
  if (abs(gama) > pnorm * tiny(1.0_r8) .and. xnorm_tmp < this%maxxnorm) then
    u = (tau - eta*ul2 - vepln*ul) / gama
    if (sqrt(xnorm_tmp**2 + abs(u)**2) > this%maxxnorm) then
      u = 0; flag = 6
    end if
  else
    u = 0; flag = 9
  end if
  xl2norm = sqrt(xl2norm**2 + abs(ul2)**2)
  xnorm   = sqrt(xl2norm**2 + abs(ul)**2 + abs(u)**2)

!!!%% Update w. Update x except if it will become too big
!!!  if (Acond < TranCond) && flag ~= flag0 && QLPiter==0  %% MINRES updates
!!!    wl2 = wl;     wl = w;
!!!    w   = (conj(v) - epln*wl2 - dlta_QLP*wl) * (1/gama_tmp);
!!!    if xnorm < maxxnorm
!!!      x = x + tau*w;
!!!    else
!!!      flag = 6;
!!!    end
!!!
!!!  else %% CS-MINRES-QLP updates
!!!    QLPiter = QLPiter + 1;
!!!    if QLPiter == 1
!!!      xl2 = zeros(n,1);
!!!      if  (iter > 1) % construct w_{k-3}, w_{k-2}, w_{k-1}
!!!	      if iter > 3
!!!	        wl2 = gamal3*wl2 + veplnl2*wl + etal*w;
!!!	      end     % w_{k-3}
!!!	      if iter > 2
!!!	        wl = gamal_QLP*wl + vepln_QLP*w;
!!!	      end     % w_{k-2}
!!!	      w = gama_QLP*w;     xl2 = x - wl*ul_QLP - w*u_QLP;
!!!      end
!!!    end
!!!    if iter == 1
!!!      wl2 = wl;      wl = conj(v)*conj(sr1);     w  = conj(v)*cr1;
!!!    elseif iter == 2
!!!      wl2 = wl;
!!!      wl  = w*cr1 + conj(v)*conj(sr1);
!!!      w   = w*sr1 - conj(v)*cr1;
!!!    else
!!!      wl2 = wl;      wl = w;               w  = wl2*sr2 - conj(v)*cr2;
!!!      wl2 = wl2*cr2  + conj(v)*conj(sr2);  v  = wl *cr1 + w*conj(sr1);
!!!      w   = wl*sr1 - w*cr1;                wl = v;
!!!    end
!!!    xl2 = xl2 + wl2*ul2;
!!!    x   = xl2 + wl*ul + w*u;
!!!  end

  if (Acond < this%TranCond .and. flag /= flag0 .and. QLPiter==0) then ! MINRES updates
    !wl2(1:n) = wl(1:n); wl(1:n) = w(1:n)
    call wl2%copy(wl)
    call wl%copy(w)
    !w(1:n) = (conjg(v(1:n)) - epln*wl2(1:n) - dlta_QLP*wl(1:n)) * (1/gama_tmp)
    call w%copy(vbar)
    call w%update(-epln, wl2, -dlta_QLP, wl)
    call w%scale(1/gama_tmp)
    if (xnorm < this%maxxnorm) then
      !x(1:n) = x(1:n) + tau*w(1:n)
      call x%update(tau, w)
    else
      flag = 6
    end if

  else ! CS-MINRES-QLP updates
    QLPiter = QLPiter + 1
    if (QLPiter == 1) then
      !xl2(1:n) = 0.0_r8
      call xl2%setzero
      if (iter > 1) then ! construct w_{k-3}, w_{k-2}, w_{k-1}
	!if (iter > 3) wl2(1:n) = gamal3*wl2(1:n) + veplnl2*wl(1:n) + etal*w(1:n) ! w_{k-3}
	if (iter > 3) call wl2%update(veplnl2, wl, etal, w, gamal3) ! w_{k-3}
	!if (iter > 2) wl(1:n) = gamal_QLP*wl(1:n) + vepln_QLP*w(1:n) ! w_{k-2}
	if (iter > 2) call wl%update(vepln_QLP, w, gamal_QLP) ! w_{k-2}
	!w(1:n) = gama_QLP*w(1:n); xl2(1:n) = x(1:n) - wl(1:n)*ul_QLP - w(1:n)*u_QLP
	call w%scale(gama_QLP)
        call xl2%copy(x)
        call xl2%update(-ul_QLP, wl, -u_QLP, w)
      end if
    end if
    if (iter == 1) then
      !wl2(1:n) = wl(1:n);      wl(1:n) = conjg(v(1:n))*conjg(sr1);     w(1:n)  = conjg(v(1:n))*cr1
      call wl2%copy(wl)
      call wl%update(conjg(sr1), vbar, zzero) ! lin comb
      call w%update(cmplx(cr1,kind=r8), vbar, zzero) ! lin comb
    else if (iter == 2) then
      !wl2(1:n) = wl(1:n);
      call wl2%copy(wl)
      !wl(1:n)  = w(1:n)*cr1 + conjg(v(1:n))*conjg(sr1)
      call wl%update(cmplx(cr1,kind=r8), w, conjg(sr1), vbar, zzero) ! lin comb
      !w(1:n)   = w(1:n)*sr1 - conjg(v(1:n))*cr1
      call w%update(cmplx(-cr1,kind=r8), vbar, sr1)
    else
      !wl2(1:n) = wl(1:n);      wl(1:n) = w(1:n);               w(1:n)  = wl2(1:n)*sr2 - conjg(v(1:n))*cr2
      call wl2%copy(wl)
      call wl%copy(w)
      call w%update(sr2, wl2, cmplx(-cr2,kind=r8), vbar, zzero) ! lin comb
      !wl2(1:n) = wl2(1:n)*cr2  + conjg(v(1:n))*conjg(sr2);  v(1:n)  = wl(1:n) *cr1 + w(1:n)*conjg(sr1)
      call wl2%update(conjg(sr2), vbar, cmplx(cr2,kind=r8))
      call v%update(cmplx(cr1,kind=r8), wl, conjg(sr1), w, zzero) ! lin comb
      call vbar%conjg(v)
      !w(1:n)   = wl(1:n)*sr1 - w(1:n)*cr1;                wl(1:n) = v(1:n)
      call w%update(sr1, wl, cmplx(-cr1,kind=r8))
      call wl%copy(v)
    end if
    !xl2(1:n) = xl2(1:n) + wl2(1:n)*ul2
    call xl2%update(ul2, wl2)
    !x(1:n)   = xl2(1:n) + wl(1:n)*ul + w(1:n)*u
    call x%copy(xl2)
    call x%update(ul, wl, u, w)
  end if

!!!  if debug
!!!    fprintf('\n\nUpdate w:\n');
!!!    fprintf('\n  w_%d     = ', iter-1 ); fprintf('%s ', num2str(wl(1:min(n,5))') );
!!!    fprintf('\n  w_%d     = ', iter );   fprintf('%s ', num2str(w(1:min(n,5))') );
!!!    fprintf('\n\nUpdate u, x and xnorm:\n');
!!!    fprintf('\n  u_%d     = %s, u_%d     = %s, u_%d     = %s',...
!!!      iter-2, num2str(ul2), iter-1, num2str(ul), iter, num2str(u) );
!!!    fprintf('\n  x_%d     = ', iter ); fprintf('%s ', num2str(x(1:min(n,5))') );
!!!    fprintf('\n  ||x_%d|| = ', iter ); fprintf('%s ', num2str(xnorm) );
!!!  end

  if (debug) then
    write(*,'(/,"Update w:",/)')
    !write(*,'("  w_",i0,"     = ",*(a,:,", "))') iter-1, (num2str(wl(j)), j=1,min(this%n,5))
    !write(*,'("  w_",i0,"     = ",*(a,:,", "))') iter, (num2str(w(j)), j=1,min(this%n,5))
    write(*,'(/,"Update u, x and xnorm:",/)')
    write(*,'("  u_",i0,"     = ",a,", u_",i0,"     = ",a,", u_",i0,"     = ",a)') &
        iter-2, num2str(ul2), iter-1, num2str(ul), iter, num2str(u)
    !write(*,'("  x_",i0,"     = ",*(a,:,", "))') iter, (num2str(x(j)), j=1,min(this%n,5))
    write(*,'("  ||x_",i0,"|| = ",a)') iter, num2str(xnorm)
  end if

!!!
!!!%% Compute the next right reflection P{k-1,k+1}
!!!  gamal_tmp = gamal;
!!!  [cr2,sr2,gamal] = SymOrtho(conj(gamal),conj(eplnn));  gamal= conj(gamal);
!!!  if debug
!!!     fprintf('\n\nCompute the next right reflection P_{%d,%d}:\n', iter-1, iter+1);
!!!     fprintf('\n  cr2_%d   = %s, sr2_%d    = %s,  gama_%d = %s', ...
!!!                  iter+1, num2str(cr2), iter+1, num2str(sr2), iter-1, num2str(gamal));
!!!  end

    gamal_tmp = gamal
    call SymOrtho(conjg(gamal),conjg(eplnn),cr2,sr2,gamal);  gamal= conjg(gamal)
    if (debug) then
       write(*,'("Compute the next right reflection P_{",i0,",",i0,"}:",/)') iter-1, iter+1
       write(*,'("  cr2_",i0,"   = ",a,", sr2_",i0,"    = ",a,",  gama_",i0," = ",a)') &
                    iter+1, num2str(cr2), iter+1, num2str(sr2), iter-1, num2str(gamal)
    end if

!!!%% Store quantities for transfering from MINRES to CS-MINRES-QLP
!!!  gamal_QLP = gamal_tmp;     vepln_QLP = vepln;     gama_QLP = gama;
!!!  ul_QLP    = ul;            u_QLP     = u;

    gamal_QLP = gamal_tmp;     vepln_QLP = vepln;     gama_QLP = gama
    ul_QLP    = ul;            u_QLP     = u

!!!%% Estimate various norms
!!!  abs_gama = abs(gama);      Anorml = Anorm;
!!!  Anorm = max([Anorm, abs(gamal), abs_gama]); %pnorm
!!!  if iter == 1
!!!    gmin   = gama;    gminl = gmin;
!!!  elseif iter > 1
!!!    gminl2 = gminl;   gminl = gmin;    gmin = min([abs(gminl2), abs(gamal), abs_gama]);
!!!  end
!!!  Acondl   = Acond;     Acond   = Anorm/gmin;
!!!  rnorml   = rnorm;     relresl = relres;
!!!  if flag ~= 9,
!!!    rnorm = abs(phi);
!!!  end
!!!  relres   = rnorm / (Anorm*xnorm + beta1);
!!!  rootl    = norm([gbar; dltan]);
!!!  Arnorml  = rnorml*rootl;
!!!  relAresl = rootl / Anorm;
!!!  if debug
!!!     fprintf('\n\nUpdate other norms:\n');
!!!     fprintf('\n  gmin_%d   = ', iter );   fprintf('%s ', num2str(gmin));
!!!     fprintf('\n  pnorm_%d  = ', iter );   fprintf('%s ', num2str(pnorm));
!!!     fprintf('\n  rnorm_%d  = ', iter );   fprintf('%s ', num2str(rnorm));
!!!     fprintf('\n  Arnorm_%d = ', iter-1);  fprintf('%s ', num2str(Arnorml));
!!!     fprintf('\n  Acond_%d  = ', iter );   fprintf('%s ', num2str(Acond));
!!!     fprintf('\n\n');
!!!   end

    abs_gama = abs(gama); Anorml = Anorm
    Anorm = max(Anorm, abs(gamal), abs_gama) ! pnorm
    if (iter == 1) then
      gmin   = abs_gama;    gminl = gmin
    else if (iter > 1) then
      gminl2 = gminl;   gminl = gmin;    gmin = min(abs(gminl2), abs(gamal), abs_gama)
    end if
    Acondl   = Acond;     Acond   = Anorm/gmin
    rnorml   = rnorm;     relresl = relres
    if (flag /= 9) rnorm = abs(phi)
    relres   = rnorm / (Anorm*xnorm + beta1)
    rootl    = sqrt(abs(gbar)**2 + abs(dltan)**2)
    Arnorml  = rnorml*rootl
    relAresl = rootl / Anorm
    if (debug) then
       write(*,'(/,"Update other norms:",/)')
       write(*,'("  gmin_",i0,"   = ",a)') iter, num2str(gmin)
       write(*,'("  pnorm_",i0,"  = ",a)') iter, num2str(pnorm)
       write(*,'("  rnorm_",i0,"  = ",a)') iter, num2str(rnorm)
       write(*,'("  Arnorm_",i0," = ",a)') iter-1, num2str(Arnorml)
       write(*,'("  Acond_",i0,"  = ",a)') iter, num2str(Acond)
       write(*,*)
    end if

!!!%% See if any of the stopping criteria are satisfied.
!!!  epsx = Anorm*xnorm*eps;
!!!  if (flag == flag0) || (flag == 9)
!!!    t1 = 1 + relres;
!!!    t2 = 1 + relAresl;
!!!    if iter     >= maxit   , flag = 8; end  % Too many itns
!!!    if Acond    >= Acondlim, flag = 7; end  % Huge Acond
!!!    if xnorm    >= maxxnorm, flag = 6; end  % xnorm exceeded its limit
!!!    if epsx     >= beta1   , flag = 5; end  % x is an eigenvector
!!!    if t2       <= 1       , flag = 4; end  % Accurate LS solution
!!!    if t1       <= 1       , flag = 3; end  % Accurate Ax=b solution
!!!    if relAresl <= rtol    , flag = 2; end  % Good enough LS solution
!!!    if relres   <= rtol    , flag = 1; end  % Good enough Ax=b solution
!!!  end

    epsx = Anorm*xnorm*epsilon(1.0_r8)
    if (flag == flag0 .or. flag == 9) then
      t1 = 1 + relres
      t2 = 1 + relAresl
      if (iter     >= this%maxit)    flag = 8 ! Too many itns
      if (Acond    >= this%Acondlim) flag = 7 ! Huge Acond
      if (xnorm    >= this%maxxnorm) flag = 6 ! xnorm exceeded its limit
      if (epsx     >= beta1)    flag = 5 ! x is an eigenvector
      if (t2       <= 1)        flag = 4 ! Accurate LS solution
      if (t1       <= 1)        flag = 3 ! Accurate Ax=b solution
      if (relAresl <= this%rtol)     flag = 2 ! Good enough LS solution
      if (relres   <= this%rtol)     flag = 1 ! Good enough Ax=b solution
    end if

!!!% The "disable" option allowed iterations to continue until xnorm
!!!% became large and x was effectively a nullvector.
!!!% We know that r will become a nullvector much sooner,
!!!% so we now disable the disable option :)
!!!
!!!%  if disable && (iter < maxit)
!!!%    flag = 0;
!!!%    if Axnorm < rtol*Anorm*xnorm
!!!%      flag = 10;
!!!%    end
!!!%  end

!!!  if flag == 2 || flag == 4 || (flag == 6 && likeLS) || flag == 7   % Possibly singular
!!!    iter  = iter - 1;
!!!    Acond = Acondl;   rnorm = rnorml;   relres = relresl;
!!!  else
!!!    if ~isempty(resvec)
!!!      resvec(iter+1) = rnorm;
!!!      Aresvec(iter)  = Arnorml;
!!!    end
!!!
!!!    if show && mod(iter-1,lines) == 0
!!!      if iter == 101
!!!        lines = 10;     headlines = 20*lines;
!!!      elseif iter == 1001
!!!        lines = 100;    headlines = 20*lines;
!!!      end
!!!      if QLPiter == 1
!!!        fprintf('%s', 'P')
!!!      else
!!!        fprintf('%s', ' ')
!!!      end
!!!      fprintf('%7g %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n',...
!!!        iter-1, rnorml, Arnorml, relresl, relAresl, Anorml, Acondl, xnorml)
!!!      if iter > 1 && mod(iter,headlines) == 1
!!!        fprintf('\n%s\n', head)
!!!      end
!!!    end
!!!  end
!!!end % while

    ! likeLS is undefined in the matlab code -- so I'll assume it was FALSE
    !if (flag == 2 .or. flag == 4 .or. (flag == 6 .and. likeLS) .or. flag == 7) then ! Possibly singular
    if (flag == 2 .or. flag == 4 .or. flag == 7) then ! Possibly singular
      iter  = iter - 1
      Acond = Acondl;   rnorm = rnorml;   relres = relresl
    else
      resvec(iter+1) = rnorm
      Aresvec(iter)  = Arnorml

      if (this%show .and. mod(iter-1,lines) == 0) then
        if (iter == 101) then
          lines = 10;     headlines = 20*lines
        else if (iter == 1001) then
          lines = 100;    headlines = 20*lines
        end if
        write(*,'(a,i7,*(es11.2))') merge('P',' ',QLPiter==1), iter-1, rnorml, Arnorml, relresl, relAresl, Anorml, Acondl, xnorml
        if (iter > 1 .and. modulo(iter,headlines) == 1) write(*,'(/,a)') head
      end if
    end if
  end do


!!!%% We have exited the main loop.
!!!if QLPiter == 1
!!!  fprintf('%s', 'P')
!!!else
!!!  fprintf('%s', ' ')
!!!end
!!!Miter = iter - QLPiter;

    if (this%show) write(*,'(a)',advance='no') merge('P', ' ', QLPiter==1)
    Miter = iter - QLPiter;

!!!%% Compute final quantities directly.
!!!r1      = b - minresxxxA(A,x,inputs{:}) + shift*x;   % r1 is workspace for residual vector
!!!rnorm   = norm(r1);
!!!Arnorm  = norm(minresxxxA(A,r1,inputs{:}) - shift*r1);
!!!xnorm   = norm(x);
!!!relres  = rnorm / (Anorm*xnorm + beta1);
!!!relAres = 0;
!!!if rnorm > realmin
!!!  relAres = Arnorm / (Anorm*rnorm);
!!!end

    call lin_op%matvec(x, r1)
    !r1(1:n) = b(1:n) - r1(1:n) + shift*x(1:n)   ! r1 is workspace for residual vector
    call r1%update((1.0_r8,0.0_r8), b, shift, x, (-1.0_r8,0.0_r8))   ! r1 is workspace for residual vector
    !rnorm = this%global_nrm2(r1)
    rnorm = r1%norm2()
    block
      !complex(r8) :: r2(size(r1))
      class(zvector), allocatable :: r2
      call r1%clone(r2)
      call lin_op%matvec(r1, r2)
      !Arnorm  = this%global_nrm2(r2(1:n) - shift*r1(1:n))
      call r2%update(-shift, r1)
      Arnorm = r2%norm2()
    end block
    !xnorm   = this%global_nrm2(x)
    xnorm   = x%norm2()
    relres  = rnorm / (Anorm*xnorm + beta1)
    relAres = 0
    if (rnorm > tiny(1.0_r8)) relAres = Arnorm / (Anorm*rnorm)

!!!if ~isempty(Aresvec)
!!!  Aresvec(iter+1) = Arnorm;
!!!  Aresvec         = Aresvec(1:iter+1);
!!!end
!!!if ~isempty(resvec)
!!!  resvec          = resvec(1:iter+1);
!!!end

    Aresvec(iter+1) = Arnorm
    Aresvec = Aresvec(1:iter+1)
    resvec = resvec(1:iter+1)

!!!if show
!!!  if rnorm > realmin
!!!    fprintf('%7g %10.2e %10.2eD%10.2e %10.2eD%10.2e %10.2e %10.2e\n\n',...
!!!            iter, rnorm, Arnorm, relres, relAres, Anorm, Acond, xnorm)
!!!  else
!!!    fprintf('%7g %10.2e %10.2eD%10.2e %10.2e             %10.2e %10.2e\n\n',...
!!!            iter, rnorm, Arnorm, relres,          Anorm, Acond, xnorm)
!!!  end
!!!
!!!  fprintf('\n')
!!!  fprintf('%s flag  =%7g  %s\n'                                , last, flag , msg(flag+2,:))
!!!  fprintf('%s iter  =%7g   (CS-MINRES%7g, CS-MINRES-QLP%7g)\n' , last, iter , Miter, QLPiter)
!!!  fprintf('%s rnorm = %11.4e     rnorm  direct = %11.4e\n'     , last, rnorm, norm(r1))
!!!  fprintf('%s                         Arnorm direct = %11.4e\n', last, Arnorm)
!!!  fprintf('%s xnorm = %11.4e     xnorm  direct = %11.4e\n'     , last, xnorm, norm(x))
!!!  fprintf('%s Anorm = %11.4e     Acond         = %11.4e\n'     , last, Anorm, Acond)
!!!end

    !direct_rnorm = this%global_nrm2(r1)
    direct_rnorm = r1%norm2()
    !direct_xnorm = this%global_nrm2(x)
    direct_xnorm = x%norm2()

    if (this%show) then
      if (rnorm > tiny(1.0_r8)) then
        write(*,'(i7,2(1x,es10.2),"D",es10.2,1x,es10.2,"D",es10.2,2(1x,es10.2),/)') iter, rnorm, Arnorm, relres, relAres, Anorm, Acond, xnorm
      else
        write(*,'(i7,2(1x,es10.2),"D",es10.2,1x,es10.2,"D",10x,2(1x,es10.2),/)') iter, rnorm, Arnorm, relres, Anorm, Acond, xnorm
      end if

      write(*,*)
      write(*,'(a," flag  =",i7,1x,a)') last, flag, msg(flag+2)
      write(*,'(a," iter  =",i7,"   (CS-MINRES",i7,", CS-MINRES-QLP",i7,")")') last, iter , Miter, QLPiter
      write(*,'(a," rnorm = ",es11.4,"     rnorm  direct = ",es11.4)') last, rnorm, direct_rnorm
      write(*,'(a,"                         Arnorm direct = ",es11.4)') last, Arnorm
      write(*,'(a," xnorm = ",es11.4,"     xnorm  direct = ",es11.4)') last, xnorm, direct_xnorm
      write(*,'(a," Anorm = ",es11.4,"     Acond         = ",es11.4)') last, Anorm, Acond
    end if

    this%flag = flag
    this%iter = iter
    this%Miter = Miter
    this%QLPiter = QLPiter
    this%relres = relres
    this%relAres = relAres
    this%Anorm = Anorm
    this%Acond = Acond
    this%xnorm = xnorm
    this%Axnorm = Axnorm
    this%resvec = resvec
    this%Aresvec = Aresvec

  end subroutine solve

!!!function [c, s, r] = SymOrtho(a, b)
!!!
!!!% SymOrtho: Stable Symmetric Householder reflection
!!!%
!!!%  USAGE:
!!!%     [c, s, r] = SymOrtho(a, b)
!!!%
!!!%  INPUTS:
!!!%    a      first element of a two-vector  [a; b]
!!!%    b      second element of a two-vector [a; b]
!!!%
!!!%  OUTPUTS:
!!!%    c      cosine(theta), where theta is the implicit angle of rotation
!!!%           (counter-clockwise) in a plane-rotation
!!!%    s      sine(theta)
!!!%    r      two-norm of [a; b]
!!!%
!!!%  DESCRIPTION:
!!!%     Stable symmetric Householder reflection that gives c and s such that
!!!%        [ c  s ][a] = [d],
!!!%        [ s -c ][b]   [0]
!!!%     where d = two-norm of vector [a, b],
!!!%        c = a / sqrt(a^2 + b^2) = a / d,
!!!%        s = b / sqrt(a^2 + b^2) = b / d.
!!!%     The implementation guards against overlow in computing sqrt(a^2 + b^2).
!!!%
!!!%  EXAMPLE:
!!!%      description
!!!%
!!!%  SEE ALSO:
!!!%     TESTSYMGIVENS.m,
!!!%     PLANEROT (MATLAB's function) --- 4 divisions while 2 would be enough,
!!!%     though not too time-consuming on modern machines
!!!%
!!!%  REFERENCES:
!!!%    Algorithm 4.9, stable *unsymmetric* Givens rotations in
!!!%     Golub and van Loan's book Matrix Computations, 3rd edition.
!!!%
!!!%  MODIFICATION HISTORY:
!!!%    10/06/2004: Replace d = norm([a,b]) by
!!!%                        d = a/c if |b| < |a| or b/s otherwise.
!!!%    10/07/2004: First two cases (either b or a == 0) rewritten to make sure
!!!%                (1) d >= 0
!!!%                (2) if [a,b] = 0, then c = 1 and s = 0 but not c = s = 0.
!!!%    09/27/2011: Change filename from SYMGIVENS2 to SYMORTHO.
!!!%    01/16/2012: Change file from SYMORTHO to SYMREFL.
!!!%
!!!%
!!!%  KNOWN BUGS:
!!!%     MM/DD/2004: description
!!!%
!!!%  AUTHORS: Sou-Cheng (Terrya) Choi, CI, University of Chicago
!!!%          Michael Saunders, SOL, Stanford University
!!!%
!!!%  CREATION DATE: 09/28/2004

  subroutine SymOrtho(a, b, c, s, r)
    complex(r8), intent(in) :: a, b
    real(r8), intent(out) :: c
    complex(r8), intent(out) :: s, r

!!!absa  = abs(a);
!!!absb  = abs(b);
!!!signa = sign(a);
!!!signb = sign(b);

    real(r8) :: absa, absb, t

    absa = abs(a)
    absb = abs(b)

!!!  %...........................
!!!  % Special cases: a or b is 0
!!!  %...........................
!!!  if b == 0
!!!    c = 1;
!!!    s = 0;
!!!    r = a;
!!!    return
!!!
!!!  elseif a == 0
!!!    c = 0;
!!!    s = 1;
!!!    r = b;
!!!    return
!!!  end

    if (b == 0) then
      c = 1.0_r8
      s = 0.0_r8
      r = a
    else if (a == 0) then
      c = 0.0_r8
      s = 1.0_r8
      r = b
    else

!!!  %...........................
!!!  % Both a and b are non-zero
!!!  %...........................
!!!  if absb > absa
!!!    t = absa/absb;
!!!    c = 1/sqrt(1+t^2); % temporary
!!!    s = c*conj(signb/signa);
!!!    c = c*t;
!!!    r = b/conj(s);
!!!  else
!!!    t = absb/absa;
!!!    c = 1/sqrt(1+t^2);
!!!    s = c*t*conj(signb/signa);
!!!    r = a/c;
!!!  end

      if (absb > absa) then
        t = absa/absb
        c = 1/sqrt(1+t*t) ! temporary
        s = c*conjg((b/absb)/(a/absa))
        c = c*t
        r = b/conjg(s)
      else
        t = absb/absa
        c = 1/sqrt(1+t*t)
        s = c*t*conjg((b/absb)/(a/absa))
        r = a/c
      end if
    end if
  end subroutine SymOrtho

  !! These auxiliary routines attempt to match the format produced by the
  !! the matlab num2str function. This is for the purpose of ease of
  !! comparison with the output of the matlab procedure.

  function dnum2str(x) result(s)
    real(r8), intent(in) :: x
    character(:), allocatable :: s
    character(16) :: tmp
    integer :: m, n
    if (x == 0) then
      s = '0'
    else
      write(tmp,'(g12.5)') x
      n = len_trim(tmp)
      m = scan(tmp(:n), ' ', back=.true.)
      s = tmp(m+1:n)
    end if
  end function

  function znum2str(x) result(s)
    complex(r8), intent(in) :: x
    character(:), allocatable :: s
    character(:), allocatable :: s_re, s_im
    s_re = num2str(x%re)
    if (x%im == 0) then
      s = s_re
    else
      s_im = num2str(x%im)
      if (verify(s_im(1:1), '+-') == 0) then
        s = s_re // s_im // 'i'
      else
        s = s_re // '+' // s_im // 'i'
      end if
    end if
  end function

end module cs_minres_solver2_type
