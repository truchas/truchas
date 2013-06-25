! -*- Mode: f90 -*-

module Usage_module
  implicit none

  private
  public :: &
       BasicOptions, &
       AllOptions

contains

  subroutine BasicOptions ()
    print *, '========'
    print *, 'UbikTest'
    print *, '========'
    print *
    print *, 'Usage: UbikTest [opts]'
    print *
    print *, 'Where [opts] are of the form keyword=value (no spaces around the =)'
    print *, 'and are detailed below'
    print *
    print *, 'NOTES:'
    print *, '  s indicates a character string'
    print *, '  n indicates an integer value'
    print *, '  x indicates a floating-point value'
    print *, '  defaults are shown in [] where applicable'
    print *
    print *, 'help=s'
    print *, '  basic - show this message and exit'
    print *, '  all   - show all options and exit'
    print *, '          (includes sample command lines)'
    print *, 'problem=n'
    print *, ' [0] - read system from file'
    print *, '  1  - Saad test problem, full storage'
    print *, '  2  - Saad test problem, ELL storage'
    print *, 'asymm=x'
    print *, '  diagonal-dominance factor for Saad problem [0.2]'
    print *, '    (relevant only when problem=1 or problem=2)'
    print *, 'unknowns=n'
    print *, '  number of unknowns for Saad problem [80]'
    print *, '    (relevant only when problem=1 or problem=2)'
    print *, 'blocksize=n'
    print *, '  blocksize for Saad problem [20]'
    print *, 'method=s'
    print *, '  cg      - conjugate gradients (CG)'
    print *, ' [gmres]  - generalized minimal residuals (GMRES)'
    print *, '  cgs     - conjugate gradients squared (CGS, BCGS)'
    print *, '  tfqmr   - transpose-free quasi-minimal residuals (TFQMR)'
    print *, '  bcgstab - stabilized bi-conjugate gradients (Bi-CGSTAB)'
    print *, '  sor     - successive over-relaxation (SOR)'
    print *, '  jacobi  - Jacobi'
    print *, '  lu      - LU decomposition'
    print *, 'pre=s'
    print *, ' [none] - no preconditioning'
    print *, '  jac   - m-step weighted Jacobi'
    print *, '  ssor  - m-step weighted symmetric SOR (SSOR)'
    print *, '  ic    - incomplete Cholesky (zero fill-in)'
    print *, '  ilu   - incomplete LU (zero fill-in)'
    print *, 'subspace=n'
    print *, '  Krylov subspace size for GMRES (ignored for other methods)'
    print *, '  0 implies automatic selection, and is the default'
    print *, 'steps=n'
    print *, '  number of steps to take when using Jacobi or SSOR preconditioning'
    print *, '  default = 1'
    print *, 'relax=x'
    print *, '  relaxation parameter for Jacobi, SOR, SSOR, IC, or ILU'
    print *, '  must be in the range (0,2), default is 1'
    print *, 'scale=s'
    print *, ' [none]'
    print *, '  row'
    print *, '  diagonal'
    print *, '  both - row, then diagonal'
    print *, 'eps=x'
    print *, '  convergence criterion [1.0e-6]'
    print *, 'itmax=n'
    print *, '  maximum number of iterations to allow'
    print *, '  0 implies automatic selection, and is the default'
    print *, 'stop=s'
    print *, '  relchg - || x - xold || / || x ||'
    print *, '  axb    - || r || / (|| A ||*|| x || + || b ||)'
    print *, ' [b]     - || r || / || b ||'
    print *, '  x      - || r || / || x ||'
    print *, '  r0     - || r || / || r0 ||'
    print *, '  r      - || r ||'
    print *, 'use_true_residual=s'
    print *, ' [no] - use available residual in stopping test'
    print *, '  yes - use || b - Ax || in stopping test'
    print *, 'verbosity=s'
    print *, '  none      - no output'
    print *, '  errors    - write out error messages only'
    print *, ' [warnings] - errors plus warning messages'
    print *, '  summary   - same as warnings, plus a one-line summary each time'
    print *, '              convergence is checked'
    print *, '  system    - same as summary, plus coefficient, preconditioner (if there'
    print *, '              is one), source, initial guess and converged solution'
    print *, '  iterates  - same as system, plus iterate and true residual (if computed)'
    print *, '              at each iteration'
    print *, '  full      - same as 4, plus vectors at each stage of the iteration,'
    print *, '              details about the preconditioner (if available), and more'
    print *, 'freq=n'
    print *, '  frequency at which to write an intermediate status summary to the tty'
    print *, '  0 implies no intermediate updates, and is the default'

    return
  end subroutine BasicOptions

  subroutine AllOptions ()

    call BasicOptions ()
    print *
    print *, 'norm=n'
    print *, '  0  - infinity (max. abs. value)'
    print *, '  1  - L-1 (sum of abs. values)'
    print *, ' [2] - L-2 for vectors, Frobenius for matrices (square root of sum of squares)'
    print *, 'sysfile=s'
    print *, '  name of file containing system [system.in]'
    print *, '    (relevant only when problem=0)'
    print *, 'luout=n'
    print *, '  logical unit for primary output [99]'
    print *, 'luerr=n'
    print *, '  logical unit for other errors [98]'
    print *, 'residual=n'
    print *, '  controls how often the true residual (b-Ax) is computed and used'
    print *, ' [0] - never'
    print *, '  1  - compute whenever a new iterate is obtained, but never replace the'
    print *, '       recursively-generated residual (r) with b-Ax'
    print *, '  2  - replace r with b-Ax every MAX(10,SQRT(N)) iterations, where N is'
    print *, '       the number of unknowns'
    print *, '  3  - compute b-Ax whenever a new iterate is obtained, but only replace'
    print *, '       r with b-Ax every MAX(10,SQRT(N)) iterations'
    print *, 'plot=s'
    print *, '  controls plotting of residual norm at each iteration'
    print *, ' [no] - store residual norms but do not automatically run Gnuplot'
    print *, '  yes - store residual norms and run Gnuplot to display'
    print *
    print *, 'SAMPLE COMMAND LINES:'
    print *
    print *, 'Saad problem, full format, CG:'
    print *, '% UbikTest problem=1 asymm=0 method=cg'
    print *
    print *, 'Saad problem, ELL format, GMRES:'
    print *, '% UbikTest problem=2 method=gmres'
    print *
    print *, 'Saad problem, ELL format, GMRES(5):'
    print *, '% UbikTest problem=2 method=gmres subspace=5'
    print *
    print *, 'Saad problem, ELL format, TFQMR, SSOR precond:'
    print *, '% UbikTest problem=2 method=tfqmr pre=ssor'
    print *
    print *, 'Saad problem, ELL format, TFQMR, SSOR(2) precond w/relaxation param 1.2:'
    print *, '% UbikTest problem=2 method=tfqmr pre=ssor relax=1.2 steps=2'

    return
  end subroutine AllOptions
end module Usage_module
