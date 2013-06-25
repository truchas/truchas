!!
!! FHYPRE
!!
!!  A Fortran interface to the Hypre package.
!!
!!    Markus Berndt <berndt@lanl.gov>
!!    Neil N. Carlson <nnc@lanl.gov>
!!
!! PROGRAMMING INTERFACE
!!
!!  This module provides a Fortran interface to a portion of the Hypre v2.x
!!  package from LLNL.  This interface is somewhat different than the optional
!!  Fortran interface that may be included in the build of Hypre.  It depends
!!  only on the Hypre C library and the C-source file "fhypre_c.c" that provides
!!  some necessary glue code.  Included is most of the IJMatrix/Vector system
!!  interface and the ParCSR BoomerAMG solver interface.
!!
!!  Except in a few cases, the names of the Fortran procedures are the same as
!!  their C counterpart but with an "f" prefixed.  They are also subroutines
!!  with the integer error return value returned as a new final argument.  The
!!  following procedure names were altered in order to remain within the 31
!!  character limit:
!!
!!    HYPRE IJMatrixSetMaxOffProcElmts  --> fHYPRE_IJMatrixSetMaxOffPValues
!!    HYPRE IJVectorSetMaxOffProcElmts  --> fHYPRE_IJVectorSetMaxOffPValues
!!    HYPRE_BoomerAMGSetStrongThreshold --> fHYPRE_BoomerAMGSetStrongThld
!!
!!  The routines have the same arguments with the following exceptions:
!!  * The MPI_Comm argument of the IJMatrixCreate and IJVectorCreate procedures
!!    is omitted.  The Hypre routines will implicitly use MPI_COMM_WORLD.
!!  * The vector, matrix and solver arguments, which in the C interface are
!!    pointers to C structures, are replaced by INTEGER(HYPRE_OBJ) arguments
!!    that are opaque handles.
!!
!!  In the C interface the BoomerAMGSetup and BoomerAMGSolve functions take
!!  "constructed" matrix and vector objects (as returned by the GetObject
!!  methods) as arguments.  This distinction has been eliminated in the
!!  Fortran interface.  The Setup and Solve subroutines accept the matrix
!!  and vector handles returned by the Create subroutines as arguments;
!!  the GetObject method is not needed.
!!
!!  The IJMatrixCreate and IJVectorCreate subroutines automatically create
!!  objects of HYPRE_PARCSR storage type; the SetObjectType methods are
!!  not needed.
!!
!!  A 0-valued handle for a Hypre matrix/vector/solver object signifies an
!!  undefined object.
!!
!! IMPLEMENTATION NOTES
!!
!!  We need to expose some of the error handling methods of Hypre (not all
!!  documented, btw).  Right now we limit ourselves to ==0 is good and !=0
!!  is bad.
!!

module fhypre

  use kinds, only: hypre_obj => c_intptr_t
  implicit none
  public
  
 !!
 !! IJMATRIX INTERFACE
 !!

  interface
    subroutine fHYPRE_IJMatrixCreate (ilower, iupper, jlower, jupper, matrix, ierr)
      use kinds
      integer :: ilower, iupper, jlower, jupper, ierr
      integer(c_intptr_t) :: matrix
    end subroutine
    subroutine fHYPRE_IJMatrixDestroy (matrix, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: ierr
    end subroutine
    subroutine fHYPRE_IJMatrixSetRowSizes (matrix, sizes, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: sizes(*), ierr
    end subroutine
    subroutine fHYPRE_IJMatrixSetDiagOffdSizes (matrix, diag_sizes, offd_sizes, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: diag_sizes(*), offd_sizes(*), ierr
    end subroutine
    subroutine fHYPRE_IJMatrixSetMaxOffPValues (matrix, max_offp_values, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: max_offp_values, ierr
    end subroutine
    subroutine fHYPRE_IJMatrixInitialize (matrix, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: ierr
    end subroutine
    subroutine fHYPRE_IJMatrixAssemble (matrix, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: ierr
    end subroutine
    subroutine fHYPRE_IJMatrixSetValues (matrix, nrows, ncols, rows, cols, values, ierr)
      use kinds
      integer(c_intptr_t) :: matrix
      integer :: nrows, ncols(*), rows(*), cols(*), ierr
      real(r8) :: values(*)
    end subroutine
  end interface

 !!
 !! IJVECTOR INTERFACE
 !!

  interface
    subroutine fHYPRE_IJVectorCreate (jlower, jupper, vector, ierr)
      use kinds
      integer :: jlower, jupper, ierr
      integer(c_intptr_t) :: vector
    end subroutine
    subroutine fHYPRE_IJVectorDestroy (vector, ierr)
      use kinds
      integer(c_intptr_t) :: vector
      integer :: ierr
    end subroutine
    subroutine fHYPRE_IJVectorSetMaxOffPValues (vector, max_offp_values, ierr)
      use kinds
      integer(c_intptr_t) :: vector
      integer :: max_offp_values, ierr
    end subroutine
    subroutine fHYPRE_IJVectorInitialize (vector, ierr)
      use kinds
      integer(c_intptr_t) :: vector
      integer :: ierr
    end subroutine
    subroutine fHYPRE_IJVectorSetValues (vector, nvalues, indices, values, ierr)
      use kinds
      integer(c_intptr_t) :: vector
      integer :: nvalues, indices(*), ierr
      real(r8) :: values(*)
    end subroutine
    subroutine fHYPRE_IJVectorGetValues (vector, nvalues, indices, values, ierr)
      use kinds
      integer(c_intptr_t) :: vector
      integer :: nvalues, indices(*), ierr
      real(r8) :: values(*)
    end subroutine
    subroutine fHYPRE_IJVectorAssemble (vector, ierr)
      use kinds
      integer(c_intptr_t) :: vector
      integer :: ierr
    end subroutine
  end interface

 !!
 !! BOOMERAMG INTERFACE
 !!

  interface
    subroutine fHYPRE_BoomerAMGCreate (solver, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGDestroy (solver, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetup (solver, A, b, x, ierr)
      use kinds
      integer(c_intptr_t) :: solver, A, b, x
      integer :: ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSolve (solver, A, b, x, ierr)
      use kinds
      integer(c_intptr_t) :: solver, A, b, x
      integer :: ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetStrongThld (solver, strong_threshold, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      real(r8) :: strong_threshold
      integer :: ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetMaxIter (solver, max_iter, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: max_iter, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetCoarsenType (solver, coarsen_type, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: coarsen_type, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetTol (solver, tol, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      real(r8) :: tol
      integer :: ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetNumSweeps (solver, num_sweeps, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: num_sweeps, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetRelaxType (solver, relax_type, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: relax_type, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetPrintLevel (solver, print_level, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: print_level, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetMaxLevels (solver, max_levels, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: max_levels, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetCycleNumSweeps (solver, num_sweeps, k, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: num_sweeps, k, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetCycleRelaxType (solver, relax_type, k, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: relax_type, k, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetCycleType (solver, cycle_type, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: cycle_type, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetDebugFlag (solver, debug, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: debug, ierr
    end subroutine
    subroutine fHYPRE_BoomerAMGSetLogging (solver, logging, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: logging, ierr
    end subroutine
  end interface

 !!
 !! PCG INTERFACE
 !!

  interface
    subroutine fHYPRE_PCGCreate (solver, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGDestroy (solver, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGSetup (solver, A, b, x, ierr)
      use kinds
      integer(c_intptr_t) :: solver, A, b, x
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGGSolve (solver, A, b, x, ierr)
      use kinds
      integer(c_intptr_t) :: solver, A, b, x
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGSetMaxIter (solver, max_iter, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: max_iter, ierr
    end subroutine
    subroutine fHYPRE_PCGGetNumIterations (solver, iters, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: iters, ierr
    end subroutine fHYPRE_PCGGetNumIterations
    subroutine fHYPRE_PCGSetTol (solver, tol, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      real(r8) :: tol
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGSetAbsTol (solver, tol, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      real(r8) :: tol
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGSetTwoNorm (solver, twonorm, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: twonorm
      integer :: ierr
    end subroutine
    subroutine fHYPRE_PCGGetFinalRelRes (solver, tol, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      real(r8) :: tol
      integer :: ierr
    end subroutine   
    subroutine fHYPRE_PCGSetPrintLevel (solver, print_level, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer :: print_level, ierr
    end subroutine
    subroutine fHYPRE_PCGSetPrecond (solver, precond, ierr)
      use kinds
      integer(c_intptr_t) :: solver
      integer(c_intptr_t) :: precond
      integer :: ierr
    end subroutine fHYPRE_PCGSetPrecond
  end interface

 !!
 !! UTILITIES INTERFACE
 !!

  interface
    subroutine fHYPRE_ClearAllErrors ()
    end subroutine
  end interface

end module fhypre
