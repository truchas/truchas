!!
!! FHYPRE
!!
!!  A Fortran interface to a subset of the Hypre linear solver library.
!!
!!  Neil N. Carlson <nnc@lanl.gov>
!!  March 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  This module provides a Fortran interface to a portion of the Hypre v2.x
!!  library from LLNL.  This interface is somewhat different than the optional
!!  Fortran interface that may be included in the build of Hypre.  It depends
!!  only on the Hypre C library.  Included is most of the IJMatrix/Vector system
!!  interface, and the ParCSR BoomerAMG and PCG solver interfaces.
!!
!!  Except in a few cases highlighted below, the Fortran procedures have the
!!  same name and signature as their C counterparts, except that the name is
!!  prefixed with an "f" and they are subroutines with the integer error return
!!  value returned as the final argument.  The C function arguments that are
!!  pointers to structures (usually hidden as a typedef) are replaced by
!!  TYPE(HYPRE_OBJ) arguments in the Fortran interface, and should be regarded
!!  as opaque handles.  Specific procedure differences follow.
!!
!!  * The MPI communicator argument has been omitted from IJVectorCreate,
!!    IJMatrixCreate, and ParCSRPCGCreate; MPI_COMM_WORLD will be used.
!!
!!  * IJVectorCreate and IJMatrixCreate create HYPRE_PARCSR type objects.
!!    A second call is required in the C interface to set the object type.
!!
!!  * In the C interface, the BoomerAMGSetup, BoomerAMGSolve, PCGSetup, and
!!    PCGSolve functions take "constructed" matrix and vector objects (as
!!    returned by the GetObject methods) as arguments.  This distinction has
!!    been eliminated in the Fortran interface. The Setup and Solve subroutines
!!    accept the matrix and vector handles returned by the Create subroutines
!!    as arguments; the GetObject method is not needed.
!!
!!  * PCGSetPrecond is implemented in the Fortran interface as
!!    PCGSetBoomerAMGPrecond.  The function pointer arguments have been omitted
!!    and the preconditioner is hardwired to BoomerAMG.
!!
!!  TYPE(HYPRE_OBJ) variables can be assigned the named constant HYPRE_NULL_OBJ
!!  to initialize them to a state that references no Hypre object, and the
!!  logical function HYPRE_ASSOCIATED() can be used to test whether such a
!!  a variable references a Hypre object as returned by one of the creation
!!  subroutines.
!!
!! IMPLEMENTATION NOTES
!!
!!  HYPRE_OBJ, HYPRE_NULL_OBJ, and HYPRE_ASSOCIATED are just aliases to C_PTR,
!!  C_NULL_PTR, and C_ASSOCIATED from the intrinsic ISO_C_BINDING module.
!!
!!  We need to expose some of the error handling methods of Hypre (not all
!!  documented, btw).  Right now we limit ourselves to ==0 is good and !=0
!!  is bad.
!!

#if defined(GNU_PR83148)
module fhypre_types
  use,intrinsic :: iso_c_binding, only: hypre_obj => c_ptr
  use,intrinsic :: iso_c_binding, only: hypre_null_obj => c_null_ptr
  use,intrinsic :: iso_c_binding, only: hypre_associated => c_associated
  public
end module fhypre_types
#endif

module fhypre

  use kinds, only: r8
  use hypre_c_binding
  use,intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
#if defined(GNU_PR83148)
  use fhypre_types
#else
  use,intrinsic :: iso_c_binding, only: hypre_obj => c_ptr
  use,intrinsic :: iso_c_binding, only: hypre_null_obj => c_null_ptr
  use,intrinsic :: iso_c_binding, only: hypre_associated => c_associated
#endif
  implicit none
  private
  
  !! Error codes
  public :: HYPRE_ERROR_GENERIC, HYPRE_ERROR_MEMORY, HYPRE_ERROR_ARG, HYPRE_ERROR_CONV

  !! Hide the face
  public :: hypre_obj, hypre_null_obj, hypre_associated

  !! IJVector interface procedures
  public :: fHYPRE_IJVectorCreate
  public :: fHYPRE_IJVectorDestroy
  public :: fHYPRE_IJVectorSetMaxOffProcElmts
  public :: fHYPRE_IJVectorInitialize
  public :: fHYPRE_IJVectorSetValues
  public :: fHYPRE_IJVectorAddToValues
  public :: fHYPRE_IJVectorGetValues
  public :: fHYPRE_IJVectorAssemble

  !! IJMatrix interface procedures
  public :: fHYPRE_IJMatrixCreate
  public :: fHYPRE_IJMatrixDestroy
  public :: fHYPRE_IJMatrixSetRowSizes
  public :: fHYPRE_IJMatrixSetDiagOffdSizes
  public :: fHYPRE_IJMatrixSetMaxOffProcElmts
  public :: fHYPRE_IJMatrixInitialize
  public :: fHYPRE_IJMatrixAssemble
  public :: fHYPRE_IJMatrixSetValues
  public :: fHYPRE_IJMatrixAddToValues

  !! BoomerAMG interface procedures
  public :: fHYPRE_BoomerAMGCreate
  public :: fHYPRE_BoomerAMGDestroy
  public :: fHYPRE_BoomerAMGSetup
  public :: fHYPRE_BoomerAMGSolve
  public :: fHYPRE_BoomerAMGSetStrongThreshold
  public :: fHYPRE_BoomerAMGSetMaxIter
  public :: fHYPRE_BoomerAMGSetCoarsenType
  public :: fHYPRE_BoomerAMGSetInterpType
  public :: fHYPRE_BoomerAMGSetTol
  public :: fHYPRE_BoomerAMGSetNumSweeps
  public :: fHYPRE_BoomerAMGSetRelaxType
  public :: fHYPRE_BoomerAMGSetPrintLevel
  public :: fHYPRE_BoomerAMGSetMaxLevels
  public :: fHYPRE_BoomerAMGSetCycleNumSweeps
  public :: fHYPRE_BoomerAMGSetCycleRelaxType
  public :: fHYPRE_BoomerAMGSetCycleType
  public :: fHYPRE_BoomerAMGSetDebugFlag
  public :: fHYPRE_BoomerAMGSetLogging
  public :: fHYPRE_BoomerAMGSetOldDefault

  !! PCG interface procedures
  public :: fHYPRE_PCGCreate
  public :: fHYPRE_PCGDestroy
  public :: fHYPRE_PCGSetup
  public :: fHYPRE_PCGSolve
  public :: fHYPRE_PCGSetPrecond
  public :: fHYPRE_PCGSetMaxIter
  public :: fHYPRE_PCGGetNumIterations
  public :: fHYPRE_PCGSetTol
  public :: fHYPRE_PCGSetAbsoluteTol
  public :: fHYPRE_PCGSetTwoNorm
  public :: fHYPRE_PCGGetFinalRelRes
  public :: fHYPRE_PCGSetPrintLevel

  !! ParCSR Hybrid interface procedures
  public :: fHYPRE_ParCSRHybridCreate
  public :: fHYPRE_ParCSRHybridDestroy
  public :: fHYPRE_ParCSRHybridSetup
  public :: fHYPRE_ParCSRHybridSolve
  public :: fHYPRE_ParCSRHybridSetTol
  public :: fHYPRE_ParCSRHybridSetAbsoluteTol
  public :: fHYPRE_ParCSRHybridSetConvergenceTol
  public :: fHYPRE_ParCSRHybridSetDSCGMaxIter
  public :: fHYPRE_ParCSRHybridSetPCGMaxIter
  public :: fHYPRE_ParCSRHybridSetSolverType
  public :: fHYPRE_ParCSRHybridSetKDim
  public :: fHYPRE_ParCSRHybridSetTwoNorm
  public :: fHYPRE_ParCSRHybridSetLogging
  public :: fHYPRE_ParCSRHybridSetPrintLevel
  public :: fHYPRE_ParCSRHybridSetStrongThreshold
  public :: fHYPRE_ParCSRHybridSetMaxLevels
  public :: fHYPRE_ParCSRHybridSetCoarsenType
  public :: fHYPRE_ParCSRHybridSetInterpType
  public :: fHYPRE_ParCSRHybridSetNumSweeps
  public :: fHYPRE_ParCSRHybridSetRelaxType
  public :: fHYPRE_ParCSRHybridGetNumIterations
  public :: fHYPRE_ParCSRHybridGetDSCGNumIterations
  public :: fHYPRE_ParCSRHybridGetPCGNumIterations
  public :: fHYPRE_ParCSRHybridGetFinalRelativeResidualNorm

  !! Miscellaneous procedures
  public :: fHYPRE_ClearAllErrors

contains

  !!!! IJVECTOR INTERFACE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fHYPRE_IJVectorCreate (jlower, jupper, vector, ierr)
    integer, intent(in) :: jlower, jupper
    type(c_ptr), intent(inout) :: vector
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_IJVectorCreate(jlower, jupper, vector)
  end subroutine

  subroutine fHYPRE_IJVectorDestroy (vector, ierr)
    type(c_ptr), intent(inout) :: vector
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorDestroy(vector)
    vector = c_null_ptr
  end subroutine

  subroutine fHYPRE_IJVectorSetMaxOffProcElmts (vector, max_offp_values, ierr)
    type(c_ptr), intent(in) :: vector
    integer, intent(in) :: max_offp_values
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorSetMaxOffProcElmts(vector, max_offp_values)
  end subroutine

  subroutine fHYPRE_IJVectorInitialize (vector, ierr)
    type(c_ptr), intent(in) :: vector
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorInitialize(vector)
  end subroutine

  subroutine fHYPRE_IJVectorSetValues (vector, nvalues, indices, values, ierr)
    type(c_ptr), intent(in) :: vector
    integer, intent(in) :: nvalues, indices(:)
    real(r8), intent(in) :: values(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorSetValues(vector, nvalues, indices, values)
  end subroutine

  subroutine fHYPRE_IJVectorAddToValues (vector, nvalues, indices, values, ierr)
    type(c_ptr), intent(in) :: vector
    integer, intent(in) :: nvalues, indices(:)
    real(r8), intent(in) :: values(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorAddToValues(vector, nvalues, indices, values)
  end subroutine

  subroutine fHYPRE_IJVectorGetValues (vector, nvalues, indices, values, ierr)
    type(c_ptr), intent(in) :: vector
    integer, intent(in) :: nvalues, indices(:)
    real(r8), intent(out) :: values(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorGetValues(vector, nvalues, indices, values)
  end subroutine

  subroutine fHYPRE_IJVectorAssemble (vector, ierr)
    type(c_ptr), intent(in) :: vector
    integer, intent(out) :: ierr
    ierr = HYPRE_IJVectorAssemble(vector)
  end subroutine

  !!!! IJMATRIX INTERFACE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fHYPRE_IJMatrixCreate (ilower, iupper, jlower, jupper, matrix, ierr)
    integer, intent(in) :: ilower, iupper, jlower, jupper
    type(c_ptr), intent(inout) :: matrix
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_IJMatrixCreate(ilower, iupper, jlower, jupper, matrix)
  end subroutine

  subroutine fHYPRE_IJMatrixDestroy (matrix, ierr)
    type(c_ptr), intent(inout) :: matrix
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixDestroy(matrix)
    matrix = c_null_ptr
  end subroutine

  subroutine fHYPRE_IJMatrixSetRowSizes (matrix, sizes, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(in)  :: sizes(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixSetRowSizes(matrix, sizes)
  end subroutine

  subroutine fHYPRE_IJMatrixSetDiagOffdSizes (matrix, diag_sizes, offd_sizes, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(in) :: diag_sizes(:), offd_sizes(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixSetDiagOffdSizes(matrix, diag_sizes, offd_sizes)
  end subroutine

  subroutine fHYPRE_IJMatrixSetMaxOffProcElmts (matrix, max_offp_values, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(in)  :: max_offp_values
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixSetMaxOffProcElmts(matrix, max_offp_values)
  end subroutine

  subroutine fHYPRE_IJMatrixInitialize (matrix, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixInitialize(matrix)
  end subroutine

  subroutine fHYPRE_IJMatrixAssemble (matrix, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixAssemble(matrix)
  end subroutine

  subroutine fHYPRE_IJMatrixSetValues (matrix, nrows, ncols, rows, cols, values, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(in)  :: nrows, ncols(:), rows(:), cols(:)
    real(r8), intent(in) :: values(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixSetValues(matrix, nrows, ncols, rows, cols, values)
  end subroutine

  subroutine fHYPRE_IJMatrixAddToValues (matrix, nrows, ncols, rows, cols, values, ierr)
    type(c_ptr), intent(in) :: matrix
    integer, intent(in)  :: nrows, ncols(:), rows(:), cols(:)
    real(r8), intent(in) :: values(:)
    integer, intent(out) :: ierr
    ierr = HYPRE_IJMatrixAddToValues(matrix, nrows, ncols, rows, cols, values)
  end subroutine

  !!!! BOOMER AMG INTERFACE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fHYPRE_BoomerAMGCreate (solver, ierr)
    type(c_ptr), intent(out) :: solver
    integer, intent(out) :: ierr
    ierr =  HYPRE_BoomerAMGCreate(solver)
  end subroutine

  subroutine fHYPRE_BoomerAMGDestroy (solver, ierr)
    type(c_ptr), intent(inout) :: solver
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGDestroy(solver)
    solver = c_null_ptr
  end subroutine

  subroutine fHYPRE_BoomerAMGSetup (solver, A, b, x, ierr)
    type(c_ptr), intent(in) :: solver, A, b, x
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_BoomerAMGSetup(solver, A, b, x)
  end subroutine

  subroutine fHYPRE_BoomerAMGSolve (solver, A, b, x, ierr)
    type(c_ptr), intent(in) :: solver, A, b, x
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_BoomerAMGSolve(solver, A, b, x)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetStrongThreshold (solver, strong_threshold, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: strong_threshold
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetStrongThreshold(solver, strong_threshold)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetMaxIter (solver, max_iter, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: max_iter
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetMaxIter(solver, max_iter)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetInterpType (solver, interp_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: interp_type
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetInterpType(solver, interp_type)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetCoarsenType (solver, coarsen_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: coarsen_type
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetCoarsenType(solver, coarsen_type)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetTol (solver, tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: tol
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetTol(solver, tol)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetNumSweeps (solver, num_sweeps, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: num_sweeps
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetNumSweeps(solver, num_sweeps)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetRelaxType (solver, relax_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: relax_type
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetRelaxType(solver, relax_type)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetPrintLevel (solver, print_level, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: print_level
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetPrintLevel(solver, print_level)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetMaxLevels (solver, max_levels, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: max_levels
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetMaxLevels(solver, max_levels)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetCycleNumSweeps (solver, num_sweeps, k, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: num_sweeps, k
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetCycleNumSweeps(solver, num_sweeps, k)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetCycleRelaxType (solver, relax_type, k, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: relax_type, k
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetCycleRelaxType(solver, relax_type, k)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetCycleType (solver, cycle_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: cycle_type
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetCycleType(solver, cycle_type)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetDebugFlag (solver, debug, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: debug
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetDebugFlag(solver, debug)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetLogging (solver, logging, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in)  :: logging
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetLogging(solver, logging)
  end subroutine

  subroutine fHYPRE_BoomerAMGSetOldDefault(solver, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(out) :: ierr
    ierr = HYPRE_BoomerAMGSetOldDefault(solver)
  end subroutine fHYPRE_BoomerAMGSetOldDefault

  !!!! PCG INTERFACE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fHYPRE_PCGCreate (solver, ierr)
    type(c_ptr), intent(inout) :: solver
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_ParCSRPCGCreate(solver)
  end subroutine

  subroutine fHYPRE_PCGDestroy (solver, ierr)
    type(c_ptr), intent(inout) :: solver
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGDestroy(solver)
    solver = c_null_ptr
  end subroutine

  subroutine fHYPRE_PCGSetup (solver, A, b, x, ierr)
    type(c_ptr), intent(in) :: solver, A, b, x
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_PCGSetup(solver, A, b, x)
  end subroutine

  subroutine fHYPRE_PCGSolve (solver, A, b, x, ierr)
    type(c_ptr), intent(in) :: solver, A, b, x
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_PCGSolve(solver, A, b, x)
  end subroutine

  subroutine fHYPRE_PCGSetPrecond (solver, precond, ierr)
    type(c_ptr), intent(in) :: solver, precond
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_PCGSetBoomerAMGPrecond(solver, precond)
  end subroutine fHYPRE_PCGSetPrecond

  subroutine fHYPRE_PCGSetMaxIter (solver, max_iter, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: max_iter
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGSetMaxIter(solver, max_iter)
  end subroutine

  subroutine fHYPRE_PCGGetNumIterations (solver, iters, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(out) :: iters, ierr
    ierr = HYPRE_PCGGetNumIterations(solver, iters)
  end subroutine fHYPRE_PCGGetNumIterations

  subroutine fHYPRE_PCGSetTol (solver, tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: tol
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGSetTol(solver, tol)
  end subroutine

  subroutine fHYPRE_PCGSetAbsoluteTol (solver, tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: tol
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGSetAbsoluteTol(solver, tol)
  end subroutine

  subroutine fHYPRE_PCGSetTwoNorm (solver, twonorm, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: twonorm
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGSetTwoNorm (solver, twonorm)
  end subroutine

  subroutine fHYPRE_PCGGetFinalRelRes (solver, tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(out) :: tol
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGGetFinalRelativeResidualNorm(solver, tol)
  end subroutine

  subroutine fHYPRE_PCGSetPrintLevel (solver, print_level, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: print_level
    integer, intent(out) :: ierr
    ierr = HYPRE_PCGSetPrintLevel(solver, print_level)
  end subroutine

  !!!! ParCSR Hybrid INTERFACE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fHYPRE_ParCSRHybridCreate (solver, ierr)
    type(c_ptr), intent(inout) :: solver
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridCreate (solver)
  end subroutine

  subroutine fHYPRE_ParCSRHybridDestroy (solver, ierr)
    type(c_ptr), intent(inout) :: solver
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridDestroy (solver)
    solver = c_null_ptr
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetup(solver, A, b, x, ierr)
    type(c_ptr), intent(in) :: solver, A, b, x
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_HybridSetup(solver, A, b, x)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSolve(solver, A, b, x, ierr)
    type(c_ptr), intent(in) :: solver, A, b, x
    integer, intent(out) :: ierr
    ierr = HYPRE_Ext_HybridSolve(solver, A, b, x)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetTol (solver, tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: tol
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetTol (solver, tol)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetAbsoluteTol (solver, tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: tol
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetAbsoluteTol (solver, tol)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetConvergenceTol (solver, cf_tol, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: cf_tol
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetConvergenceTol (solver, cf_tol)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetDSCGMaxIter (solver, dscg_max_its, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: dscg_max_its
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetDSCGMaxIter (solver, dscg_max_its)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetPCGMaxIter (solver, pcg_max_its, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: pcg_max_its
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetPCGMaxIter (solver, pcg_max_its)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetSolverType (solver, solver_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: solver_type
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetSolverType (solver, solver_type)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetKDim (solver, k_dim, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: k_dim
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetKDim (solver, k_dim)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetTwoNorm (solver, two_norm, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: two_norm
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetTwoNorm (solver, two_norm)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetLogging (solver, logging, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: logging
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetLogging (solver, logging)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetPrintLevel (solver, print_level, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: print_level
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetPrintLevel (solver, print_level)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetStrongThreshold (solver, strong_threshold, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(in) :: strong_threshold
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetStrongThreshold (solver, strong_threshold)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetMaxLevels (solver, max_levels, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: max_levels
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetMaxLevels (solver, max_levels)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetCoarsenType (solver, coarsen_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: coarsen_type
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetCoarsenType (solver, coarsen_type)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetInterpType (solver, interp_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: interp_type
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetInterpType (solver, interp_type)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetNumSweeps (solver, num_sweeps, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: num_sweeps
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetNumSweeps (solver, num_sweeps)
  end subroutine

  subroutine fHYPRE_ParCSRHybridSetRelaxType (solver, relax_type, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(in) :: relax_type
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridSetRelaxType (solver, relax_type)
  end subroutine

  subroutine fHYPRE_ParCSRHybridGetNumIterations (solver, num_its, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(out) :: num_its
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridGetNumIterations (solver, num_its)
  end subroutine

  subroutine fHYPRE_ParCSRHybridGetDSCGNumIterations (solver, dscg_num_its, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(out) :: dscg_num_its
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridGetDSCGNumIterations (solver, dscg_num_its)
  end subroutine

  subroutine fHYPRE_ParCSRHybridGetPCGNumIterations (solver, pcg_num_its, ierr)
    type(c_ptr), intent(in) :: solver
    integer, intent(out) :: pcg_num_its
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridGetPCGNumIterations (solver, pcg_num_its)
  end subroutine

  subroutine fHYPRE_ParCSRHybridGetFinalRelativeResidualNorm (solver, norm, ierr)
    type(c_ptr), intent(in) :: solver
    real(r8), intent(out) :: norm
    integer, intent(out) :: ierr
    ierr = HYPRE_ParCSRHybridGetFinalRelativeResidualNorm (solver, norm)
  end subroutine

  !!!! MISCELLANEOUS PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fHYPRE_ClearAllErrors ()
    integer :: ierr
    ierr = HYPRE_ClearAllErrors ()  ! ignore return code
  end subroutine

end module fhypre
