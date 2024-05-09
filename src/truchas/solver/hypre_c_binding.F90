!!
!! HYPRE_C_BINDING
!!
!! Raw bindings to a subset of the Hypre library C interface.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!!  This uses the C interoperability features of F2003 to interface directly
!!  to Hypre's C API functions, and it replaces the fragile F77-callable C glue
!!  code from fhypre_c.c
!!
!!  Bindings to just a small part of the Hypre library are defined; only the
!!  things used by Truchas.
!!
!!  The various typedefs used in the C function signatures, like HYPRE_Solver
!!  and HYPRE_IJMatrix, are structure pointers and on the Fortran side should
!!  be considered as c_ptr type arguments.  These are simply opaque handles
!!  that get passed between Hypre functions.  Some C arguments are passed by
!!  value and some by reference, but those details can be ignored from the
!!  Fortran side because the specification of the interfaces ensures that the
!!  right thing will happen.
!!
!!  This module also defines bindings to a small number of extended Hypre
!!  C functions from hypre_ext.c  Those functions are not part of Hypre, but
!!  provide some higher-level functionality that was more easily implemented
!!  in C than Fortran.
!!

module hypre_c_binding

  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
  implicit none
  private
  
  integer(c_int), parameter, public :: HYPRE_ERROR_GENERIC = 1
  integer(c_int), parameter, public :: HYPRE_ERROR_MEMORY  = 2
  integer(c_int), parameter, public :: HYPRE_ERROR_ARG     = 4
  integer(c_int), parameter, public :: HYPRE_ERROR_CONV    = 256

  !! Functions from the IJVector interface.
  public :: HYPRE_IJVectorDestroy
  public :: HYPRE_IJVectorSetMaxOffProcElmts
  public :: HYPRE_IJVectorInitialize
  public :: HYPRE_IJVectorSetValues
  public :: HYPRE_IJVectorAddToValues
  public :: HYPRE_IJVectorGetValues
  public :: HYPRE_IJVectorAssemble

  !! Functions from the IJMatrix interface.
  public :: HYPRE_IJMatrixDestroy
  public :: HYPRE_IJMatrixSetRowSizes
  public :: HYPRE_IJMatrixSetDiagOffdSizes
  public :: HYPRE_IJMatrixSetMaxOffProcElmts
  public :: HYPRE_IJMatrixInitialize
  public :: HYPRE_IJMatrixAssemble
  public :: HYPRE_IJMatrixSetValues
  public :: HYPRE_IJMatrixAddToValues

  !! Functions from the BoomerAMG interface.
  public :: HYPRE_BoomerAMGCreate
  public :: HYPRE_BoomerAMGDestroy
  public :: HYPRE_BoomerAMGSetStrongThreshold
  public :: HYPRE_BoomerAMGSetMaxIter
  public :: HYPRE_BoomerAMGSetCoarsenType
  public :: HYPRE_BoomerAMGSetInterpType
  public :: HYPRE_BoomerAMGSetTol
  public :: HYPRE_BoomerAMGSetNumSweeps
  public :: HYPRE_BoomerAMGSetRelaxType
  public :: HYPRE_BoomerAMGSetPrintLevel
  public :: HYPRE_BoomerAMGSetMaxLevels
  public :: HYPRE_BoomerAMGSetCycleNumSweeps
  public :: HYPRE_BoomerAMGSetCycleRelaxType
  public :: HYPRE_BoomerAMGSetCycleType
  public :: HYPRE_BoomerAMGSetDebugFlag
  public :: HYPRE_BoomerAMGSetLogging
  public :: HYPRE_BoomerAMGSetOldDefault

  !! Functions from the PCG interface
  public :: HYPRE_PCGDestroy
  public :: HYPRE_PCGSetMaxIter
  public :: HYPRE_PCGGetNumIterations
  public :: HYPRE_PCGSetTol
  public :: HYPRE_PCGSetAbsoluteTol
  public :: HYPRE_PCGSetTwoNorm
  public :: HYPRE_PCGGetFinalRelativeResidualNorm
  public :: HYPRE_PCGSetPrintLevel

  !! Functions from the ParCSR Hybrid interface
  public :: HYPRE_ParCSRHybridCreate
  public :: HYPRE_ParCSRHybridDestroy
  public :: HYPRE_ParCSRHybridSetTol
  public :: HYPRE_ParCSRHybridSetAbsoluteTol
  public :: HYPRE_ParCSRHybridSetConvergenceTol
  public :: HYPRE_ParCSRHybridSetDSCGMaxIter
  public :: HYPRE_ParCSRHybridSetPCGMaxIter
  public :: HYPRE_ParCSRHybridSetSolverType
  public :: HYPRE_ParCSRHybridSetKDim
  public :: HYPRE_ParCSRHybridSetTwoNorm
  public :: HYPRE_ParCSRHybridSetLogging
  public :: HYPRE_ParCSRHybridSetPrintLevel
  public :: HYPRE_ParCSRHybridSetStrongThreshold
  public :: HYPRE_ParCSRHybridSetMaxLevels
  public :: HYPRE_ParCSRHybridSetCoarsenType
  public :: HYPRE_ParCSRHybridSetInterpType
  public :: HYPRE_ParCSRHybridSetNumSweeps
  public :: HYPRE_ParCSRHybridSetRelaxType
  public :: HYPRE_ParCSRHybridGetNumIterations
  public :: HYPRE_ParCSRHybridGetDSCGNumIterations
  public :: HYPRE_ParCSRHybridGetPCGNumIterations
  public :: HYPRE_ParCSRHybridGetFinalRelativeResidualNorm

  !! Miscelaneous functions
  public :: HYPRE_ClearAllErrors
  public :: HYPRE_Initialize
  public :: HYPRE_Finalize

  !! Functions from hypre_ext.c
  public :: HYPRE_Ext_IJVectorCreate
  public :: HYPRE_Ext_IJMatrixCreate
  public :: HYPRE_Ext_BoomerAMGSetup
  public :: HYPRE_Ext_BoomerAMGSolve
  public :: HYPRE_Ext_ParCSRPCGCreate
  public :: HYPRE_Ext_PCGSetup
  public :: HYPRE_Ext_PCGSolve
  public :: HYPRE_Ext_PCGSetBoomerAMGPrecond
  public :: HYPRE_Ext_HybridSetup
  public :: HYPRE_Ext_HybridSolve

 !!
 !! IJVECTOR INTERFACES
 !!

  interface
    ! See HYPRE_Ext_IJVectorCreate below; wraps HYPRE_IJVectorCreate.
    function HYPRE_IJVectorDestroy(vector) &
        result(ierr) bind(c, name="HYPRE_IJVectorDestroy")
      import c_ptr, c_int
      type(c_ptr), value :: vector
      integer(c_int) :: ierr
    end function
    function HYPRE_IJVectorSetMaxOffProcElmts(vector, max_offp_values) &
        result(ierr) bind(c, name="HYPRE_IJVectorSetMaxOffProcElmts")
      import c_ptr, c_int
      type(c_ptr), value :: vector
      integer(c_int), value :: max_offp_values
      integer(c_int) :: ierr
    end function
    function HYPRE_IJVectorInitialize(vector) &
        result(ierr) bind(c, name="HYPRE_IJVectorInitialize")
      import c_ptr, c_int
      type(c_ptr), value :: vector
      integer(c_int) :: ierr
    end function
    function HYPRE_IJVectorSetValues(vector, nvalues, indices, values) &
        result(ierr) bind(c, name="HYPRE_IJVectorSetValues")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: vector
      integer, value :: nvalues
      integer(c_int), intent(in) :: indices(*)
      real(c_double), intent(in) :: values(*)
      integer(c_int) :: ierr
    end function
    function HYPRE_IJVectorAddToValues(vector, nvalues, indices, values) &
        result(ierr) bind(c, name="HYPRE_IJVectorAddToValues")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: vector
      integer, value :: nvalues
      integer(c_int), intent(in) :: indices(*)
      real(c_double), intent(in) :: values(*)
      integer(c_int) :: ierr
    end function
    function HYPRE_IJVectorGetValues(vector, nvalues, indices, values) &
        result(ierr) bind(c, name="HYPRE_IJVectorGetValues")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: vector
      integer(c_int), value :: nvalues
      integer(c_int), intent(in) :: indices(*)
      real(c_double), intent(out) :: values(*)
      integer(c_int) :: ierr
    end function
    function HYPRE_IJVectorAssemble(vector) &
        result(ierr) bind(c, name="HYPRE_IJVectorAssemble")
      import c_ptr, c_int
      type(c_ptr), value :: vector
      integer(c_int) :: ierr
    end function
  end interface

 !!
 !! IJMATRIX INTERFACES
 !!

  interface
    ! See HYPRE_Ext_IJMatrixCreate below; wraps HYPRE_IJMatrixCreate.
    function HYPRE_IJMatrixDestroy(matrix) &
        result(ierr) bind(c, name="HYPRE_IJMatrixDestroy")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixSetRowSizes(matrix, sizes) &
        result(ierr) bind(c, name="HYPRE_IJMatrixSetRowSizes")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int), intent(in) :: sizes(*)
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixSetDiagOffdSizes(matrix, diag_sizes, offd_sizes) &
        result(ierr) bind(c, name="HYPRE_IJMatrixSetDiagOffdSizes")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int), intent(in) :: diag_sizes(*), offd_sizes(*)
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixSetMaxOffProcElmts(matrix, max_offp_values) &
        result(ierr) bind(c, name="HYPRE_IJMatrixSetMaxOffProcElmts")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int), value :: max_offp_values
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixInitialize(matrix) &
        result(ierr) bind(c, name="HYPRE_IJMatrixInitialize")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixAssemble(matrix) &
        result(ierr) bind(c, name="HYPRE_IJMatrixAssemble")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixSetValues(matrix, nrows, ncols, rows, cols, values) &
        result(ierr) bind(c, name="HYPRE_IJMatrixSetValues")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: matrix
      integer(c_int), value :: nrows
      integer(c_int), intent(in) :: ncols(*), rows(*), cols(*)
      real(c_double) :: values(*)
      integer(c_int) :: ierr
    end function
    function HYPRE_IJMatrixAddToValues(matrix, nrows, ncols, rows, cols, values) &
        result(ierr) bind(c, name="HYPRE_IJMatrixAddToValues")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: matrix
      integer(c_int), value :: nrows
      integer(c_int), intent(in) :: ncols(*), rows(*), cols(*)
      real(c_double) :: values(*)
      integer(c_int) :: ierr
    end function
  end interface

 !!
 !! BOOMERAMG INTERFACES
 !!

  interface
    function HYPRE_BoomerAMGCreate(solver) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGCreate")
      import c_ptr, c_int
      type(c_ptr) :: solver
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGDestroy(solver) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGDestroy")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) :: ierr
    end function
    ! See HYPRE_Ext_BoomerAMGSetup below; wraps HYPRE_BoomerAMGSetup.
    ! See HYPRE_Ext_BoomerAMGSolve below; wraps HYPRE_BoomerAMGSolve.
    function HYPRE_BoomerAMGSetStrongThreshold (solver, strong_threshold) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetStrongThreshold")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: strong_threshold
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetMaxIter (solver, max_iter) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetMaxIter")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: max_iter
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetCoarsenType (solver, coarsen_type) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetCoarsenType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: coarsen_type
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetInterpType (solver, interp_type) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetInterpType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: interp_type
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetTol (solver, tol) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetTol")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: tol
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetNumSweeps (solver, num_sweeps) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetNumSweeps")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) , value :: num_sweeps
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetRelaxType (solver, relax_type) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetRelaxType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: relax_type
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetPrintLevel (solver, print_level) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetPrintLevel")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: print_level
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetMaxLevels (solver, max_levels) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetMaxLevels")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: max_levels
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetCycleNumSweeps (solver, num_sweeps, k) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetCycleNumSweeps")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: num_sweeps, k
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetCycleRelaxType (solver, relax_type, k) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetCycleRelaxType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: relax_type, k
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetCycleType (solver, cycle_type) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetCycleType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: cycle_type
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetDebugFlag (solver, debug) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetDebugFlag")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: debug
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetLogging (solver, logging) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetLogging")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: logging
      integer(c_int) :: ierr
    end function
    function HYPRE_BoomerAMGSetOldDefault(solver) &
        result(ierr) bind(c, name="HYPRE_BoomerAMGSetOldDefault")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) :: ierr
    end function
  end interface

 !!
 !! PCG INTERFACES
 !!

  interface
    ! See HYPRE_Ext_ParCSRPCGCreate below; wraps HYPRE_ParCSRPCGCreate.
    function HYPRE_PCGDestroy (solver) &
        result(ierr) bind(c, name="HYPRE_ParCSRPCGDestroy")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) :: ierr
    end function
    ! See HYPRE_Ext_PCGSetup below.
    ! See HYPRE_Ext_PCGSolve below.
    function HYPRE_PCGSetMaxIter(solver, max_iter) &
        result(ierr) bind(c, name="HYPRE_PCGSetMaxIter")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: max_iter
      integer(c_int) :: ierr
    end function
    function HYPRE_PCGGetNumIterations(solver, iters) &
        result(ierr) bind(c, name="HYPRE_PCGGetNumIterations")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), intent(out) :: iters
      integer(c_int) :: ierr
    end function HYPRE_PCGGetNumIterations
    function HYPRE_PCGSetTol(solver, tol) &
        result(ierr) bind(c, name="HYPRE_PCGSetTol")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: tol
      integer(c_int) :: ierr
    end function
    function HYPRE_PCGSetAbsoluteTol(solver, tol) &
        result(ierr) bind(c, name="HYPRE_PCGSetAbsoluteTol")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: tol
      integer(c_int) :: ierr
    end function
    function HYPRE_PCGSetTwoNorm(solver, twonorm) &
        result(ierr) bind(c, name="HYPRE_PCGSetTwoNorm")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: twonorm
      integer(c_int) :: ierr
    end function
    function HYPRE_PCGGetFinalRelativeResidualNorm(solver, tol) &
        result(ierr) bind(c, name="HYPRE_PCGGetFinalRelativeResidualNorm")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), intent(out) :: tol
      integer(c_int) :: ierr
    end function
    function HYPRE_PCGSetPrintLevel(solver, print_level) &
        result(ierr) bind(c, name="HYPRE_PCGSetPrintLevel")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: print_level
      integer(c_int) :: ierr
    end function
  end interface

 !!
 !! ParCSR Hybrid SOLVER INTERFACES
 !!

  interface
    function HYPRE_ParCSRHybridCreate (solver) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridCreate")
      import c_ptr, c_int
      type(c_ptr), intent(inout) :: solver
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridDestroy (solver) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridDestroy")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) :: ierr
    end function
    ! HYPRE_ParCSRHybridSetup: use HYPRE_Ext_HybridSetup below.
    ! HYPRE_ParCSRHybridSolve: use HYPRE_Ext_HybridSolve below.
    function HYPRE_ParCSRHybridSetTol (solver, tol) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetTol")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: tol
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetAbsoluteTol (solver, tol) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetAbsoluteTol")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: tol
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetConvergenceTol (solver, cf_tol) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetConvergenceTol")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: cf_tol
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetDSCGMaxIter (solver, dscg_max_its) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetDSCGMaxIter")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: dscg_max_its
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetPCGMaxIter (solver, pcg_max_its) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetPCGMaxIter")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: pcg_max_its
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetSolverType (solver, solver_type) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetSolverType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: solver_type
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetKDim (solver, k_dim) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetKDim")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: k_dim
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetTwoNorm (solver, two_norm) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetTwoNorm")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: two_norm
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetLogging (solver, logging) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetLogging")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: logging
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetPrintLevel (solver, print_level) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetPrintLevel")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: print_level
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetStrongThreshold (solver, strong_threshold) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetStrongThreshold")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), value :: strong_threshold
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetMaxLevels (solver, max_levels) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetMaxLevels")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: max_levels
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetCoarsenType (solver, coarsen_type) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetCoarsenType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: coarsen_type
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetInterpType (solver, interp_type) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetInterpType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: interp_type
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetNumSweeps (solver, num_sweeps) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetNumSweeps")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: num_sweeps
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridSetRelaxType (solver, relax_type) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridSetRelaxType")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), value :: relax_type
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridGetNumIterations (solver, num_its) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridGetNumIterations")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), intent(out) :: num_its
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridGetDSCGNumIterations (solver, dscg_num_its) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridGetDSCGNumIterations")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), intent(out) :: dscg_num_its
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridGetPCGNumIterations (solver, pcg_num_its) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridGetPCGNumIterations")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int), intent(out) :: pcg_num_its
      integer(c_int) :: ierr
    end function
    function HYPRE_ParCSRHybridGetFinalRelativeResidualNorm (solver, norm) &
        result(ierr) bind(c, name="HYPRE_ParCSRHybridGetFinalRelativeResidualNorm")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: solver
      real(c_double), intent(out) :: norm
      integer(c_int) :: ierr
    end function
  end interface

 !!
 !! HYPRE UTILITIES
 !!

  interface
    function HYPRE_ClearAllErrors() &
        result(ierr) bind(c, name="HYPRE_ClearAllErrors")
      import c_int
      integer(c_int) :: ierr
    end function
    function HYPRE_Initialize() &
        result(ierr) bind(c, name="HYPRE_Initialize")
      import c_int
      integer(c_int) :: ierr
    end function
    function HYPRE_Finalize() &
        result(ierr) bind(c, name="HYPRE_Finalize")
      import c_int
      integer(c_int) :: ierr
    end function
  end interface

  !!
  !! HYPRE_EXT -- our own higher-level extensions to Hypre.
  !!

  interface
    function HYPRE_Ext_IJVectorCreate(jlower, jupper, vector) &
        result(ierr) bind(c, name="HYPRE_Ext_IJVectorCreate")
      import c_ptr, c_int
      integer(c_int), value :: jlower, jupper
      type(c_ptr) :: vector
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_IJMatrixCreate(ilower, iupper, jlower, jupper, matrix) &
        result(ierr) bind(c, name="HYPRE_Ext_IJMatrixCreate")
      import c_ptr, c_int
      integer(c_int), value :: ilower, iupper, jlower, jupper
      type(c_ptr) :: matrix
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_BoomerAMGSetup(solver, A, b, x) &
        result(ierr) bind(c, name="HYPRE_Ext_BoomerAMGSetup")
      import c_ptr, c_int
      type(c_ptr), value :: solver, A, b, x
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_BoomerAMGSolve(solver, A, b, x) &
        result(ierr) bind(c, name="HYPRE_Ext_BoomerAMGSolve")
      import c_ptr, c_int
      type(c_ptr), value :: solver, A, b, x
      integer(c_int) :: ierr
    end function
    !! PCG Interfaces
    function HYPRE_Ext_ParCSRPCGCreate(solver) &
        result(ierr) bind(c, name="HYPRE_Ext_ParCSRPCGCreate")
      import c_ptr, c_int
      type(c_ptr), intent(inout) :: solver
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_PCGSetup(solver, A, b, x) &
        result(ierr) bind(c, name="HYPRE_Ext_PCGSetup")
      import c_ptr, c_int
      type(c_ptr), value :: solver, A, b, x
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_PCGSolve(solver, A, b, x) &
        result(ierr) bind(c, name="HYPRE_Ext_PCGSolve")
      import c_ptr, c_int
      type(c_ptr), value :: solver, A, b, x
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_PCGSetBoomerAMGPrecond(solver, precond) &
        result(ierr) bind(c, name="HYPRE_Ext_PCGSetBoomerAMGPrecond")
      import c_ptr, c_int
      type(c_ptr), value :: solver, precond
      integer(c_int) :: ierr
    end function
    !! ParCSRHybrid Interfaces
    function HYPRE_Ext_HybridSetup(solver, A, b, x) &
        result(ierr) bind(c, name="HYPRE_Ext_HybridSetup")
      import c_ptr, c_int
      type(c_ptr), value :: solver, A, b, x
      integer(c_int) :: ierr
    end function
    function HYPRE_Ext_HybridSolve(solver, A, b, x) &
        result(ierr) bind(c, name="HYPRE_Ext_HybridSolve")
      import c_ptr, c_int
      type(c_ptr), value :: solver, A, b, x
      integer(c_int) :: ierr
    end function
  end interface

end module hypre_c_binding
