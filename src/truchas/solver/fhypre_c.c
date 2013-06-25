/*
 *  fhypre_c.c -- C glue code between the Hypre package's C API and
 *  the Fortran API implemented in fhypre.F90.  The function return
 *  values converted to final void function arguments.  MPI communicator
 *  hardwired to MPI_COMM_WORLD and not exposed.  Structure pointer
 *  arguments converted to void pointers.  Arguments passed by value
 *  converted to pass-by-reference.  Function names prepended with
 *  an "f" and shortened in a few cases.  Fortran name mangling assumes
 *  lower case, single appended underscore, but other manglings can be
 *  implemented with the addition or modification of a couple macros.
 *  Some sequences of function calls have been condensed into single
 *  calls; see the IJMatrixCreate, IJVectorCreate, BoomerAMGSetup/Solve.
 */

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_ls.h"

#include <FortranCInterface_names.h>

/* Fortran-callable function names */
#define FHYPRE_IJMATRIXCREATE            TR_ROUTINE_GLOBAL_(fhypre_ijmatrixcreate,FHYPRE_IJMATRIXCREATE)
#define FHYPRE_IJMATRIXDESTROY           TR_ROUTINE_GLOBAL_(fhypre_ijmatrixdestroy, FHYPRE_IJMATRIXDESTROY)
#define FHYPRE_IJMATRIXSETROWSIZES       TR_ROUTINE_GLOBAL_(fhypre_ijmatrixsetrowsizes,FHYPRE_IJMATRIXSETROWSIZES)
#define FHYPRE_IJMATRIXSETDIAGOFFDSIZES  TR_ROUTINE_GLOBAL_(fhypre_ijmatrixsetdiagoffdsizes,FHYPRE_IJMATRIXSETDIAGOFFDSIZES)
#define FHYPRE_IJMATRIXSETMAXOFFPVALUES  TR_ROUTINE_GLOBAL_(fhypre_ijmatrixsetmaxoffpvalues,FHYPRE_IJMATRIXSETMAXOFFPVALUES)
#define FHYPRE_IJMATRIXINITIALIZE        TR_ROUTINE_GLOBAL_(fhypre_ijmatrixinitialize,FHYPRE_IJMATRIXINITIALIZE)
#define FHYPRE_IJMATRIXSETVALUES         TR_ROUTINE_GLOBAL_(fhypre_ijmatrixsetvalues,FHYPRE_IJMATRIXSETVALUES)
#define FHYPRE_IJMATRIXASSEMBLE          TR_ROUTINE_GLOBAL_(fhypre_ijmatrixassemble,FHYPRE_IJMATRIXASSEMBLE)
  
#define FHYPRE_IJVECTORCREATE            TR_ROUTINE_GLOBAL_(fhypre_ijvectorcreate,FHYPRE_IJVECTORCREATE)
#define FHYPRE_IJVECTORDESTROY           TR_ROUTINE_GLOBAL_(fhypre_ijvectordestroy,FHYPRE_IJVECTORDESTROY)
#define FHYPRE_IJVECTORSETMAXOFFPVALUES  TR_ROUTINE_GLOBAL_(fhypre_ijvectorsetmaxoffpvalues,FHYPRE_IJVECTORSETMAXOFFPVALUES)
#define FHYPRE_IJVECTORINITIALIZE        TR_ROUTINE_GLOBAL_(fhypre_ijvectorinitialize,FHYPRE_IJVECTORINITIALIZE)
#define FHYPRE_IJVECTORSETVALUES         TR_ROUTINE_GLOBAL_(fhypre_ijvectorsetvalues,FHYPRE_IJVECTORSETVALUES)
#define FHYPRE_IJVECTORGETVALUES         TR_ROUTINE_GLOBAL_(fhypre_ijvectorgetvalues,FHYPRE_IJVECTORGETVALUES)
#define FHYPRE_IJVECTORASSEMBLE          TR_ROUTINE_GLOBAL_(fhypre_ijvectorassemble,FHYPRE_IJVECTORASSEMBLE)

#define FHYPRE_BOOMERAMGCREATE           TR_ROUTINE_GLOBAL_(fhypre_boomeramgcreate,FHYPRE_BOOMERAMGCREATE)
#define FHYPRE_BOOMERAMGDESTROY          TR_ROUTINE_GLOBAL_(fhypre_boomeramgdestroy,FHYPRE_BOOMERAMGDESTROY)
#define FHYPRE_BOOMERAMGSETUP            TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetup,FHYPRE_BOOMERAMGSETUP)
#define FHYPRE_BOOMERAMGSOLVE            TR_ROUTINE_GLOBAL_(fhypre_boomeramgsolve,FHYPRE_BOOMERAMGSOLVE)
#define FHYPRE_BOOMERAMGSETMAXITER       TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetmaxiter,FHYPRE_BOOMERAMGSETMAXITER)
#define FHYPRE_BOOMERAMGSETPRINTLEVEL    TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetprintlevel,FHYPRE_BOOMERAMGSETPRINTLEVEL)
#define FHYPRE_BOOMERAMGSETRELAXTYPE     TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetrelaxtype,FHYPRE_BOOMERAMGSETRELAXTYPE)
#define FHYPRE_BOOMERAMGSETCOARSENTYPE   TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetcoarsentype,FHYPRE_BOOMERAMGSETCOARSENTYPE)
#define FHYPRE_BOOMERAMGSETNUMSWEEPS     TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetnumsweeps,FHYPRE_BOOMERAMGSETNUMSWEEPS)
#define FHYPRE_BOOMERAMGSETMAXLEVELS     TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetmaxlevels,FHYPRE_BOOMERAMGSETMAXLEVELS)
#define FHYPRE_BOOMERAMGSETTOL           TR_ROUTINE_GLOBAL_(fhypre_boomeramgsettol,FHYPRE_BOOMERAMGSETTOL)
#define FHYPRE_BOOMERAMGSETSTRONGTHLD    TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetstrongthld,FHYPRE_BOOMERAMGSETSTRONGTHLD)
#define FHYPRE_BOOMERAMGSETCYCLENUMSWEEPS TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetcyclenumsweeps,FHYPRE_BOOMERAMGSETCYCLENUMSWEEPS)
#define FHYPRE_BOOMERAMGSETCYCLERELAXTYPE TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetcyclerelaxtype,FHYPRE_BOOMERAMGSETCYCLERELAXTYPE)
#define FHYPRE_BOOMERAMGSETCYCLETYPE     TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetcycletype,FHYPRE_BOOMERAMGSETCYCLETYPE)
#define FHYPRE_BOOMERAMGSETDEBUGFLAG     TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetdebugflag,FHYPRE_BOOMERAMGSETDEBUGFLAG)
#define FHYPRE_BOOMERAMGSETLOGGING       TR_ROUTINE_GLOBAL_(fhypre_boomeramgsetlogging,FHYPRE_BOOMERAMGSETLOGGING)

#define FHYPRE_PCGCREATE                 TR_ROUTINE_GLOBAL_(fhypre_pcgcreate,FHYPRE_PCGCREATE)
#define FHYPRE_PCGDESTROY                TR_ROUTINE_GLOBAL_(fhypre_pcgdestroy,FHYPRE_PCGDESTROY)
#define FHYPRE_PCGSETUP                  TR_ROUTINE_GLOBAL_(fhypre_pcgsetup,FHYPRE_PCGSETUP)
#define FHYPRE_PCGSOLVE                  TR_ROUTINE_GLOBAL_(fhypre_pcgsolve,FHYPRE_PCGSOLVE)
#define FHYPRE_PCGSETMAXITER             TR_ROUTINE_GLOBAL_(fhypre_pcgsetmaxiter,FHYPRE_PCGSETMAXITER)
#define FHYPRE_PCGSETPRINTLEVEL          TR_ROUTINE_GLOBAL_(fhypre_pcgsetprintlevel,FHYPRE_PCGSETPRINTLEVEL)
#define FHYPRE_PCGSETTOL                 TR_ROUTINE_GLOBAL_(fhypre_pcgsettol,FHYPRE_PCGSETTOL)
#define FHYPRE_PCGSETABSTOL              TR_ROUTINE_GLOBAL_(fhypre_pcgsetabstol,FHYPRE_PCGSETABSTOL)
#define FHYPRE_PCGSETTWONORM             TR_ROUTINE_GLOBAL_(fhypre_pcgsettwonorm,FHYPRE_PCGSETTWONORM)
#define FHYPRE_PCGGETFINALRELRES         TR_ROUTINE_GLOBAL_(fhypre_pcggetfinalrelres,FHYPRE_PCGGETFINALRELRES)
#define FHYPRE_PCGSETPRECOND             TR_ROUTINE_GLOBAL_(fhypre_pcgsetprecond,FHYPRE_PCGSETPRECOND)
#define FHYPRE_PCGGETNUMITERATIONS       TR_ROUTINE_GLOBAL_(fhypre_pcggetnumiterations,FHYPRE_PCGGETNUMITERATIONS)

#define FHYPRE_CLEARALLERRORS            TR_ROUTINE_GLOBAL_(fhypre_clearallerrors,FHYPRE_CLEARALLERRORS)

/*** IJMatrix system interface ************************************************/

void
FHYPRE_IJMATRIXCREATE(int *ilower, int *iupper, int *jlower, int *jupper,
                      void **matrix, int *ierr)
{
  *ierr = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, *ilower, *iupper, *jlower, *jupper,
                               (HYPRE_IJMatrix*) matrix);
  *ierr = HYPRE_IJMatrixSetObjectType((HYPRE_IJMatrix) *matrix, HYPRE_PARCSR);
}

void
FHYPRE_IJMATRIXDESTROY(void **matrix, int *ierr)
{
  *ierr = HYPRE_IJMatrixDestroy((HYPRE_IJMatrix) *matrix);
}

void
FHYPRE_IJMATRIXSETROWSIZES(void **matrix, int *sizes, int *ierr)
{
  *ierr = HYPRE_IJMatrixSetRowSizes((HYPRE_IJMatrix) *matrix, sizes);
}

void
FHYPRE_IJMATRIXSETDIAGOFFDSIZES(void **matrix, int *diag_sizes, int *offdiag_sizes, int *ierr)
{
  *ierr = HYPRE_IJMatrixSetDiagOffdSizes((HYPRE_IJMatrix) *matrix, diag_sizes, offdiag_sizes);
}

void
FHYPRE_IJMATRIXSETMAXOFFPVALUES(void **matrix, int *max_offp_values, int *ierr)
{
  *ierr = HYPRE_IJMatrixSetMaxOffProcElmts((HYPRE_IJMatrix) *matrix, *max_offp_values);
}

void
FHYPRE_IJMATRIXINITIALIZE(void **matrix, int *ierr)
{
  *ierr = HYPRE_IJMatrixInitialize((HYPRE_IJMatrix) *matrix);
}

void
FHYPRE_IJMATRIXSETVALUES(void **matrix, int *nrows, int *ncols, int *rows, int *cols, double *values, int *ierr)
{
  *ierr = HYPRE_IJMatrixSetValues((HYPRE_IJMatrix) *matrix, *nrows, ncols, rows, cols, values);
}

void
FHYPRE_IJMATRIXASSEMBLE(void **matrix, int *ierr)
{
  *ierr = HYPRE_IJMatrixAssemble((HYPRE_IJMatrix) *matrix);
}

/*** IJVector system interface ************************************************/

void
FHYPRE_IJVECTORCREATE(int *jlower, int *jupper, void **vector, int *ierr)
{
  *ierr = HYPRE_IJVectorCreate(MPI_COMM_WORLD, *jlower, *jupper, (HYPRE_IJVector*) vector);
  *ierr = HYPRE_IJVectorSetObjectType((HYPRE_IJVector) *vector, HYPRE_PARCSR);
}

void
FHYPRE_IJVECTORDESTROY(void **vector, int *ierr)
{
  *ierr = HYPRE_IJVectorDestroy((HYPRE_IJVector) *vector);
}

void
FHYPRE_IJVECTORSETMAXOFFPVALUES(void **vector, int *max_offp_values, int *ierr)
{
  *ierr = HYPRE_IJVectorSetMaxOffProcElmts((HYPRE_IJVector) *vector, *max_offp_values);
}

void
FHYPRE_IJVECTORINITIALIZE(void **vector, int *ierr)
{
  *ierr = HYPRE_IJVectorInitialize((HYPRE_IJVector) *vector);
}

void
FHYPRE_IJVECTORSETVALUES(void **vector, int *nvalues, int *indices, double *values, int *ierr)
{
  *ierr = HYPRE_IJVectorSetValues((HYPRE_IJVector) *vector, *nvalues, indices, values);
}

void
FHYPRE_IJVECTORGETVALUES(void **vector, int *nvalues, int *indices, double *values, int *ierr)
{
  *ierr = HYPRE_IJVectorGetValues((HYPRE_IJVector) *vector, *nvalues, indices, values);
}

void
FHYPRE_IJVECTORASSEMBLE(void **vector, int *ierr)
{
  *ierr = HYPRE_IJVectorAssemble((HYPRE_IJVector) *vector);
}

/*** ParCSR BoomerAMG interface ***********************************************/

void
FHYPRE_BOOMERAMGCREATE(void **solver, int *ierr)
{
  *ierr = HYPRE_BoomerAMGCreate((HYPRE_Solver*) solver);
}

void
FHYPRE_BOOMERAMGDESTROY(void **solver, int *ierr)
{
  *ierr = HYPRE_BoomerAMGDestroy((HYPRE_Solver) *solver);
}

void
FHYPRE_BOOMERAMGSETUP(void **solver, void **A, void **b, void **x, int *ierr)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix) *A, (void**) &par_A);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *b, (void**) &par_b);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *x, (void**) &par_x);
  *ierr = HYPRE_BoomerAMGSetup((HYPRE_Solver) *solver, par_A, par_b, par_x);
}

void
FHYPRE_BOOMERAMGSOLVE(void **solver, void **A, void **b, void **x, int *ierr)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix) *A, (void**) &par_A);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *b, (void**) &par_b);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *x, (void**) &par_x);
  *ierr = HYPRE_BoomerAMGSolve((HYPRE_Solver) *solver, par_A, par_b, par_x);
}

void
FHYPRE_BOOMERAMGSETSTRONGTHLD(void **solver, double *strong_threshold, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetStrongThreshold((HYPRE_Solver) *solver, *strong_threshold);
}

void
FHYPRE_BOOMERAMGSETMAXITER(void **solver, int *max_iter, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetMaxIter((HYPRE_Solver) *solver, *max_iter);
}

void
FHYPRE_BOOMERAMGSETCOARSENTYPE(void **solver, int *coarsen_type, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetCoarsenType((HYPRE_Solver) *solver, *coarsen_type);
}

void
FHYPRE_BOOMERAMGSETTOL(void **solver, double *tol, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetTol((HYPRE_Solver) *solver, *tol);
}

void
FHYPRE_BOOMERAMGSETNUMSWEEPS(void **solver, int *num_sweeps, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver) *solver, *num_sweeps);
}

void
FHYPRE_BOOMERAMGSETRELAXTYPE(void **solver, int *relax_type, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetRelaxType((HYPRE_Solver) *solver, *relax_type);
}

void
FHYPRE_BOOMERAMGSETPRINTLEVEL(void **solver, int *print_level, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetPrintLevel((HYPRE_Solver) *solver, *print_level);
}

void
FHYPRE_BOOMERAMGSETMAXLEVELS(void **solver, int *max_levels, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetMaxLevels((HYPRE_Solver) *solver, *max_levels);
}

void
FHYPRE_BOOMERAMGSETCYCLENUMSWEEPS(void **solver, int *num_sweeps, int *k, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetCycleNumSweeps((HYPRE_Solver) *solver, *num_sweeps, *k);
}

void
FHYPRE_BOOMERAMGSETCYCLERELAXTYPE(void **solver, int *relax_type, int *k, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetCycleRelaxType((HYPRE_Solver) *solver, *relax_type, *k);
}

void
FHYPRE_BOOMERAMGSETCYCLETYPE(void **solver, int *cycle_type, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetCycleType((HYPRE_Solver) *solver, *cycle_type);
}

void
FHYPRE_BOOMERAMGSETDEBUGFLAG(void **solver, int *debug, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetDebugFlag((HYPRE_Solver) *solver, *debug);
}

void
FHYPRE_BOOMERAMGSETLOGGING(void **solver, int *logging, int *ierr)
{
  *ierr = HYPRE_BoomerAMGSetDebugFlag((HYPRE_Solver) *solver, *logging);
}

/*** ParCSR PCG interface ***********************************************/

void
FHYPRE_PCGCREATE(void **solver, int *ierr)
{
  *ierr = HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, (HYPRE_Solver*) solver);
}

void
FHYPRE_PCGDESTROY(void **solver, int *ierr)
{
  *ierr = HYPRE_ParCSRPCGDestroy((HYPRE_Solver) *solver);
}

void
FHYPRE_PCGSETUP(void **solver, void **A, void **b, void **x, int *ierr)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix) *A, (void**) &par_A);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *b, (void**) &par_b);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *x, (void**) &par_x);
  *ierr = HYPRE_ParCSRPCGSetup((HYPRE_Solver) *solver, par_A, par_b, par_x);
}

void
FHYPRE_PCGSOLVE(void **solver, void **A, void **b, void **x, int *ierr)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix) *A, (void**) &par_A);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *b, (void**) &par_b);
  HYPRE_IJVectorGetObject((HYPRE_IJVector) *x, (void**) &par_x);
  *ierr = HYPRE_ParCSRPCGSolve((HYPRE_Solver) *solver, par_A, par_b, par_x);
}

void
FHYPRE_PCGSETMAXITER(void **solver, int *max_iter, int *ierr)
{
  *ierr = HYPRE_PCGSetMaxIter((HYPRE_Solver) *solver, *max_iter);
}

void
FHYPRE_PCGSETTOL(void **solver, double *tol, int *ierr)
{
  *ierr = HYPRE_PCGSetTol((HYPRE_Solver) *solver, *tol);
}

void
FHYPRE_PCGSETABSTOL(void **solver, double *tol, int *ierr)
{
  *ierr = HYPRE_PCGSetAbsoluteTol((HYPRE_Solver) *solver, *tol);
}

void
FHYPRE_PCGSETTWONORM(void **solver, int *twonorm, int *ierr)
{
  *ierr = HYPRE_PCGSetTwoNorm((HYPRE_Solver) *solver, *twonorm);
}

void
FHYPRE_PCGGETFINALRELRES(void **solver, double *tol, int *ierr)
{
  *ierr = HYPRE_PCGGetFinalRelativeResidualNorm((HYPRE_Solver) *solver, tol);
}

void
FHYPRE_PCGSETPRINTLEVEL(void **solver, int *print_level, int *ierr)
{
  *ierr = HYPRE_PCGSetPrintLevel((HYPRE_Solver) *solver, *print_level);
}

void
FHYPRE_PCGSETPRECOND(void **solver, void **precond, int *ierr)
{
  *ierr = HYPRE_ParCSRPCGSetPrecond( (HYPRE_Solver) *solver, 
                                HYPRE_BoomerAMGSolve,
                                HYPRE_BoomerAMGSetup,
                               (HYPRE_Solver) *precond ); 
}

void
FHYPRE_PCGGETNUMITERATIONS(void **solver, int *iters, int *ierr)
{
  *ierr = HYPRE_PCGGetNumIterations( (HYPRE_Solver) *solver, iters);
}

/*** HYPRE_Utilities ************************************************/

FHYPRE_CLEARALLERRORS()
{
  int ierr = HYPRE_ClearAllErrors();
}
