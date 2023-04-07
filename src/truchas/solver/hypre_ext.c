/*==============================================================================

  This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.

==============================================================================*/

/*
 *  This is a completely reworked version of fhypre_c.c, which contained
 *  C glue code between the Hypre package's C API and the Fortran API
 *  found in fhypre.F90.  The C interoperability features of Fortran 2003
 *  have rendered most of that unnecessary.  What remains here are a few
 *  high-level extensions that are more easily implemented in C than Fortran.
 *  There are wrappers for some object creation functions that hardwire the
 *  MPI communicator to MPI_COM_WORLD.  The rest wrap boiler plate calls
 *  required by the solvers.
 *
 *  Neil N. Carlson <nnc@lanl.gov> March 2014
 */

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_ls.h"

int
HYPRE_Ext_IJMatrixCreate(int ilower, int iupper, int jlower, int jupper, HYPRE_IJMatrix *matrix)
{
  int ierr;
  ierr = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, jlower, jupper, matrix);
  if (ierr) return ierr;
  ierr = HYPRE_IJMatrixSetObjectType(*matrix, HYPRE_PARCSR);
  return ierr;
}

int
HYPRE_Ext_IJVectorCreate(int jlower, int jupper, HYPRE_IJVector *vector)
{
  int ierr;
  ierr = HYPRE_IJVectorCreate(MPI_COMM_WORLD, jlower, jupper, vector);
  if (ierr) return ierr;
  ierr = HYPRE_IJVectorSetObjectType(*vector, HYPRE_PARCSR);
  return ierr;
}

int
HYPRE_Ext_BoomerAMGSetup(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_BoomerAMGSolve(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_ParCSRPCGCreate(HYPRE_Solver *solver)
{
  return HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, solver);
}

/* NNC: I don't recall the background for the next two functions and how they
 * came to be.  The ParCSRPCGSetup/Solve functions aren't documented in the
 * reference manual.  One would think the documented HYPRE_PCGSetup/Solve
 * functions should replace these, but they don't seem to work (I just tried). */

int
HYPRE_Ext_PCGSetup(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_ParCSRPCGSetup(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_PCGSolve(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_ParCSRPCGSolve(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_PCGSetBoomerAMGPrecond(HYPRE_Solver solver, HYPRE_Solver precond)
{
  return HYPRE_ParCSRPCGSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
}

int
HYPRE_Ext_HybridSetup(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_ParCSRHybridSetup(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_HybridSolve(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_ParCSRHybridSolve(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_AMSSetup(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_AMSSetup(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_AMSSolve(HYPRE_Solver solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix par_A;
  HYPRE_ParVector par_b, par_x;
  HYPRE_IJMatrixGetObject(A, (void**) &par_A);
  HYPRE_IJVectorGetObject(b, (void**) &par_b);
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  return HYPRE_AMSSolve(solver, par_A, par_b, par_x);
}

int
HYPRE_Ext_AMSSetDiscreteGradient(HYPRE_Solver solver, HYPRE_IJMatrix G)
{
  HYPRE_ParCSRMatrix par_G;
  HYPRE_IJMatrixGetObject(G, (void**) &par_G);
  return HYPRE_AMSSetDiscreteGradient(solver, par_G);
}

int
HYPRE_Ext_AMSSetCoordinateVectors(HYPRE_Solver solver,
                                  HYPRE_IJVector x, HYPRE_IJVector y, HYPRE_IJVector z)
{
  HYPRE_ParVector par_x, par_y, par_z;
  HYPRE_IJVectorGetObject(x, (void**) &par_x);
  HYPRE_IJVectorGetObject(y, (void**) &par_y);
  HYPRE_IJVectorGetObject(z, (void**) &par_z);
  return HYPRE_AMSSetCoordinateVectors(solver, par_x, par_y, par_z);
}

int
HYPRE_Ext_AMSSetEdgeConstantVectors(HYPRE_Solver solver,
                                    HYPRE_IJVector Gx, HYPRE_IJVector Gy, HYPRE_IJVector Gz)
{
  HYPRE_ParVector par_Gx, par_Gy, par_Gz;
  HYPRE_IJVectorGetObject(Gx, (void**) &par_Gx);
  HYPRE_IJVectorGetObject(Gy, (void**) &par_Gy);
  HYPRE_IJVectorGetObject(Gz, (void**) &par_Gz);
  return HYPRE_AMSSetEdgeConstantVectors(solver, par_Gx, par_Gy, par_Gz);
}

int
HYPRE_Ext_AMSSetAlphaPoissonMatrix(HYPRE_Solver solver, HYPRE_IJMatrix A_alpha)
{
  HYPRE_ParCSRMatrix par_A_alpha;
  HYPRE_IJMatrixGetObject(A_alpha, (void**) &par_A_alpha);
  return HYPRE_AMSSetAlphaPoissonMatrix(solver, par_A_alpha);
}

int
HYPRE_Ext_AMSSetBetaPoissonMatrix(HYPRE_Solver solver, HYPRE_IJMatrix A_beta)
{
  HYPRE_ParCSRMatrix par_A_beta;
  HYPRE_IJMatrixGetObject(A_beta, (void**) &par_A_beta);
  return HYPRE_AMSSetBetaPoissonMatrix(solver, par_A_beta);
}

int
HYPRE_Ext_AMSSetInteriorNodes(HYPRE_Solver solver, HYPRE_IJVector interior_nodes)
{
  HYPRE_ParVector par_interior_nodes;
  HYPRE_IJVectorGetObject(interior_nodes, (void**) &par_interior_nodes);
  return HYPRE_AMSSetInteriorNodes(solver, par_interior_nodes);
}
