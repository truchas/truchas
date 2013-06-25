/****************************************************************************\
 *                                                                          *
 *   Copyright (c) 1995, 2000, 2005 Sandia Corporation.                           *
 *   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, *
 *   the U.S. Government retains certain rights in this software.           *
 *   For more info, see the README file in the top-level directory.         * 
 *                                                                          *
\****************************************************************************/

/*
@(#)++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@(#)
@(#)    $RCSfile: VF_RadSolve.c,v $
@(#)    $Revision: 1.6 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_RadSolve.c,v $
@(#)
@(#)    DESCRIPTION:  solver interface to solve the radiosity matrix eqn.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double *radq_g=NULL;
static double *radj_g=NULL;
static double *tsurf_g=NULL;
static double *eps_g=NULL;

#include "vf.h"

void 
VF_RadSolve(int encl,
            double *radq, double *tsurf, double *eps, double sigma, 
            double tol, int max_iter, int *num_iter, 
            int method, int dynamic, int debug)
{
  int    i, ii, max_surfaces;
  static int    max_size=0;
  VFenclosure   *enclosure=VF_GetEnclosure(encl);

  max_surfaces = VF_MaxPatches();
  if (max_size==0) {
    max_size = max_surfaces;
    radq_g   = VF_Newd(max_size);
    tsurf_g  = VF_Newd(max_size);
    eps_g    = VF_Newd(max_size);
  } else if (max_surfaces>max_size) {
    max_size = max_surfaces;
    radq_g   = VF_ReNewd(radq_g ,max_size);
    tsurf_g  = VF_ReNewd(tsurf_g,max_size);
    eps_g    = VF_ReNewd(eps_g  ,max_size);
  }
  VF_Gather_RadiosityVectors(tsurf,eps,radq,tsurf_g,eps_g,radq_g);
  switch (method) {
  case VF_RADSOLVE_CG:
    VF_RadSolveCG(radq_g, NULL, tsurf_g, eps_g, sigma, 
                  tol, max_iter, num_iter, debug);
    break;
  case VF_RADSOLVE_GMRES:
    VF_RadSolveGMRES(radq_g, NULL, tsurf_g, eps_g, sigma, 
                     tol, max_iter, num_iter, debug);
    break;
  case VF_RADSOLVE_AZTEC_CG:
    VF_RadSolveAZTEC(radq_g, NULL, tsurf_g, eps_g, sigma, 
                     tol, max_iter, num_iter, dynamic, method, debug);
    break;
  case VF_RADSOLVE_AZTEC_GMRES:
    VF_RadSolveAZTEC(radq_g, NULL, tsurf_g, eps_g, sigma,
                     tol, max_iter, num_iter, dynamic, method, debug);
    break;
  }
  for (i=0; i<enclosure->host_npatches; i++) {
    ii       = enclosure->host2vflib_map[i];
    radq[i]  = radq_g[ii];
    tsurf[i] = tsurf_g[ii];
    eps[i]   = eps_g[ii];
  }
}

void 
VF_RadSolveAux(int encl,
               double *radq, double *radj, double *tsurf, double *eps, 
               double sigma, double tol, int max_iter, int *num_iter, 
               int method, int dynamic, int debug)
{
  int    i, ii, max_surfaces;
  static int    max_size=0;
  VFenclosure   *enclosure=VF_GetEnclosure(encl);

  max_surfaces = VF_MaxPatches();
  if (max_size==0) {
    max_size = max_surfaces;
    radq_g   = VF_Newd(max_size);
    radj_g   = VF_Newd(max_size);
    tsurf_g  = VF_Newd(max_size);
    eps_g    = VF_Newd(max_size);
  } else if (max_surfaces>max_size) {
    max_size = max_surfaces;
    radq_g   = VF_ReNewd(radq_g ,max_size);
    radj_g   = VF_ReNewd(radq_g ,max_size);
    tsurf_g  = VF_ReNewd(tsurf_g,max_size);
    eps_g    = VF_ReNewd(eps_g  ,max_size);
  }
  VF_Gather_RadiosityVectorsAux(tsurf,eps,radq,radj,tsurf_g,eps_g,radq_g,radj_g);
  switch (method) {
  case VF_RADSOLVE_CG:
    VF_RadSolveCG(radq_g, radj_g, tsurf_g, eps_g, sigma, 
                  tol, max_iter, num_iter, debug);
    break;
  case VF_RADSOLVE_GMRES:
    VF_RadSolveGMRES(radq_g, radj_g, tsurf_g, eps_g, sigma, 
                     tol, max_iter, num_iter, debug);
    break;
  case VF_RADSOLVE_AZTEC_CG:
    VF_RadSolveAZTEC(radq_g, radj_g, tsurf_g, eps_g, sigma, 
                     tol, max_iter, num_iter, dynamic, method, debug);
    break;
  case VF_RADSOLVE_AZTEC_GMRES:
    VF_RadSolveAZTEC(radq_g, radj_g, tsurf_g, eps_g, sigma,
                     tol, max_iter, num_iter, dynamic, method, debug);
    break;
  }
  for (i=0; i<enclosure->host_npatches; i++) {
    ii       = enclosure->host2vflib_map[i];
    radq[i]  = radq_g[ii];
    radj[i]  = radj_g[ii];
    tsurf[i] = tsurf_g[ii];
    eps[i]   = eps_g[ii];
  }
}

void VF_RadSolve_Residual(double eps[], double rhs[], double sol[], double r[])
{
  int i;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  VF_RadSolve_MatVecMul(eps, sol, r);
  for (i=0; i<enclosure->npatches_l; i++){
    r[i] = rhs[enclosure->row[i].global_index]-r[i];
  }
}

void VF_RadSolve_MatVecMul(double eps[], double x[], double b[])
{
  int    i, j, k, n, n_l, nonzeros, *index;
  float  *vf;
  double sum, rho;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (n=0; n<n_l; n++, row++) {
    sum      = 0.0;
    i        = row->global_index;
    rho      = 1.0-eps[i];
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      k     = index[j];
      sum += vf[j]*x[k];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      k     = index[j];
      sum += vf[j]*x[k];
    }
    b[n] = x[i]-sum*rho;
  }
}

void VF_RadSolve_Residual_Scale(double eps[], double rhs[], double x[], 
                                double b[], double s[])
{
  int i;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  VF_RadSolve_MatVecMul_Scale(eps, x, b, s);
  for (i=0; i<enclosure->npatches_l; i++){
    b[i] = rhs[enclosure->row[i].global_index]-b[i];
  }
}

void VF_RadSolve_MatVecMul_Scale(double eps[], double x[], 
                                 double b[],   double s[])
{
  int    i, j, k, n, n_l, nonzeros, *index;
  float  *vf;
  double sum, rho;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (n=0; n<n_l; n++, row++) {
    sum      = 0.0;
    i        = row->global_index;
    rho      = 1.0-eps[i];
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      k     = index[j];
      sum += vf[j]*x[k]*s[k];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      k     = index[j];
      sum += vf[j]*x[k]*s[k];
    }
    b[n] = row->area*s[i]*(x[i]*s[i]-sum*rho);
  }
}

void 
VF_RadSolveCleanUp(void)
{
  VF_Free(radq_g);
  VF_Free(tsurf_g);
  VF_Free(eps_g);
  VF_RadSolveCleanupCG();
  VF_RadSolveCleanupGMRES();
  VF_RadSolveCleanupAZTEC();
}
