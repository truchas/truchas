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
@(#)    $RCSfile: VF_RadSolveCG.c,v $
@(#)    $Revision: 1.7 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_RadSolveCG.c,v $
@(#)
@(#)    DESCRIPTION:  CG solver to solve the radiosity matrix equation.
@(#)    The non-Aztec version is based on the CG solver from the KRYSOLV
@(#)    library.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

static double *z=NULL;
static double *r=NULL;
static double *p=NULL;
static double *ap=NULL;
static double *scale=NULL;

void
VF_RadSolveCleanupCG(void)
{
  VF_Free(z);
  VF_Free(r);
  VF_Free(p);
  VF_Free(ap);
  VF_Free(scale);
}

void
VF_RadSolveCG(double *radq, double *radj, double *tsurf, double *eps, 
              double sigma, double tol, int max_iter, int *num_iter, int debug)
{
  int    i, ii, n_l, n_g;
  int    converged=FALSE, its=0;
  double t4, flux, L2_norm, L2_norm0;
  double start, end, elapsed, time1, *sol, *rhs;
  double beta, p_ap_dot, alpha, nalpha, r_z_dot, r_z_dot_old;
  double tmin, tmax, emin, emax;
  VFrow  *row;
  static int g_size=0, l_size=0;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  /*=====================================================*/
  /* SOLVE THE RADIOSITY MATRIX EQUATION WITH PCG SOLVER */
  /*=====================================================*/
  if (debug>=VF_OUTPUT_SUMMARY) {
    start = VF_Clock();
    if (VFLIB_Rank==0) {
      printf("\n\n");
      printf("  **********************************");
      printf("**********************************\n");
      printf("   R A D I O S I T Y    M A T R I X    S O L V E R\n");
      printf("  **********************************");
      printf("**********************************\n");
      printf("\n");
      printf("   Solving radiosity equations with CG for enclosureID <%s>\n",enclosure->id);
    }
  }
  n_g = enclosure->npatches_g;
  n_l = enclosure->npatches_l;
  VF_GetDPbuffer0_ptr(&sol);
  VF_GetDPbuffer1_ptr(&rhs);
  /*============================*/
  /* ALLOCATE SOME LOCAL MEMORY */
  /*============================*/
  if (g_size==0) {
    g_size = n_g;
    p      = VF_Newd(n_g);
    scale  = VF_Newd(n_g);
  } else if (n_g>g_size) {
    g_size = n_g;
    p      = VF_ReNewd(p, n_g);
    scale  = VF_ReNewd(scale, n_g);
  }
  if (l_size==0) {
    l_size = n_l;
    z      = VF_Newd(n_l);
    r      = VF_Newd(n_l);
    ap     = VF_Newd(n_l);
  } else if (n_l>l_size) {
    l_size = n_l;
    z      = VF_ReNewd(z,  n_l);
    r      = VF_ReNewd(r,  n_l);
    ap     = VF_ReNewd(ap, n_l);
  }
  /*============================================*/
  /* FIND MIN/MAX TEMPERATURES AND EMMISIVITIES */
  /*============================================*/
  tmin = tsurf[0];
  tmax = tsurf[0];
  emin = eps[0];
  emax = eps[0];
  for (i=0; i<n_g; i++) {
    tmin = MIN(tmin,tsurf[i]);
    tmax = MAX(tmax,tsurf[i]);
    emin = MIN(emin,eps[i]);
    emax = MAX(emax,eps[i]);
  }
  /*=================================*/
  /* INITIALIZE SOLUTION, RHS VECTOR */
  /*=================================*/
  row = enclosure->row;
  for (i=0; i<n_l; i++, row++) {
    ii   = row->global_index;
    t4   = sigma*tsurf[ii]*tsurf[ii]*tsurf[ii]*tsurf[ii];
    if (eps[ii]==0.0) {
      sol[ii] = t4;
      rhs[ii] = 0.0;
    } else {
      sol[ii] = t4+radq[ii]*(1.0-eps[ii])/eps[ii];
      rhs[ii] = eps[ii]*t4*row->area;
    }
    if (debug>=VF_OUTPUT_DEBUG_0) {
      printf("%d (%d):  %g  %g  %g  %g  %g\n",i,ii,tsurf[ii],eps[ii],radq[ii],rhs[ii],sol[ii]);
    }
    scale[ii] = 1.0/sqrt(row->area);
    rhs[ii]  *= scale[ii];
    sol[ii]  /= scale[ii];
  }
  VF_ExchangeDouble(rhs);
  VF_ExchangeDouble(sol);
  VF_ExchangeDouble(scale);
  /*============================*/
  /* CALCULATE INITIAL RESIDUAL */
  /*============================*/
  VF_RadSolve_Residual_Scale(eps,rhs,sol,r,scale);  
  L2_norm0 = sqrt(Ddot(n_l,r,r));
  L2_norm  = L2_norm0;
  for (i=0; i<n_l; i++) {
    z[i] = r[i];
  }
  r_z_dot = Ddot(n_l,r,z);
  if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
    printf("\n   Initial residual = %e    tol = %e\n",L2_norm0,tol);
    if (L2_norm0>=tol) {
      printf("      iter: %4d   residual = %.6e\n",its,L2_norm);
    }
  }
  if (L2_norm0<tol) converged = TRUE;
  /*================*/
  /* ITERATION LOOP */
  /*================*/
  if (debug>=VF_OUTPUT_SUMMARY) {
    time1 = VF_Clock();
  }
  beta = 0.0;
  for (i=0; i<n_g; i++) {
    p[i] = 0.0;
  }
  while (!converged && (its<max_iter)) {
    its++;
    row = enclosure->row;
    for (i=0; i<n_l; i++, row++) {
      ii    = row->global_index;
      p[ii] = z[i] + beta*p[ii];
    } 
    VF_RadSolve_MatVecMul_Scale(eps,p,ap,scale); 
    p_ap_dot = Ddot1(enclosure,p,ap);
    alpha    = r_z_dot/p_ap_dot; 
    nalpha   = -alpha;
    Daxpy1(enclosure,alpha,p,sol); 
    Daxpy(n_l,nalpha,ap,r);
    Dcopy(n_l,r,z);
    r_z_dot_old = r_z_dot;
    r_z_dot     = Ddot(n_l,r,z);
    beta        = r_z_dot/r_z_dot_old;
    L2_norm     = sqrt(Ddot(n_l,r,r));
    if (L2_norm<tol) converged = TRUE;
    if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
      printf("      iter: %4d   residual = %.6e\n",its,L2_norm);
      fflush(stdout);
    }
  }
  *num_iter = its;
  /*=========================================================*/
  /* FINISHED WITH RADIOSITY MATRIX SOLVE, CHECK CONVERGENCE */
  /*=========================================================*/
  if (*num_iter==max_iter) {
    if (VFLIB_Rank==0) {
      printf("\n");
      printf("   *** WARNING ***\n");
      printf("   NO CONVERGENCE FOR RADIATION PROBLEM\n");
      printf("   IN ENCLOSURE <%s> AFTER %d ITERATIONS\n\n",
             enclosure->id,max_iter);
    }
  } else if (*num_iter<0) { 
    if (VFLIB_Rank==0) {
      printf("\n");
      printf("   *** WARNING ***\n");
      printf("   DIVERGING SOLUTION FOR RADIATION PROBLEM\n");
      printf("   IN ENCLOSURE <%s> AFTER %d ITERATIONS\n\n",
             enclosure->id,-(*num_iter));
    }
  } else {
    /*====================================*/
    /* RADIOSITY MATRIX SOLVE SUCCESSFUL, */
    /* COMPUTE THE FINAL HEAT FLUXES      */
    /*====================================*/
    row = enclosure->row;
    for (i=0; i<n_l; i++, row++) {
      ii       = row->global_index;
      sol[ii] *= scale[ii];
    } 
    VF_ExchangeDouble(sol);
    VF_ComputeFluxes(radq, radj, sol, &flux, debug);
  }
  /*==========================================*/
  /* PRINT SOME STATISTICS IF IN VERBOSE MODE */
  /*==========================================*/
  if (debug>=VF_OUTPUT_SUMMARY && VFLIB_Rank==0) {
    end     = VF_Clock();
    time1   = time1-start;
    elapsed = end-start;
    printf("\n");
    printf("   Total iterations:    %d\n",*num_iter);
    printf("   Total elapsed time:  %.2f (sec.)\n",elapsed);
    printf("     CG setup time:     %.2f (sec.)\n",time1);
    printf("     CG solve time:     %.2f (sec.)\n",elapsed-time1);
    printf("   Residual L2 norm:    %g\n",L2_norm);
    printf("   Flux Integration:    %g\n",flux);
    printf("   Minimum Emissivity:  %g\n",emin);
    printf("   Maximum Emissivity:  %g\n",emax);
    printf("   Minimum Temperature: %g\n",tmin);
    printf("   Maximum Temperature: %g\n",tmax);
    printf("\n\n");
  }
  if (*num_iter==max_iter) *num_iter = -max_iter;
}
