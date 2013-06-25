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
@(#)    $RCSfile: VF_RadSolveGMRES.c,v $
@(#)    $Revision: 1.6 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_RadSolveGMRES.c,v $
@(#)
@(#)    DESCRIPTION:  GMRES solver to solve the radiosity matrix equation.
@(#)    The non-Aztec version is based on the GMRES solver from the KRYSOLV
@(#)    library.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

static int    kspace=20; 
static double **ss=NULL;
static double **hh=NULL;
static double *c=NULL;
static double *r=NULL;
static double *s=NULL;
static double *rs=NULL; 
static double *dots=NULL;
static double *temp=NULL;
  
void
VF_RadSolveCleanupGMRES(void)
{
  int i;
  VF_Free(c);
  VF_Free(r);
  VF_Free(s);
  VF_Free(rs);
  VF_Free(dots);
  VF_Free(temp);
  if (hh!=NULL) {
    for (i=0; i<kspace+1; i++){
      VF_Free(hh[i]);
    }
    VF_Free(hh);
  }
  if (ss!=NULL) {
    for (i=0; i<kspace+1; i++){
      VF_Free(ss[i]);
    }
    VF_Free(ss);
  }
}

void
VF_RadSolveGMRES(double *radq, double *radj, double *tsurf, double *eps, 
                 double sigma, double tol, int max_iter, int *num_iter, int debug)
{
  int    i, ii, i1, j, k, k1, n_l, n_g;
  int    converged=FALSE;
  int    its=0, ortho_flag=1;
  double t4, qsum, L2_norm, L2_norm0, ro, scale, t;
  double start, end, elapsed, time1, *sol, *rhs;
  double gam, epsmac=1.0e-12, *ptr_ss, *ptr_r;
  double tmin, tmax, emin, emax;
  VFrow  *row;
  static int g_size=0, l_size=0;
  VFenclosure *enclosure=VF_CurrentEnclosure();

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
      printf("   Solving radiosity equations with GMRES for enclosureID <%s>\n",
             enclosure->id);
    }
  }
  n_g = enclosure->npatches_g;
  n_l = enclosure->npatches_l;
  VF_GetDPbuffer0_ptr(&sol);
  VF_GetDPbuffer1_ptr(&rhs);
  /*============================*/
  /* ALLOCATE SOME LOCAL MEMORY */
  /*============================*/
  if (dots==NULL) {
    c    = VF_Newd(kspace);
    s    = VF_Newd(kspace);
    dots = VF_Newd(kspace+1);
    rs   = VF_Newd(kspace+1);
    ss   = VF_NewdPtr(kspace+1);
    hh   = VF_NewdPtr(kspace+1);
    for (j=0; j<kspace+1; j++){
      hh[j] = VF_Newd(kspace);
    }
  }
  if (g_size==0) {
    g_size = n_g;
    temp   = VF_Newd(n_g);
  } else if (n_g>g_size) {
    g_size = n_g;
    temp   = VF_ReNewd(temp,n_g);
  }
  if (l_size==0) {
    l_size = n_l;
    for (j=0; j<kspace+1; j++){
      ss[j] = VF_Newd(n_l);
    }
    r = VF_Newd(n_l);
  } else if (n_l>l_size) {
    l_size = n_l;
    for (j=0; j<kspace+1; j++){
      ss[j] = VF_ReNewd(ss[j],n_l);
    }
    r = VF_ReNewd(r,n_l);
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
  /*====================================*/
  /* INITIALIZE SOLUTION AND RHS VECTOR */
  /*====================================*/
  row = enclosure->row;
  for (i=0; i<n_l; i++, row++) {
    ii = row->global_index;
    t4 = sigma*tsurf[ii]*tsurf[ii]*tsurf[ii]*tsurf[ii];
    if (eps[ii]==0.0) {
      sol[ii] = t4;
      rhs[ii] = 0.0;
    } else {
      sol[ii] = t4+radq[ii]*(1.0-eps[ii])/eps[ii];
      rhs[ii] = eps[ii]*t4;
    }
    if (debug>=VF_OUTPUT_DEBUG_0) {
      printf("%d (%d):  %g  %g  %g  %g  %g\n",i,ii,tsurf[ii],eps[ii],radq[ii],rhs[ii],sol[ii]);
    }
  }
  VF_ExchangeDouble(rhs);
  VF_ExchangeDouble(sol);
  /*============================*/
  /* CALCULATE INITIAL RESIDUAL */
  /*============================*/
  VF_RadSolve_Residual(eps,rhs,sol,r); 
  L2_norm0 = sqrt(Ddot(n_l,r,r));
  L2_norm  = L2_norm0;
  if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
    printf("\n   Initial residual = %e    tol = %e\n",L2_norm0,tol);
    if (L2_norm0>=tol) {
      printf("      iter: %4d   residual = %.6e\n",its,L2_norm);
    }
  }
  if (L2_norm0<tol) converged = TRUE;
  for (i=0; i<n_l; i++) {
    r[i] = rhs[enclosure->row[i].global_index];
  }
  /*==========================*/
  /* OUTER RESTART GMRES LOOP */
  /*==========================*/
  if (debug>=VF_OUTPUT_SUMMARY) {
    time1 = VF_Clock();
  }
  while (converged!=TRUE && (its<max_iter)) {
    VF_RadSolve_MatVecMul(eps,sol,ss[0]);      
    ptr_ss = ss[0];
    ptr_r  = r;
    for (j=n_l; j>0; j--){
      *ptr_ss = *ptr_r++ - *ptr_ss;
      ptr_ss++;
    }
    /*================================================*/
    /* INITIAL RESIDUAL & SCALED CONVERGENCE CRITERIA */
    /*================================================*/ 
    ro = sqrt(Ddot(n_l,ss[0],ss[0]));
    if (ro==0.0e+0) break;
    scale = 1.0 / ro;
    Dscal(n_l,scale,ss[0]);
    /*=================================================*/ 
    /* INITIALIZE 1ST TERM OF RHS OF HESSENBERG SYSTEM */
    /*=================================================*/ 
    rs[0] = ro;
    i     = 0;
    /*====================*/
    /* INNER ARNOLDI LOOP */
    /*====================*/
    do {
      its++;
      i1 = i + 1;
      for (j=0; j<n_l; j++) {
        temp[enclosure->row[j].global_index] = ss[i][j];
      }
      VF_RadSolve_MatVecMul(eps,temp,ss[i1]);
     /*=================================*/
      /* GRAMM SCHMIDT ORTHOGONILIZATION */
      /*=================================*/
      if (ortho_flag==0) {
        /*==============*/
        /* CLASSICAL GS */
        /*==============*/
        for (j=0; j<=i; j++) {
          dots[j]  = Ddot(n_l,ss[j],ss[i1]);
          hh[j][i] = dots[j];
          t        = -dots[j] ;
          Daxpy(n_l,t,ss[j],ss[i1]);
        }
      } else if (ortho_flag==1) {
        /*=============*/
        /* MODIFIED GS */
        /*=============*/
        for (j=0; j<=i; j++) {
          t        = Ddot(n_l,ss[j],ss[i1]);
          hh[j][i] = t;
          Daxpy(n_l,-t,ss[j],ss[i1]);
        }
      }
      /*==================*/
      /* NORMALIZE VECTOR */
      /*==================*/
      t         = sqrt(Ddot(n_l,ss[i1],ss[i1]));
      hh[i1][i] = t;
      scale     = 1.0/t;
      Dscal(n_l,scale,ss[i1]);
      /*==============================================*/
      /* UPDATE FACTORIZATION OF HH BY PLANE ROTATION */
      /*==============================================*/
      if (i!=0) {
        for (k=1; k<=i; k++) {
          k1        = k-1;
          t         = hh[k1][i];
          hh[k1][i] = c[k1]*t+s[k1]*hh[k][i];
          hh[k][i]  = -s[k1]*t+c[k1]*hh[k][i];
        }
      }
      gam = sqrt(hh[i][i]*hh[i][i]+hh[i1][i]*hh[i1][i]);
      if (gam==0.0e+0) gam = epsmac;
      /*===============================*/
      /* DETERMINE NEXT PLANE ROTATION */
      /*===============================*/
      c[i]   = hh[i][i]/gam;
      s[i]   = hh[i1][i]/gam;
      rs[i1] = -s[i]*rs[i];
      rs[i]  = c[i]*rs[i];
      /*============================================*/
      /* DETERMINE RESIDUAL NORM & TEST CONVERGENCE */
      /*============================================*/
      hh[i][i] = c[i]*hh[i][i] + s[i]*hh[i1][i];
      ro       = fabs(rs[i1]);
      i++;       /*subspace dim. counter dim(K)=i-1 */
      if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
        printf("      iter: %4d   residual = %e\n",its,ro);
      }
    } while ((i<kspace) && (ro>tol) && (its<max_iter));
    /*=================================================*/
    /* SOLVE UPPER TRIANGULAR SYSTEM, COMPUTE SOLUTION */
    /*=================================================*/
    i--;          /* set i = Krylov subspace dimension */
    rs[i] /= hh[i][i];
    for (ii=1; ii<=i; ii++) {
      k  = i - ii;
      k1 = k + 1;
      t  = rs[k];
      for (j=k1; j<=i; j++){
        t -= hh[k][j]*rs[j];
      }
      rs[k] = t/hh[k][k];
    }
    /*====================================*/
    /* DONE WITH BACK SUBSTITUTION, FORM  */
    /* LINEAR COMBINATION TO GET SOLUTION */
    /*====================================*/
    for (j=0; j<=i; j++) {
      t = rs[j];
      Daxpy2(enclosure,t,ss[j],sol);
    }
    if (ro<tol) {
      VF_ExchangeDouble(sol);
      VF_RadSolve_Residual(eps,rhs,sol,ss[0]);
      L2_norm = sqrt(Ddot(n_l,ss[0],ss[0]));
      if (L2_norm<tol) {
        if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
          printf("\n      Real residual is %12.5e\n", L2_norm);
        }
        converged = TRUE;
      } else{
        if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
          printf("\n      Real residual check failed (%12.5e)\n",
                 L2_norm/L2_norm0);
        }
        converged = FALSE;
      }
    }
  }
    
  if (its==max_iter) {
    VF_RadSolve_Residual(eps,rhs,sol,ss[0]);
    L2_norm = sqrt(Ddot(n_l,ss[0],ss[0]));
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
    VF_ComputeFluxes(radq, radj, sol, &qsum, debug);
  }
  /*==========================================*/
  /* PRINT SOME STATISTICS IF IN VERBOSE MODE */
  /*==========================================*/
  if (debug>=VF_OUTPUT_SUMMARY) {
    end     = VF_Clock();
    time1   = time1-start;
    elapsed = end-start;
    if (VFLIB_Rank==0) {
      printf("\n");
      printf("   Total iterations:    %d\n",*num_iter);
      printf("   Total elapsed time:  %.2f (sec.)\n",elapsed);
      printf("     GMRES setup time:  %.2f (sec.)\n",time1);
      printf("     GMRES solve time:  %.2f (sec.)\n",elapsed-time1);
      printf("   Residual L2 norm:    %g\n",L2_norm);
      printf("   Flux Integration:    %g\n",qsum);
      printf("   Minimum Emissivity:  %g\n",emin);
      printf("   Maximum Emissivity:  %g\n",emax);
      printf("   Minimum Temperature: %g\n",tmin);
      printf("   Maximum Temperature: %g\n",tmax);
      printf("\n\n");
    }
  }
  if (*num_iter==max_iter) *num_iter = -max_iter;
}
