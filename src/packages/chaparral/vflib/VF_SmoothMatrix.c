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
@(#)    $RCSfile: VF_SmoothMatrix.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_SmoothMatrix.c,v $
@(#)
@(#)    DESCRIPTION:  Least squares smoothing of the viewfactor matrix
@(#)    using lagrange multipliers and a CG iterative solver with Jacobi
@(#)    preconditioning.  An iterative solver with the weighting matrix 
@(#)    built "on-the-fly" is used instead of a direct solver so that
@(#)    additional memory doesn't have to be allocated to store the 
@(#)    weighting matrix.
@(#)
@(#)    REFERENCE:    "Least-Squares Smoothing of Direct-Exchange Areas in
@(#)                  Zonal Analysis," M.E. Larsen and J.R. Howell, J. Heat
@(#)                  Transfer, 108, 1986.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vf.h"

void   VF_MatrixStats(int*, int*, int*);
void   ResidualSM(double wt, double d[], 
                  double b[], double u[], double r[]);
void   MatVecMulSM0(double wt, 
                    double d[], double x[], double *y);
void   MatVecMulSM1(double wt, 
                    double d[], double x[], double *y);
void   MatVecMulSM2(double wt, 
                    double d[], double x[], double *y);
void   MatVecMulSM3(double wt, 
                    double d[], double x[], double *y);
void   MatVecMulSM4(double wt, 
                    double d[], double x[], double *y);
int    UpdateMatrix0(double wt, double s[]);
int    UpdateMatrix1(double wt, double s[]);
int    UpdateMatrix2(double wt, double s[]);
int    UpdateMatrix3(double wt, double s[]);
int    UpdateMatrix4(double wt, double s[]);
void   ScaleSM0(double scale[], double wt);
void   ScaleSM1(double scale[], double wt);
void   ScaleSM2(double scale[], double wt);
void   ScaleSM3(double scale[], double wt);
void   ScaleSM4(double scale[], double wt);
void   ZeroNegativeEntries(VFrow*);

typedef void ScaleSM_FN(double*, double);
typedef void MatVecMulSM_FN(double, double*, double*, double*);
typedef int  UpdateMatrix_FN(double, double*);

static ScaleSM_FN      *ScaleSM;
static MatVecMulSM_FN  *MatVecMulSM;
static UpdateMatrix_FN *UpdateMatrix;

void VF_SmoothMatrix(int encl, double wt, double tol, int max_iter, 
                     int symmetric, int output)
{
  int    i, ii, nrows_g, nrows_l, sym_method, bad=0;
  int    pass=0, max_passes=5, done=0;
  int    knt0, knt1, knt2, converged, iter;
  double L2_norm0, L2_norm, rowsum, area;
  double beta, p_ap_dot, alpha, nalpha;
  double *sol, *rhs, r_z_dot, r_z_dot_old;
  double mean0, variance0, mean1, variance1, mean2, variance2;
  double *z, *r, *p, *ap, *scale;
  double sum00, sum01, sum02;
  double sum0, min0, max0, sd0;
  double sum1, min1, max1, sd1;
  double sum2, min2, max2, sd2;
  double clock0, clock1;
  VFrow  *row;
  VFenclosure *enclosure=VF_GetEnclosure(encl);
    
  clock0 = VF_Clock();
  VF_OutputSmoothingBanner(wt, tol, max_iter);
  if (fabs(wt-floor(wt+0.0000001))<0.0000001) {
    ii = (int)(wt+0.1);
  } else {
    ii = 0;
  }
  switch (ii) {
  case 1:
    ScaleSM      = (void (*)(double *,double))ScaleSM1;
    MatVecMulSM  = (void (*)(double ,double *,double *,double *))MatVecMulSM1;
    UpdateMatrix = (int  (*)(double ,double *))UpdateMatrix1;
    break;
  case 2:
    ScaleSM      = (void (*)(double *,double))ScaleSM2;
    MatVecMulSM  = (void (*)(double ,double *,double *,double *))MatVecMulSM2;
    UpdateMatrix = (int  (*)(double ,double *))UpdateMatrix2;
    break;
  case 3:
    ScaleSM      = (void (*)(double *,double))ScaleSM3;
    MatVecMulSM  = (void (*)(double ,double *,double *,double *))MatVecMulSM3;
    UpdateMatrix = (int  (*)(double ,double *))UpdateMatrix3;
    break;
  case 4:
    ScaleSM      = (void (*)(double *,double))ScaleSM4;
    MatVecMulSM  = (void (*)(double ,double *,double *,double *))MatVecMulSM4;
    UpdateMatrix = (int  (*)(double ,double *))UpdateMatrix4;
    break;
  default:
    ScaleSM      = (void (*)(double *,double))ScaleSM0;
    MatVecMulSM  = (void (*)(double ,double *,double *,double *))MatVecMulSM0;
    UpdateMatrix = (int  (*)(double ,double *))UpdateMatrix0;
    break;
  }
  /*==============================================*/
  /* USE JACOBI PRECONDITIONED CONJUGATE GRADIENT */
  /* (PCG) ITERATIVE SOLVER TO SOLVE THE MATRIX   */
  /*==============================================*/
  if (enclosure->smoothed) {
    if (output>=VF_OUTPUT_SUMMARY && VFLIB_Rank==0) {
      printf("     Matrix has already been previously smoothed\n");
    }
    return;
  }
  VF_GetDPbuffer0_ptr(&rhs);
  VF_GetDPbuffer1_ptr(&sol);
  nrows_g    = enclosure->npatches_g;
  nrows_l    = enclosure->npatches_l;
  z          = VF_Newd(nrows_l);
  r          = VF_Newd(nrows_l);
  p          = VF_Newd(nrows_g);
  ap         = VF_Newd(nrows_l);
  scale      = VF_Newd(nrows_g);
  sym_method = symmetric;
  if (output>=VF_OUTPUT_VERBOSE) {
    for (row=enclosure->row, sum00=0.0, i=0; i<nrows_l; i++, row++) {
      sum00 += row->raw_rowsum;
    }
    VF_GlobalSumDouble(&sum00);
    VF_ReciprocityStats(&mean0, &variance0);
  }
  /*=======================================================*/
  /* MAKE SURE THE MATRIX IS SYMMETRIC WITH RESPECT TO     */
  /* ITS SPARSITY PATTERN.  NEED TO DO THIS BECAUSE THE    */
  /* STOCHASTIC BEHAVIOR OF THE HEMICUBE METHOD (JITTERING */
  /* AND ALIASING) AND RAYCAST METHOD CAN LEAD TO CASES    */
  /* WHERE F(I,J) AND F(J,I) ARE NOT BOTH ZERO OR BOTH     */
  /* NONZERO.  IN THESE CASES, F(I,J) WILL BE VERY CLOSE   */
  /* TO ZERO SO WE WILL JUST ZERO OUT THE NONZERO TERM.    */
  /*=======================================================*/
  VF_MakeMatrixSymmetric(encl, sym_method, output);
  if (output>=VF_OUTPUT_VERBOSE) {
    for (row=enclosure->row, sum01=0.0, i=0; i<nrows_l; i++, row++) {
      sum01 += row->sym_rowsum;
    }
    VF_GlobalSumDouble(&sum01);
    VF_ReciprocityStats(&mean1, &variance1);
  }
  /*============================*/
  /* INITIALIZE SOLUTION VECTOR */
  /*============================*/
  for (row=enclosure->row, i=0; i<nrows_l; i++, row++) {
    ii      = row->global_index;
    area    = row->area;
    rowsum  = row->sym_rowsum;
    sol[ii] = 0.0;
  }
  do {
    iter      = 0;
    converged = FALSE;
    /*===========================================*/
    /* FIND JACOBI PRECONDITIONING SCALE FACTORS */
    /*===========================================*/
    ScaleSM(scale, wt);
    /*=====================================*/
    /* INITIALIZE RHS AND SOLUTION VECTORS */
    /*=====================================*/
    for (row=enclosure->row, i=0; i<nrows_l; i++, row++) {
      ii       = row->global_index;
      area     = row->area;
      rowsum   = row->sym_rowsum;
      rhs[ii]  = area*(1.0-rowsum)*scale[ii];
      sol[ii] /= scale[ii];
    }
    VF_ExchangeDouble(rhs);
    VF_ExchangeDouble(sol);
    VF_ExchangeDouble(scale);
    /*=================================*/
    /* CALCULATE INITIAL TRUE RESIDUAL */
    /*=================================*/
    ResidualSM(wt,scale,rhs,sol,r);
    for (row=enclosure->row, i=0; i<nrows_l; i++, row++) {
      ii    = row->global_index;
      ap[i] = r[i]/scale[ii];
    }
    L2_norm0 = sqrt(Ddot(nrows_l,ap,ap));
    L2_norm  = L2_norm0;
    for (i=0; i<nrows_l; i++) {
      z[i] = r[i];
    }
    r_z_dot = Ddot(nrows_l,r,z);
    if (output>=VF_OUTPUT_DEBUG_0 && VFLIB_Rank==0) {
      printf("\n");
      printf("\t\titer: %4d\tresidual = %.6e\n",iter,L2_norm0);
      fflush(stdout);
    }
    if (L2_norm0<tol) converged = TRUE;
    /*================*/
    /* ITERATION LOOP */
    /*================*/
    beta = 0.0;
    while (converged!=TRUE && (iter<max_iter)) {
      iter++;
      for (row=enclosure->row, i=0; i<nrows_l; i++, row++) {
        ii    = row->global_index;
        p[ii] = z[i] + beta*p[ii];
      }
      MatVecMulSM(wt,scale,p,ap);
      p_ap_dot = Ddot1(enclosure,p,ap);
      alpha    = r_z_dot/p_ap_dot;
      nalpha   = -alpha;
      Daxpy1(enclosure,alpha,p,sol);
      Daxpy(nrows_l,nalpha,ap,r);
      Dcopy(nrows_l,r,z);
      r_z_dot_old = r_z_dot;
      r_z_dot     = Ddot(nrows_l,r,z);
      beta        = r_z_dot/r_z_dot_old;
      for (row=enclosure->row, i=0; i<nrows_l; i++, row++) {
        ii    = row->global_index;
        ap[i] = r[i]/scale[ii];
      }
      L2_norm = sqrt(Ddot(nrows_l,ap,ap));
      if (L2_norm<tol) {
        for (row=enclosure->row, i=0; i<nrows_l; i++, row++) {
          ii       = row->global_index;
          sol[ii] *= scale[ii];
        }
        VF_ExchangeDouble(sol);
        converged = TRUE;
      }
      if (output>=VF_OUTPUT_DEBUG_0 && VFLIB_Rank==0) {
        printf("\t\titer: %4d\tresidual = %.6e\n",iter,L2_norm);
        fflush(stdout);
      }
    }
    /*===============================================*/
    /* FINISHED WITH MATRIX SOLVE, CHECK CONVERGENCE */
    /*===============================================*/
    if (iter==max_iter) {
      /*=======================================================*/
      /* NO CONVERGENCE, MAXIMUM NUMBER OF ITERATIONS EXCEEDED */
      /*=======================================================*/
      done = 1;
      if (VFLIB_Rank==0) {
        printf("\n");
        printf("\t*** WARNING ***\n");
        printf("\tNO CONVERGENCE FOR VF MATRIX SMOOTHING\n");
        printf("\tIN ENCLOSURE %d AFTER %d ITERATIONS,\n",
               encl,max_iter);
        printf("\tNO SMOOTHING PERFORMED\n");
        printf("\tTolerance = %g\n",tol);
        printf("\tL2_norm0  = %g\n",L2_norm0);
        printf("\tL2_norm   = %g\n",L2_norm);
        printf("\n");
      }
    } else if (iter==0) {
      /*======================================*/
      /* CONVERGENCE, NO ITERATIONS PERFORMED */
      /*======================================*/
      done = 1;
      if (VFLIB_Rank==0) {
        printf("     no smoothing performed, initial residual (%g) < tol\n",
               L2_norm0);
      }
    } else {
      /*=============*/
      /* CONVERGENCE */
      /*=============*/
      done = 1;
      /*===========================*/
      /* SET NEW VIEWFACTOR VALUES */
      /*===========================*/
      bad = UpdateMatrix(wt, sol);
      if (bad>0) {
        if (output>=VF_OUTPUT_VERBOSE) {
          VF_MatrixStats(&knt0, &knt1, &knt2);
          if (VFLIB_Rank==0) {
            if (pass==0) printf("\n");
            printf("     Pass %d ..........\n",pass);
            printf("       iterations = %d\n",iter);
            printf("       L2_norm0   = %g\n",L2_norm0);
            printf("       L2_norm    = %g\n",L2_norm);
            printf("\n");
            printf("       %d negative viewfactors detected after smoothing\n",bad);
            if (pass<max_passes) {
              printf("       zeroing out negative viewfactors and making another pass\n");
            } else {
              printf("       maximum number of passes (%d) exceeded\n",max_passes);
            }  
            printf("       Nonzero lower triangular entries = %d\n",knt1);
            printf("       Nonzero upper triangular entries = %d\n",knt2);
            printf("\n");
            fflush(stdout);
          }
        }
        pass++;
        if (pass>=max_passes) {
          done = 1;
        } else {
          done = 0;
        }
      } else {
        pass++;
      }
    }
  } while (!done);
  if (output>=VF_OUTPUT_SUMMARY) {
    if (output>=VF_OUTPUT_VERBOSE) {
      for (sum02=0.0, row=enclosure->row, i=0; i<nrows_l; i++, row++) {
        sum02 += row->smooth_rowsum;
      }
      VF_GlobalSumDouble(&sum02);
      VF_RowsumStats(&min0, &max0, &sum0, &sd0,
                     &min1, &max1, &sum1, &sd1,
                     &min2, &max2, &sum2, &sd2);
      VF_ReciprocityStats(&mean2, &variance2);
    }
    clock1 = VF_Clock();
    enclosure->time_vf_smooth = (clock1-clock0)-enclosure->time_vf_symmetry;
    if (VFLIB_Rank==0) {
      if (output>=VF_OUTPUT_SUMMARY) {
        printf("     Number of passes     = %d\n",pass);
        printf("     Number of iterations = %d\n",iter);
        printf("     Elapsed time         = %.2f\n",enclosure->time_vf_smooth);
      } else if (output>=VF_OUTPUT_VERBOSE) {
        printf("     row sum total should be                   = %.6f\n",(float)nrows_g);
        if (sym_method==VF_SYMMETRIC_NONE) {
          printf("     row sum total before smoothing            = %.6f\n",sum00);
        } else {
          printf("     row sum total before forcing symmetry     = %.6f\n",sum00);
          printf("     row sum total after forcing symmetry      = %.6f\n",sum01);
        }
        printf("     row sum total after smoothing             = %.6f\n",sum02);
        printf("     Raw rowsum error min                      = %11.4e\n",min0);
        printf("     Raw rowsum error max                      = %11.4e\n",max0);
        printf("     Raw rowsum error mean                     = %11.4e +/- %10.4e\n",sum0,sd0);
        if (sym_method!=VF_SYMMETRIC_NONE) {
          printf("     Symmetric rowsum error min                = %11.4e\n",min1);
          printf("     Symmetric rowsum error max                = %11.4e\n",max1);
          printf("     Symmetric rowsum error mean               = %11.4e +/- %10.4e\n",sum1,sd1);
        }
        printf("     Smoothed rowsum error min                 = %11.4e\n",min2);
        printf("     Smoothed rowsum error max                 = %11.4e\n",max2);
        printf("     Smoothed rowsum error mean                = %11.4e +/- %10.4e\n",sum2,sd2);
        if (sym_method==VF_SYMMETRIC_NONE) {
          printf("     reciprocity error before smoothing        = %11.4e +/- %10.4e\n",mean0,variance0);
        } else {
          printf("     reciprocity error before forcing symmetry = %11.4e +/- %10.4e\n",mean0,variance0);
          printf("     reciprocity error after forcing symmetry  = %11.4e +/- %10.4e\n",mean1,variance1);
        }
        printf("     reciprocity error after smoothing         = %11.4e +/- %10.4e\n",mean2,variance2);
        printf("     number of iterations to smooth            = %d\n",iter);
        printf("     true residual                             = %g\n",L2_norm);
        printf("     Number of smoothing passes                = %d\n",pass);
      }
      printf("\n");
      if (bad>0) {
        printf("     WARNING!!  Negative viewfactors detected!\n\n");
      }
      if (!converged) {
        printf("     WARNING!!  Solution is not converged!\n\n");
      }
    }
  }
  enclosure->smoothed = 1;
  /*=====================*/
  /* FREE SCRATCH MEMORY */
  /*=====================*/
  VF_Free(ap);
  VF_Free(p);
  VF_Free(z);
  VF_Free(r);
  VF_Free(scale);
}

void ScaleSM0(double scale[], double wt)
{
  int    i, j, ii, nrows_l, nonzeros;
  float  *vf;
  double vfa, vfi, area;
  double wii, wij, wsum;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  /*===========================================*/
  /* FIND JACOBI PRECONDITIONING SCALE FACTORS */
  /*===========================================*/
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (i=0; i<nrows_l; i++, row++) {
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = pow(vfa,wt);
      wsum += wij;
    }
    vf       = row->array1.data;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = pow(vfa,wt);
      wsum += wij;
    }
    vfa = vfi*area;
    wii = pow(vfa,wt)+wsum;
    if (wii==0.0) {
      scale[ii] = 1.0;
    } else {
      scale[ii] = 1.0/sqrt(wii);
    }
  }
}

void ScaleSM1(double scale[], double wt)
{
  int    i, j, ii, nrows_l, nonzeros;
  float  *vf;
  double vfa, vfi, area;
  double wii, wij, wsum;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  /*===========================================*/
  /* FIND JACOBI PRECONDITIONING SCALE FACTORS */
  /*===========================================*/
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (i=0; i<nrows_l; i++, row++) {
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa;
      wsum += wij;
    }
    vf       = row->array1.data;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa;
      wsum += wij;
    }
    vfa = vfi*area;
    wii = vfa+wsum;
    if (wii==0.0) {
      scale[ii] = 1.0;
    } else {
      scale[ii] = 1.0/sqrt(wii);
    }
  }
}

void ScaleSM2(double scale[], double wt)
{
  int    i, j, ii, nrows_l, nonzeros;
  float  *vf;
  double vfa, vfi, area;
  double wii, wij, wsum;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  /*===========================================*/
  /* FIND JACOBI PRECONDITIONING SCALE FACTORS */
  /*===========================================*/
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (i=0; i<nrows_l; i++, row++) {
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa;
      wsum += wij;
    }
    vf       = row->array1.data;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa;
      wsum += wij;
    }
    vfa = vfi*area;
    wii = vfa*vfa+wsum;
    if (wii==0.0) {
      scale[ii] = 1.0;
    } else {
      scale[ii] = 1.0/sqrt(wii);
    }
  }
}

void ScaleSM3(double scale[], double wt)
{
  int    i, j, ii, nrows_l, nonzeros;
  float  *vf;
  double vfa, vfi, area;
  double wii, wij, wsum;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  /*===========================================*/
  /* FIND JACOBI PRECONDITIONING SCALE FACTORS */
  /*===========================================*/
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (i=0; i<nrows_l; i++, row++) {
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa;
      wsum += wij;
    }
    vf       = row->array1.data;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa;
      wsum += wij;
    }
    vfa = vfi*area;
    wii = vfa*vfa*vfa+wsum;
    if (wii==0.0) {
      scale[ii] = 1.0;
    } else {
      scale[ii] = 1.0/sqrt(wii);
    }
  }
}

void ScaleSM4(double scale[], double wt)
{
  int    i, j, ii, nrows_l, nonzeros;
  float  *vf;
  double vfa, vfi, area;
  double wii, wij, wsum;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  /*===========================================*/
  /* FIND JACOBI PRECONDITIONING SCALE FACTORS */
  /*===========================================*/
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (i=0; i<nrows_l; i++, row++) {
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa*vfa;
      wsum += wij;
    }
    vf       = row->array1.data;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa*vfa;
      wsum += wij;
    }
    vfa = vfi*area;
    wii = vfa*vfa*vfa*vfa+wsum;
    if (wii==0.0) {
      scale[ii] = 1.0;
    } else {
      scale[ii] = 1.0/sqrt(wii);
    }
  }
}

void ResidualSM(double wt, double d[], 
                double rhs[], double sol[], double r[])
{
  int i;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  MatVecMulSM(wt, d, sol, r);
  for (i=0; i<enclosure->npatches_l; i++) {
    r[i] = rhs[enclosure->row[i].global_index] - r[i];
  }
}

void MatVecMulSM0(double wt, 
                  double d[], double x[], double b[])
{
  int    i, ii, j, n_l, nonzeros, *index;
  float  *vf, vfi;
  double sum, wsum, wii, wij, vfa, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (i=0; i<n_l; i++, row++) {
    sum  = 0.0;
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = pow(vfa,wt);
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = pow(vfa,wt);
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vfa  = vfi*area;
    wii  = pow(vfa,wt)+wsum;
    b[i] = sum+wii*d[ii]*x[ii]*d[ii];
  }
}

void MatVecMulSM1(double wt, 
                  double d[], double x[], double b[])
{
  int    i, ii, j, n_l, nonzeros, *index;
  float  *vf, vfi;
  double sum, wsum, wii, wij, vfa, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (i=0; i<n_l; i++, row++) {
    sum  = 0.0;
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vfa  = vfi*area;
    wii  = vfa+wsum;
    b[i] = sum+wii*d[ii]*x[ii]*d[ii];
  }
}

void MatVecMulSM2(double wt, 
                  double d[], double x[], double b[])
{
  int    i, ii, j, n_l, nonzeros, *index;
  float  *vf, vfi;
  double sum, wsum, wii, wij, vfa, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (i=0; i<n_l; i++, row++) {
    sum  = 0.0;
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vfa  = vfi*area;
    wii  = vfa*vfa+wsum;
    b[i] = sum+wii*d[ii]*x[ii]*d[ii];
  }
}

void MatVecMulSM3(double wt, 
                  double d[], double x[], double b[])
{
  int    i, ii, j, n_l, nonzeros, *index;
  float  *vf, vfi;
  double sum, wsum, wii, wij, vfa, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (i=0; i<n_l; i++, row++) {
    sum  = 0.0;
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vfa  = vfi*area;
    wii  = vfa*vfa*vfa+wsum;
    b[i] = sum+wii*d[ii]*x[ii]*d[ii];
  }
}

void MatVecMulSM4(double wt, 
                  double d[], double x[], double b[])
{
  int    i, ii, j, n_l, nonzeros, *index;
  float  *vf, vfi;
  double sum, wsum, wii, wij, vfa, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  n_l = enclosure->npatches_l;
  row = enclosure->row;
  VF_ExchangeDouble(x);
  for (i=0; i<n_l; i++, row++) {
    sum  = 0.0;
    wsum = 0.0;
    ii   = row->global_index;
    vfi  = row->diagonal;
    area = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa*vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa   = vf[j]*area;
      wij   = vfa*vfa*vfa*vfa;
      wsum += wij;
      if (index[j]!=ii) sum += wij*d[ii]*x[index[j]]*d[index[j]];
    }
    vfa  = vfi*area;
    wii  = vfa*vfa*vfa*vfa+wsum;
    b[i] = sum+wii*d[ii]*x[ii]*d[ii];
  }
}

int UpdateMatrix0(double wt, double s[])
{
  int    i, ii, j, bad=0, bad1;
  int    nrows_l, nonzeros, *index;
  float  *vf;
  double vfa, wij, rsum, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (ii=0; ii<nrows_l; ii++, row++) {
    bad1     = 0;
    rsum     = 0.0;
    i        = row->global_index;
    area     = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = pow(vfa,wt);
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = pow(vfa,wt);
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    if (bad1) ZeroNegativeEntries(row);
    row->smooth_rowsum = rsum;
    bad += bad1;
  }
  VF_GlobalSumInt(&bad);
  return bad;
}


int UpdateMatrix1(double wt, double s[])
{
  int    i, ii, j, bad=0, bad1;
  int    nrows_l, nonzeros, *index;
  float  *vf;
  double vfa, wij, rsum, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (ii=0; ii<nrows_l; ii++, row++) {
    bad1     = 0;
    rsum     = 0.0;
    i        = row->global_index;
    area     = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    if (bad1) ZeroNegativeEntries(row);
    row->smooth_rowsum = rsum;
    bad += bad1;
  }
  VF_GlobalSumInt(&bad);
  return bad;
}

int UpdateMatrix2(double wt, double s[])
{
  int    i, ii, j, bad=0, bad1;
  int    nrows_l, nonzeros, *index;
  float  *vf;
  double vfa, wij, rsum, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (ii=0; ii<nrows_l; ii++, row++) {
    bad1     = 0;
    rsum     = 0.0;
    i        = row->global_index;
    area     = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa*vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum += vf[j];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa*vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum += vf[j];
    }
    if (bad1) ZeroNegativeEntries(row);
    row->smooth_rowsum = rsum;
    bad += bad1;
  }
  VF_GlobalSumInt(&bad);
  return bad;
}

int UpdateMatrix3(double wt, double s[])
{
  int    i, ii, j, bad=0, bad1;
  int    nrows_l, nonzeros, *index;
  float  *vf;
  double vfa, wij, rsum, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (ii=0; ii<nrows_l; ii++, row++) {
    bad1     = 0;
    rsum     = 0.0;
    i        = row->global_index;
    area     = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa*vfa*vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa*vfa*vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    if (bad1) ZeroNegativeEntries(row);
    row->smooth_rowsum = rsum;
    bad += bad1;
  }
  VF_GlobalSumInt(&bad);
  return bad;
}

int UpdateMatrix4(double wt, double s[])
{
  int    i, ii, j, bad=0, bad1;
  int    nrows_l, nonzeros, *index;
  float  *vf;
  double vfa, wij, rsum, area;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  nrows_l = enclosure->npatches_l;
  row     = enclosure->row;
  for (ii=0; ii<nrows_l; ii++, row++) {
    bad1     = 0;
    rsum     = 0.0;
    i        = row->global_index;
    area     = row->area;
    vf       = row->array0.data;
    index    = row->array0.index;
    nonzeros = row->array0.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa*vfa*vfa*vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    vf       = row->array1.data;
    index    = row->array1.index;
    nonzeros = row->array1.cnt;
    for (j=0; j<nonzeros; j++) {
      vfa    = vf[j]*area;
      wij    = vfa*vfa*vfa*vfa;
      vf[j] += (s[i]+s[index[j]])*wij/area;
      if (vf[j]<0.0) {
        vf[j] = 0.0;
        bad1++;
      }
      rsum  += vf[j];
    }
    if (bad1) ZeroNegativeEntries(row);
    row->smooth_rowsum = rsum;
    bad += bad1;
  }
  VF_GlobalSumInt(&bad);
  return bad;
}

void
ZeroNegativeEntries(VFrow *row)
{
  int j, k, cnt, nonzeros, *index;
  float *vf;

  cnt      = 0;
  k        = 0;
  vf       = row->array0.data;
  index    = row->array0.index;
  nonzeros = row->array0.cnt;
  for (j=0; j<nonzeros; j++) {
    if (vf[j]!=0.0) {
      if (j!=k) {
        vf[k]    = vf[j];
        index[k] = index[j];
      }
      cnt++;
      k++;
    }
  }
  row->array0.cnt = cnt;
            
  cnt      = 0;
  k        = 0;
  vf       = row->array1.data;
  index    = row->array1.index;
  nonzeros = row->array1.cnt;
  for (j=0; j<nonzeros; j++) {
    if (vf[j]!=0.0) {
      if (j!=k) {
        vf[k]    = vf[j];
        index[k] = index[j];
      }
      cnt++;
      k++;
    }
  }
  row->array1.cnt = cnt;
}

