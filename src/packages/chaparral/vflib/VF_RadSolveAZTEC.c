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
@(#)    $RCSfile: VF_RadSolveAZTEC.c,v $
@(#)    $Revision: 1.4 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_RadSolveAZTEC.c,v $
@(#)
@(#)    DESCRIPTION:  AZTEC solver interface to solve the radiosity matrix 
@(#)    equation.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

#ifdef AZTEC

# include <az_aztec.h>

void sort_update_array(int*, int*, int);
void cheap_find_local_indices(int, int*, int*, int**, int*);
void cheap_transform(int*, int**, int*, double*, int*, int**, int**, int**, int);
void VF_ReorderAztecMatrix(int, int*, double*, int*, int*, int);

static int    *update=NULL;
static int    *update_index=NULL;
static int    *bindx=NULL;
static double *sol=NULL;
static double *rhs=NULL;
static double *val=NULL;
  
void VF_RadsolveCleanupAZTEC(void)
{
  VF_Free(update);
  VF_Free(update_index);
  VF_Free(bindx);
  VF_Free(sol);
  VF_Free(rhs);
  VF_Free(val);
}

void
VF_RadSolveAZTEC(double *radq, double *radj, double *tsurf, double *eps, 
                 double sigma, double tol, int max_iter, int *num_iter, 
		 int dynamic, int solver, int debug)
{
  int    i, ii, j, jj, n_l, n_g, n, *index, *b_tmp, az_solver;
  int    new_size, MSRnonzeros, N_update, iout, ierr, n_tmp;
  int    options[AZ_OPTIONS_SIZE], proc_config[AZ_PROC_SIZE];
  float  *vf, vf_ii, *v_tmp;
  double rho, t4, qsum, L2_norm, *gsum;
  double start, end, elapsed=0.0, time1=0.0, time2=0.0, time3=0.0;
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];
  double tmin, tmax, emin, emax;
  static int    dat_size=0, vec_size=0, mat_size=0;
  VFenclosure   *enclosure=VF_CurrentEnclosure();

  switch (solver) {
  case VF_RADSOLVE_AZTEC_CG:
    az_solver = AZ_cg;
    break;
  case VF_RADSOLVE_AZTEC_GMRES:
    az_solver = AZ_gmres;
    break;
  }
  if (debug) {
    start = VF_Clock();
    if (VFLIB_Rank==0) {
      printf("\n\n");
      printf("  **********************************");
      printf("**********************************\n");
      printf("   R A D I O S I T Y    M A T R I X    S O L V E R\n");
      printf("  **********************************");
      printf("**********************************\n");
      printf("\n");
      printf("   Solving radiosity equations with AZTEC_");
      switch (solver) {
      case VF_RADSOLVE_AZTEC_CG:
        printf("CG");
        az_solver = AZ_cg;
        break;
      case VF_RADSOLVE_AZTEC_GMRES:
        printf("GMRES");
        az_solver = AZ_gmres;
        break;
      default:
        printf("unknown");
        break;
      }
      printf(" for enclosure <%s>\n",enclosure->id);
      fflush(stdout);
    }
  }
#ifndef VF_NO_MPI
  AZ_set_proc_config(proc_config,VFLIB_Comm);
#else
  AZ_set_proc_config(proc_config,AZ_NOT_MPI);
#endif
  AZ_defaults(options, params);
  n_g      = enclosure->npatches_g;
  n_l      = enclosure->npatches_l;
  N_update = n_l;
  iout     = MAX(0,debug-1);  /* -3 for everything */
  VF_GetDPbuffer0_ptr(&gsum);
  VF_GetSPbuffer1_ptr(&v_tmp);
  VF_GetINTbuffer_ptr(&b_tmp);

  /*========================================*/
  /* COUNT NON-ZEROES FOR MSR MATRIX SIZING */
  /*========================================*/
  for (MSRnonzeros=0, i=0; i<n_l; i++) {
    MSRnonzeros += enclosure->row[i].array0.cnt;
    MSRnonzeros += enclosure->row[i].array1.cnt;
    if (enclosure->row[i].diagonal != 0.0) MSRnonzeros--;
  }
  /*=======================================*/
  /* ALLOCATE MEMORY FOR THE UPDATE VECTOR */
  /*=======================================*/
  new_size = N_update;
  if (dat_size==0) {
    dat_size     = new_size;
    update       = VF_Newi(dat_size);
    update_index = VF_Newi(dat_size);
  } else if (new_size>dat_size) {
    dat_size     = new_size;
    update       = VF_ReNewi(update,      dat_size);
    update_index = VF_ReNewi(update_index,dat_size);
  }
  /*====================================*/
  /* ALLOCATE MEMORY FOR THE MSR MATRIX */
  /*====================================*/
  new_size = MSRnonzeros+N_update+10;
  if (mat_size==0) {
    mat_size = new_size;
    bindx    = VF_Newi(mat_size);
    val      = VF_Newd(mat_size);
  } else if (new_size>mat_size) {
    mat_size = new_size;
    bindx    = VF_ReNewi(bindx,mat_size);
    val      = VF_ReNewd(val,  mat_size);
  }
  /*=============================*/
  /* CONSTRUCT THE UPDATE VECTOR */
  /*=============================*/
  for (i=0; i<N_update; i++) {
    update[i] = enclosure->row[i].global_index;
  }
  /*==========================*/
  /* CONSTRUCT THE MSR MATRIX */
  /*==========================*/
  bindx[0] = N_update+1;
  for (i=0; i<N_update; i++) {
    ii     = update[i];
    n      = bindx[i];
    vf_ii  = enclosure->row[i].diagonal;
    rho    = 1.0-eps[ii];
    val[i] = 1.0-rho*vf_ii;
    n_tmp  = 0;
    vf     = enclosure->row[i].array0.data;
    index  = enclosure->row[i].array0.index;
    for (j=0; j<enclosure->row[i].array0.cnt; j++) {
      if (index[j]!=ii) {
        v_tmp[n_tmp] = -rho*vf[j];
        b_tmp[n_tmp] = index[j];
        n_tmp++;
      }
    }
    vf    = enclosure->row[i].array1.data;
    index = enclosure->row[i].array1.index;
    for (j=0; j<enclosure->row[i].array1.cnt; j++) {
      if (index[j]!=ii) {
        v_tmp[n_tmp] = -rho*vf[j];
        b_tmp[n_tmp] = index[j];
        n_tmp++;
      }
    }
    VF_SortSparseArrayAux(b_tmp,v_tmp,n_tmp);
    for (j=0; j<n_tmp; j++) {
      val[n]   = v_tmp[j];
      bindx[n] = b_tmp[j];
      n++;
    }
    bindx[i+1] = n;
  }
  /*=======================================*/
  /* REORDER THE MATRIX BASED ON INTERNAL  */
  /* AND EXTERNAL NODES AND THEN TRANSFORM */
  /* THE GLOBALLY INDEXED MATRIX TO A      */
  /* LOCALLY INDEXED MATRIX AND REORDER    */
  /*=======================================*/
  if (debug) {
    time2 = VF_Clock();
  }
  if (enclosure->data_org==NULL) {
    AZ_transform(proc_config, &(enclosure->external), bindx,
                 val, update, &(enclosure->update_index), 
                 &(enclosure->external_index), 
                 &(enclosure->data_org), N_update,
                 NULL, NULL, NULL, NULL, AZ_MSR_MATRIX);
  } else {
    cheap_transform(proc_config, &(enclosure->external), bindx, 
                    val, update, &(enclosure->update_index), 
                    &(enclosure->external_index), 
                    &(enclosure->data_org), N_update);
  }
  /*==================================*/
  /* SET AZTEC OPTIONS AND PARAMETERS */
  /*==================================*/
  options[AZ_solver]   = az_solver;
  options[AZ_orthog]   = AZ_modified;  /*AZ_classic;*/
  options[AZ_scaling]  = AZ_none;      /*AZ_sym_diag;*/
  options[AZ_precond]  = AZ_none;      /*AZ_Jacobi;*/
  options[AZ_poly_ord] = 2;
  options[AZ_conv]     = AZ_noscaled;  /*AZ_r0;*/
  options[AZ_output]   = iout;
  options[AZ_max_iter] = max_iter;
  options[AZ_kspace]   = 20; 
  params[AZ_tol]       = tol;
  if (debug) {
    time2 = VF_Clock()-time2;
  }
  /*====================================================*/
  /* ALLOCATE MEMORY FOR THE SOLN VECTOR AND RHS VECTOR */
  /*====================================================*/
  new_size = N_update+enclosure->data_org[AZ_N_external];
  if (vec_size==0) {
    vec_size = new_size;
    rhs      = VF_Newd(vec_size);
    sol      = VF_Newd(vec_size);
  } else if (new_size>vec_size) {
    vec_size = new_size;
    rhs      = VF_ReNewd(rhs,vec_size);
    sol      = VF_ReNewd(sol,vec_size);
  }
  /*============================================*/
  /* CONSTRUCT INITIAL SOLUTION AND RHS VECTORS */
  /*============================================*/
  for (i=0; i<N_update; i++) {
    ii = update[i];
    jj = enclosure->update_index[i];
    if (i==0) {
      tmin = tsurf[ii];
      tmax = tsurf[ii];
      emin = eps[ii];
      emax = eps[ii];
    } else {
      tmin = MIN(tmin,tsurf[ii]);
      tmax = MAX(tmax,tsurf[ii]);
      emin = MIN(emin,eps[ii]);
      emax = MAX(emax,eps[ii]);
    }
    t4 = sigma*tsurf[ii]*tsurf[ii]*tsurf[ii]*tsurf[ii];
    if (eps[ii]==0.0) {
      sol[jj] = t4;
      rhs[jj] = 0.0;
    } else {
      sol[jj] = t4+radq[ii]*(1.0-eps[ii])/eps[ii];
      rhs[jj] = eps[ii]*t4;
    }
  }
  VF_MinDouble(&tmin, 0);
  VF_MinDouble(&tmax, 0);
  VF_MinDouble(&emin, 0);
  VF_MinDouble(&emax, 0);
  if (debug) {
    time3 = VF_Clock();
    time1 = (time3-start)-time2;
  }
  /*=====================================*/
  /* SOLVE THE RADIOSITY MATRIX EQUATION */
  /*=====================================*/
  AZ_solve(sol,rhs,options,params,NULL,bindx,NULL,NULL,NULL,
           val,enclosure->data_org,status,proc_config);
  L2_norm   = status[AZ_rec_r];
  *num_iter = status[AZ_its];
  /*=========================================================*/
  /* FINISHED WITH RADIOSITY MATRIX SOLVE, CHECK CONVERGENCE */
  /*=========================================================*/
  if (*num_iter==max_iter) {
    if (VFLIB_Rank==0) {
      printf("\n");
      printf("   *** WARNING ***\n");
      printf("   NO CONVERGENCE FOR RADIATION PROBLEM\n");
      printf("   IN ENCLOSURE %s AFTER %3d ITERATIONS\n\n",
             enclosure->id,max_iter);
    }
  } else if (*num_iter<0) { 
    if (VFLIB_Rank==0) {
      printf("\n");
      printf("   *** WARNING ***\n");
      printf("   DIVERGING SOLUTION FOR RADIATION PROBLEM\n");
      printf("   IN ENCLOSURE %s AFTER %3d ITERATIONS\n\n",
             enclosure->id,-(*num_iter));
    }
  } else {
    /*=============================*/
    /* CONVERT TO GLOBAL NUMBERING */
    /*=============================*/
    for (i=0; i<N_update; i++) {
      gsum[update[i]] = sol[enclosure->update_index[i]];
    }
    /*====================================*/
    /* RADIOSITY MATRIX SOLVE SUCCESSFUL, */
    /* COMPUTE THE FINAL HEAT FLUXES      */
    /*====================================*/
    VF_ExchangeDouble(gsum);
    VF_ComputeFluxes(radq, radj, gsum, &qsum, debug);
  }
  /*====================================================*/
  /* FREE UP SOME MEMORY IF THIS IS A DYNAMIC ENCLOSURE */
  /*====================================================*/
  if (dynamic) {
    VF_Free(enclosure->data_org);
    VF_Free(enclosure->update_index);
    VF_Free(enclosure->external);
    VF_Free(enclosure->external_index);
  }
  /*==========================================*/
  /* PRINT SOME STATISTICS IF IN VERBOSE MODE */
  /*==========================================*/
  if (debug) {
    end     = VF_Clock();
    time3   = end-time3;
    elapsed = end-start;
    if (VFLIB_Rank==0) {
      printf("\n");
      printf("   Dynamic flag:        %d\n",dynamic);
      printf("   Total iterations:    %d\n",*num_iter);
      printf("   Total elapsed time:  %.2f (sec.)\n",elapsed);
      printf("     AZTEC setup time:  %.2f (sec.)\n",time1);
      printf("     AZTEC xform time:  %.2f (sec.)\n",time2);
      printf("     AZTEC solve time:  %.2f (sec.)\n",time3);
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

void sort_update_array(int array1[], int array2[], int nelements)
{
  int i,j,k,l,ir,rra,rri;

  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        k   = --l-1;
        rra = array1[k];
        rri = array2[k];
      } else {
        rra          = array1[ir-1];
        rri          = array2[ir-1];
        array1[ir-1] = array1[0];
        array2[ir-1] = array2[0];
        if (--ir == 1) {
          array1[0] = rra;
          array2[0] = rri;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (array1[j-1]<array1[j])) ++j;
        if (rra<array1[j-1]) {
          array1[i-1] = array1[j-1];
          array2[i-1] = array2[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      array1[i-1] = rra;
      array2[i-1] = rri;
    }
  }
}

void cheap_find_local_indices(int N_update, int bindx[], int update[],
                              int **external, int *N_external)
{
  int  j, k;
  int *bins,shift;
  int  start,end;

  /************************** execution begins ****************************/

  /* set up some bins so that we will be able to use AZ_quick_find() */

  bins = VF_Newi(N_update / 4 + 10);

  AZ_init_quick_find(update, N_update, &shift, bins);

  /*
   *  Compute the location of the first and last
   *  column index that is stored in the bindx[].
   */

  start = bindx[0]; 
  end   = bindx[N_update];
  for (j=start; j<end; j++) {
    k = AZ_quick_find(bindx[j], update, N_update,shift,bins);
    if (k != -1) {
      bindx[j] = k;
    } else {
      k        = AZ_find_index(bindx[j],*external,*N_external);
      bindx[j] = k + N_update;
    }
  }
  VF_Free(bins);
}

void cheap_transform(int proc_config[], int *external[], int bindx[],
                     double val[], int update[], int *update_index[],
                     int *extern_index[], int *data_org[], int N_update)
{
  int         N_extern;   /* Number of pts needed by this processor for
                             matrix-vector multiply but not updated by this
                             processor.  */

  N_extern = (*data_org)[AZ_N_external];
  
  cheap_find_local_indices(N_update, bindx, update, external, &N_extern);
  /*
  MWG_reorder_matrix(N_update, bindx, val, *update_index,
                     *extern_index,N_extern);
  */
  AZ_reorder_matrix(N_update, bindx, val, *update_index, *extern_index,
                    NULL, NULL, NULL, N_extern, NULL, AZ_ALL, AZ_MSR_MATRIX);
}

void VF_ReorderAztecMatrix(int N_update, int bindx[], double val[],
                           int update_index[], int extern_index[],int N_external)

/*******************************************************************************

  Reorder the matrix so that it corresponds to the new ordering given by
  'update_index' and 'extern_index'.

  IMPORTANT: This routine assumes that update_index[] contains two sequencies
  of numbers that are ordered but intertwined. For example,

  update_index:  4 5 0 6 1 2 3 7

  seq 1 =>    0   1 2 3

  seq 2 =>4 5   6 7

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  val,
  bindx:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  update_index:    update_index[i] gives the local numbering of global point
                   'update[i]'.

  extern_index:    extern_index[i] gives the local numbering of global point
                   'external[i]'.

  N_external:      Number of external elements on this processor.

*******************************************************************************/
{
  int   start, end;
  int   val_length, indx_length;
  int  *temp;
  int   i, j;
  char *yo = "MWG_reorder_matrix: ";

  /**************************** execution begins ******************************/

  start = N_update+1;      /* first nonzero offdiag */
  end   = bindx[N_update]; /* last nonzero          */
  
  /*
   * Change column indices (bindx) to reflect new ordering depending depending
   * on whether or not a point is internal or external.
   */

  for (i = start; i < end; i++) {
    if (bindx[i] < N_update) bindx[i] = update_index[bindx[i]];
    else                     bindx[i] = extern_index[bindx[i] - N_update];
  }

  /* reorder rows */


  /* We move the rows in four steps:
   *  1) sort the first N_update values of 'val'.
   *  2) sort the first N_update values of 'bindx'.
   *     We do this by first converting the ptrs to values
   *     representing the number of nonzero off diagonals.
   *  3) sort the off diagonal column indices.
   *  4) sort the off diagonal matrix nozeros.
   */

  j = bindx[0];
  AZ_convert_ptrs_to_values(bindx, N_update);
  AZ_sortqlists((char *) &(bindx[N_update + 1]), bindx, update_index,
                end - N_update - 1, sizeof(int), N_update);
  AZ_sortqlists((char *) &(val[N_update + 1]), bindx, update_index,
                end - N_update - 1, sizeof(double), N_update);
  AZ_sortqlists((char *) val, 0, update_index, N_update, sizeof(double),
                N_update);
  AZ_sortqlists((char *) bindx, 0, update_index, N_update,
                sizeof(int), N_update);
  AZ_convert_values_to_ptrs(bindx, N_update, j);   
}

#else

void
VF_RadSolveAZTEC(double *radq, double *radj, double *tsurf, double *eps, double sigma,
                 double tol, int max_iter, int *num_iter, int dynamic, 
                 int solver, int debug)
{
  if (VFLIB_Rank==0) {
    printf("FATAL ERROR:  This version of CHAPARRAL has not been built\n");
    printf("              with support for the AZTEC solvers\n");
    VF_Exit(-1);
  }
}

void
VF_RadSolveCleanupAZTEC(void)
{}

#endif
