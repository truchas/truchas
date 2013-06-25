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
@(#)    $RCSfile: VF_CalcVF_MonteCarlo.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_MonteCarlo.c,v $
@(#)
@(#)    DESCRIPTION:  Use monte-carlo method to numerically 
@(#)    evaluate the viewfactor between two surfaces.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vf.h"

double VF_CalcVF_MonteCarlo (Poly *poly_i, Poly *poly_j, 
                             CandidateList *candidates, 
                             int ncandidates)
{
  int     i, n, ns, nv, nw=1000, vis, check=0;
  double  rsq, cos_theta_i, cos_theta_j, n_dot_d;
  double  tmax, area_j, delta_ff, ff, vf=0.0, Nt, tmp;
  double  w_sum1=0.0, w_sum2=0.0, w_mean, w_stddev, Nw, viz;
  Vector  q, r;
  Ray     ray;
  static  Point  *sample_i=NULL, *sample_j=NULL;
  static  double *wff1=NULL, *wff2=NULL;
  static  int    max_samples=0;
  Adaptive *adaptive;
  VFenclosure *enclosure=VF_CurrentEnclosure();
    
  adaptive = &(enclosure->adaptive);
  if (wff1==NULL) wff1 = VF_Newd(nw);
  if (wff2==NULL) wff2 = VF_Newd(nw);
  Nw = (double)nw;
  ns = adaptive->montecarlo.sampling.n;         
  if (max_samples=0 || sample_i==NULL || sample_j==NULL) {
    max_samples = ns;
    sample_i    = (Point *)VF_Newv(ns*sizeof(Point));
    sample_j    = (Point *)VF_Newv(ns*sizeof(Point));
  } else if (max_samples<ns) {
    max_samples = ns;
    sample_i    = (Point *)VF_ReNewv(sample_i,ns*sizeof(Point));
    sample_j    = (Point *)VF_ReNewv(sample_j,ns*sizeof(Point));
  }
  if (adaptive->montecarlo.sampling.method==VF_RANDOM_SAMPLE ||
      adaptive->montecarlo.sampling.method==VF_HALTON_SAMPLE) check=1;
  VF_SamplesUVtoXYZ(poly_i, &(adaptive->montecarlo.sampling), 
                    adaptive->montecarlo.uv_samples, sample_i);
  VF_SamplesUVtoXYZ(poly_j, &(adaptive->montecarlo.sampling), 
                    adaptive->montecarlo.uv_samples, sample_j);
  area_j = VF_PolyArea(poly_j);
  for (i=0, nv=0, n=0; n<ns; n++) {
    V3_Sub(&sample_j[n],&sample_i[n],&r);
    ray.O = sample_i[n];
    ray.D = r;
    V3_Normalize(&ray.D,tmp);
    n_dot_d = V3_Dot(&(poly_j->normal), &ray.D);
    if (fabs(n_dot_d)<VF_RAY_TEST) {
      vis = 1;
    } else {
      V3_Sub(&(poly_j->p[0]), &ray.O, &q);
      tmax = V3_Dot(&(poly_j->normal), &q)/n_dot_d;
      vis  = VF_RayPatchListTest(&ray, n, tmax, candidates, ncandidates);
    }
    if (!vis) {
      nv++;
      rsq         = V3_Dot(&r,&r);
      cos_theta_i = V3_Dot(&ray.D,&(poly_i->normal));
      cos_theta_j = V3_Dot(&ray.D,&(poly_j->normal));
      delta_ff    = -(cos_theta_i*cos_theta_j)/rsq;
      if (delta_ff<0.0) {
        delta_ff = 0.0;
      }
    } else {
      delta_ff = 0.0;
    }
    vf += delta_ff;
    Nt  = (double)(n+1);
    if (check) {
      ff = area_j*vf/(Nt*M_PI);
      if (n>=nw) {
        w_sum1 += (ff-wff1[i]);
        w_sum2 += (ff*ff-wff2[i]);
        wff1[i] = ff;
        wff2[i] = ff*ff;
        if (nv>10 && ff>0.0) {
          w_stddev = sqrt((Nw*w_sum2-w_sum1*w_sum1)/(Nw*(Nw-1.0)));
          if (w_stddev/ff<adaptive->montecarlo.tol1 ||
              w_stddev   <adaptive->montecarlo.tol2) break;
        }
      } else {
        w_sum1 += ff;
        w_sum2 += ff*ff;
        wff1[i] = ff;
        wff2[i] = ff*ff;
      }
      i++; if (i==nw) i=0;
      if (n>=10000 && nv==0) break;
    }
  }
  if (nv==n+1) {
    ff = VF_CalcVF_Unoccluded(poly_i, poly_j);
  } else {
    ff = area_j*vf/(Nt*M_PI);
  }
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1) {
    viz    = (double)(100*nv)/Nt;
    w_mean = w_sum1/Nw;
    printf("     %7d samples:   ff = %.6e\n",n,ff);
    if (nv>0) {
      printf("                       wff = %.6e\n",w_mean);
      printf("                       wsd = %.6e\n",w_stddev);
      printf("                     wres1 = %.6e\n",w_stddev/ff);
      printf("                     wres2 = %.6e\n",w_stddev/w_mean);
    }
    printf("                       vis = %6.2f%%  %d/%d\n",viz,nv,n+1);
  }
  return (ff);
}
