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
@(#)    $RCSfile: VF_CalcVF_DoubleArea.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_DoubleArea.c,v $
@(#)
@(#)    DESCRIPTION:  Use double area method to numerically 
@(#)    evaluate the viewfactor between two surfaces.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vf.h"

double 
VF_CalcVF_DoubleArea (Poly *poly_i, Poly *poly_j, 
                      CandidateList *candidates, int ncandidates, 
                      int nseg_i, int nseg_j)
{
  int    i1, i2, j1, j2, n, nsi1, nsi2, nsj1, nsj2, vis;
  double tmax, ff, cf, tmp, aspect, vf=0.0, n_dot_d;
  Ray    ray;
  Poly   subpoly_i, subpoly_j;
  Point  center_i, center_j;
  Vector q;

  switch (poly_i->np) {
  case 2:
    nsi1 = nseg_i;
    nsi2 = 1;
    break;
  case 3:
    nsi1 = nseg_i;
    nsi2 = nseg_i;
    break;
  case 4:
    aspect = VF_QuadAspectRatio(poly_i);
    if (aspect>=1.0) {
      nsi1 = (int)((double)nseg_i*aspect);
      nsi2 = nseg_i;
    } else {
      nsi1 = nseg_i;
      nsi2 = (int)((double)nseg_i/aspect);
    }
    break;
  }
  switch (poly_j->np) {
  case 2:
    nsj1 = nseg_j;
    nsj2 = 1;
    break;
  case 3:
    nsj1 = nseg_j;
    nsj2 = nseg_j;
    break;
  case 4:
    aspect = VF_QuadAspectRatio(poly_j);
    if (aspect>=1.0) {
      nsj1 = (int)((double)nseg_j*aspect);
      nsj2 = nseg_j;
    } else {
      nsj1 = nseg_j;
      nsj2 = (int)((double)nseg_j/aspect);
    }
    break;
  }
  for (vf=0.0, n=0, i1=0; i1<nsi1; i1++) {
    if (poly_i->np==3) nsi2 = i1*2+1;
    for (i2=0; i2<nsi2; i2++) {
      VF_SubPoly(i1,i2,nsi1,nsi1,poly_i,&subpoly_i);
      VF_PolyCenter(&subpoly_i,&center_i);
      ray.O = center_i;
      for (ff=0.0, j1=0; j1<nsj1; j1++) {
        if (poly_j->np==3) nsj2 = j1*2+1;
        for (j2=0; j2<nsj2; j2++, n++) {
          VF_SubPoly(j1,j2,nsj1,nsj1,poly_j,&subpoly_j);
#define NO_HARDCORE
#ifdef HARDCORE
          cf = VF_Visibility(&subpoly_i, &subpoly_j, 
                             candidates, ncandidates)*
          VF_Unoccluded(&subpoly_i, &subpoly_j);
#else
          VF_PolyCenter(&subpoly_j,&center_j);
          V3_Sub(&center_j,&center_i,&(ray.D));
          V3_Normalize(&ray.D,tmp);
          n_dot_d = V3_Dot(&(poly_j->normal), &ray.D);
          if (fabs(n_dot_d)<VF_RAY_TEST) {
            vis = 1;
          } else {
            V3_Sub(&(subpoly_j.p[0]), &ray.O, &q);
            tmax = V3_Dot(&(poly_j->normal), &q)/n_dot_d;
            vis  = VF_RayPatchListTest(&ray, n, tmax, candidates, ncandidates);
          }
          if (vis) {
            cf = 0.0;
          } else {
             cf = VF_CalcVF_Unoccluded(&subpoly_i, &subpoly_j);
          }
#endif
          if (cf>0.0) ff += cf;
        }
      }
      if (ff>0.0) {
        vf += VF_PolyArea(&subpoly_i)*ff;
      }
    }
  }
  vf /= VF_PolyArea(poly_i);
  return (vf);
}
