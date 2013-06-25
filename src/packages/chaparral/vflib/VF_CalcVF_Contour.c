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
@(#)    $RCSfile: VF_CalcVF_Contour.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_Contour.c,v $
@(#)
@(#)    DESCRIPTION:  Use line contour method to calculate the viewfactor
@(#)    bewtween 2 unoccluded polygons.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

double 
VF_CalcVF_Contour(Poly *poly_i, Poly *poly_j, int nseg_i, int nseg_j)
{
  int    i, ii, j, jj, ni, nj, np_i, np_j;
  double scale_i, scale_j, r, vdot, vf=0.0;
  Point  sample_i, sample_j, v_i, v_j;

  np_i = poly_i->np;
  np_j = poly_j->np;
  VF_PolyCenter(poly_i, &sample_i);
  VF_PolyCenter(poly_j, &sample_j);
  scale_i = 1.0/(double)(nseg_i);
  scale_j = 1.0/(double)(nseg_j);
  for (vf=0.0, ii=1, i=0; i<np_i; i++, ii=(i+1)%np_i) {
    V3_Sub(&(poly_i->p[ii]),&(poly_i->p[i]),&v_i);
    V3_Scale(&v_i,scale_i);
    V3_Init(&v_i,&sample_i);
    V3_Scale(&sample_i,0.5);
    V3_Add(&sample_i,&(poly_i->p[i]),&sample_i);
    for (ni=0; ni<nseg_i; ni++) {
      for (jj=1, j=0; j<np_j; j++, jj=(j+1)%np_j) {
        V3_Sub(&(poly_j->p[jj]),&(poly_j->p[j]),&v_j);
        V3_Scale(&v_j,scale_j);
        V3_Init(&v_j,&sample_j);
        V3_Scale(&sample_j,0.5);
        V3_Add(&sample_j,&(poly_j->p[j]),&sample_j);
        vdot = V3_Dot(&v_i,&v_j);
        for (nj=0; nj<nseg_j; nj++) {
          r   = V3_DistanceBetween2Points(&sample_i,&sample_j);
          vf += log(r) * vdot;
          V3_Add(&sample_j, &v_j, &sample_j);
        }
      }
      V3_Add(&sample_i, &v_i, &sample_i);
    }
  }
  vf /= 2.0*M_PI*VF_PolyArea(poly_i);
  return (vf);
}
