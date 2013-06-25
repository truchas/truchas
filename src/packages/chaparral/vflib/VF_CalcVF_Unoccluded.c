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
@(#)    $RCSfile: VF_CalcVF_Unoccluded.c,v $
@(#)    $Revision: 1.4 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_Unoccluded.c,v $
@(#)
@(#)    DESCRIPTION:  Calculate the unoccluded view factor, FF(I,J)
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
 int isnan(double ff) {return(0);}
#endif

#include "vf.h"

double VF_CalcVF_Unoccluded(Poly *poly_i, Poly *poly_j)
{
  int        i, j, nseg_i, nseg_j, npoly_i, npoly_j;
  double     area_I, area_i, area_j;
  double     r, r_over_di, r_over_dj;
  double     FIJ, Fij, ff;
  Point      center_i, center_j;
  Poly       *poly_ii, *poly_jj;
  POLYstack  PolyStack_I, PolyStack_J;
  VFtopology *topology=VF_CurrentTopology();
        
  if (topology->geom==VF_2Dplanar) {
    ff = VF_CalcVF_Hottel(poly_i, poly_j);
    if (ff>0.0) FIJ = ff;
  } else {
    VF_InitPolyStack(&PolyStack_I);
    VF_InitPolyStack(&PolyStack_J);
    if (poly_i->np<=4) {
      npoly_i = 1;
      VF_PushPolyStack(&PolyStack_I, *poly_i);
    } else {
      npoly_i = 2;
      VF_PushSplitPolyStack(&PolyStack_I, poly_i);
    }
    if (poly_j->np<=4) {
      npoly_j = 1;
    } else {
      npoly_j = 2;
    }
    for (FIJ=0.0, area_I=0.0, i=0; i<npoly_i; i++) {
      VF_PopPolyStackPtr(&PolyStack_I, poly_ii);
      area_i  = VF_PolyArea(poly_ii);
      area_I += area_i;
      if (npoly_j==1) {
        VF_PushPolyStack(&PolyStack_J, *poly_j);
      } else {
        VF_PushSplitPolyStack(&PolyStack_J, poly_j);
      }
      for (Fij=0.0, j=0; j<npoly_j; j++) {
        VF_PopPolyStackPtr(&PolyStack_J, poly_jj);
        if (VF_SharedPolyEdge(poly_ii, poly_jj)) {
          ff = VF_CalcVF_Analytic(poly_ii, poly_jj);
          if (isnan(ff)) {
            ff = VF_CalcVF_Gauss(poly_ii, poly_jj, 10, 10);
            /*...
            printf("Shared poly edge\n");
            printf("1) VF_CalcVF_Analytic(poly_ii, poly_jj) = NaN\n");
            PrintPoly(poly_ii,"poly_ii");
            PrintPoly(poly_jj,"poly_jj");
            printf("   ff = %g\n",ff);
            ...*/
          }
        } else {
          nseg_i = 1;
          nseg_j = 1;
          area_j = VF_PolyArea(poly_jj);
          VF_PolyCenter(poly_ii, &center_i);
          VF_PolyCenter(poly_jj, &center_j);
          r         = V3_DistanceBetween2Points(&center_i,&center_j);
          r_over_di = r/(2.0*sqrt(area_i/M_PI));
          r_over_dj = r/(2.0*sqrt(area_j/M_PI));
          if (r_over_di < 100.0 ) nseg_i = 2;
          if (r_over_di <  10.0 ) nseg_i = 3;
          if (r_over_di <   5.0 ) nseg_i = 4;
          if (r_over_di <   2.0 ) nseg_i = 5;
          if (r_over_di <   1.0 ) nseg_i = 6;
          if (r_over_dj < 100.0 ) nseg_j = 2;
          if (r_over_dj <  10.0 ) nseg_j = 3;
          if (r_over_dj <   5.0 ) nseg_j = 4;
          if (r_over_dj <   2.0 ) nseg_j = 5;
          if (r_over_dj <   1.0 ) nseg_j = 6;
          if (nseg_i*nseg_j < 9) {
            ff = VF_CalcVF_Gauss(poly_ii, poly_jj, nseg_i, nseg_j);
          } else {
            ff = VF_CalcVF_Analytic(poly_ii, poly_jj);
            if (isnan(ff)) {
              ff = VF_CalcVF_Gauss(poly_ii, poly_jj, 10, 10);
              /*...
              printf("2) VF_CalcVF_Analytic(poly_ii, poly_jj) = NaN\n");
              PrintPoly(poly_ii,"poly_ii");
              PrintPoly(poly_jj,"poly_jj");
              printf("   ff = %g\n",ff);
              ...*/
            }
          }
        }
        if (ff>0.0) Fij += ff;
      }
      FIJ += area_i*Fij;
    }
    FIJ /= area_I;
  }
  return (FIJ);
}

