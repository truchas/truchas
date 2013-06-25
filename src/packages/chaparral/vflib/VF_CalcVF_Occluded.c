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
@(#)    $RCSfile: VF_CalcVF_Occluded.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_Occluded.c,v $
@(#)
@(#)    DESCRIPTION:  Calculate the partially occluded view factor, FF(I,J)
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

double VF_CalcVF_Occluded(Poly *poly_i, Poly *poly_j, CandidateList *candidates,
                          int ncandidates, double visibility)
{
  char   *name[]={"MonteCarlo", "DoubleArea", "SimpleRatio", "Hemicube"};
  int    nseg_i=1, nseg_j=1, method=0;
  double vf=0.0, ff;
  /*.....
  double r, r_over_di, r_over_dj, area_i, area_j;
  Point3 center_i, center_j;
  .....*/
  VFenclosure *enclosure=VF_CurrentEnclosure();
  VFtopology  *topology=VF_CurrentTopology();

  if (topology->geom==VF_2Dplanar) {
    ff = VF_CalcVF_Hottel(poly_i, poly_j);
    if (ff>0.0) vf = visibility*ff;
  } else {
    /*.....
    area_i = VF_PolyArea(poly_i);
    area_j = VF_PolyArea(poly_j);
    VF_PolyCenter(poly_i, &center_i);
    VF_PolyCenter(poly_j, &center_j);
    r         = V3_DistanceBetween2Points(&center_i,&center_j);
    r_over_di = r/(2.0*sqrt(area_i/M_PI));
    r_over_dj = r/(2.0*sqrt(area_j/M_PI));
    if (r_over_di < 100.00 ) nseg_i = 2;
    if (r_over_di <   5.00 ) nseg_i = 3;
    if (r_over_di <   2.00 ) nseg_i = 4;
    if (r_over_di <   1.00 ) nseg_i = 5;
    if (r_over_di <   0.50 ) nseg_i = 6;
    if (r_over_di <   0.30 ) nseg_i = 7;
    if (r_over_di <   0.20 ) nseg_i = 8;
    if (r_over_di <   0.10 ) nseg_i = 9;
    if (r_over_di <   0.05 ) nseg_i = 10;
    if (r_over_dj < 100.00 ) nseg_j = 2;
    if (r_over_dj <   5.00 ) nseg_j = 3;
    if (r_over_dj <   2.00 ) nseg_j = 4;
    if (r_over_dj <   1.00 ) nseg_j = 5;
    if (r_over_dj <   0.50 ) nseg_j = 6;
    if (r_over_dj <   0.30 ) nseg_j = 7;
    if (r_over_dj <   0.20 ) nseg_j = 8;
    if (r_over_dj <   0.10 ) nseg_j = 9;
    if (r_over_dj <   0.05 ) nseg_j = 10;
    if (nseg_i>3 && nseg_j>3) {
      method = 1;
    }
    .....*/
    if (topology->debug_level>=VF_OUTPUT_DEBUG_2) {
      printf("           VF_CalcVF_Occluded, method = %s, ni= %d, nj = %d\n",
             name[method],nseg_i,nseg_j);
    }
    switch (method) {
    case 0:
      vf = VF_CalcVF_MonteCarlo(poly_i, poly_j, candidates, ncandidates);
      break;
    case 1:
      vf = VF_CalcVF_DoubleArea(poly_i, poly_j, candidates, ncandidates, 
                         nseg_i, nseg_j);
      break;
    case 2:
      vf = visibility*VF_CalcVF_Unoccluded(poly_i, poly_j);
      break;
    case 3:
      if (enclosure->adaptive.hc_mode) {
        vf = enclosure->adaptive.hemicube.vf[0];
      } else {
        vf = VF_CalcVF_Hemicube(poly_i, poly_j, candidates, ncandidates);
      }
      break;
    }
    if (topology->debug_level>=VF_OUTPUT_DEBUG_2) {
      printf("             vf = %g\n",vf);
    }
  }
  if (vf<0.0) vf = 0.0;
  return (vf);
}
