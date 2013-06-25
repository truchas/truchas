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
@(#)    $RCSfile: VF_CalcVF_HemicubeRow.c,v $
@(#)    $Revision: 1.4 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_HemicubeRow.c,v $
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

double VF_CalcVF_HemicubeRow (Facet *facet)
{
  int         k, iseg, jseg, sub_divide;
  int         facet_i, subfacet_i, nsub_polys_i;
  double      *vf, *VF, *VFptr;
  double      area_I, area_i, min_sep;
  double      eff_diameter, sub_area, new_diameter;
  double      min_dist, aspect, new_area;
  ViewPort    view;
  Poly        *poly_i;
  POLYstack   PolyStack_I;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();
  Hemicube    *hemicube=&(enclosure->hemicube);

  facet_i    = facet->index;
  area_I     = VF_FacetArea(facet);
  vf         = hemicube->vf;
  min_sep    = hemicube->min_separation;
  sub_divide = hemicube->sub_divide;
  VF_GetDPbuffer0_ptr(&VF);
  hemicube->min_equiv_radius = MIN(hemicube->min_equiv_radius, sqrt(area_I/M_PI));
  VF_FacetToPolyStack(facet, &PolyStack_I, &nsub_polys_i);
  for (area_I=0.0, subfacet_i=0; subfacet_i<nsub_polys_i; subfacet_i++) {
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1 && nsub_polys_i>1) {
      printf("       Sub poly_i %d of %d\n",subfacet_i,nsub_polys_i);
      fflush(stdout);
    }
    VF_PopPolyStackPtr(&PolyStack_I, poly_i);
    area_i	 = VF_PolyArea(poly_i);
    area_I	+= area_i;
    eff_diameter = 2.0*sqrt(area_i/M_PI);
    min_dist	 = VF_FindMinSeperationDist(poly_i, facet_i, subfacet_i);
    hemicube->min_distance = MIN(min_dist,hemicube->min_distance);
    if (min_dist<min_sep*eff_diameter) {
      if (poly_i->np == 4) {
        new_diameter = min_dist/min_sep;
        aspect       = VF_QuadAspectRatio(poly_i);
        if (aspect>=1.0) {
          sub_area = area_i/aspect;
          new_area = M_PI*0.25*new_diameter*new_diameter;
          jseg     = (int)(sqrt((double)(sub_area/new_area))+0.5);
          iseg     = (int)((double)jseg*aspect);
        } else {
          aspect   = 1.0/aspect;
          sub_area = area_i/aspect;
          new_area = M_PI*0.25*new_diameter*new_diameter;
          iseg     = (int)(sqrt((double)(sub_area/new_area))+0.5);
          jseg     = (int)((double)iseg*aspect);
        }
        iseg           = MAX(1,iseg);
        jseg           = MAX(1,jseg);
        hemicube->imax = MAX(hemicube->imax,iseg);
        hemicube->jmax = MAX(hemicube->jmax,jseg);
        if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0 && (iseg>1||jseg>1)) {
          printf("     Facet %d: want (%d,%d) subdivisions",facet_i,iseg,jseg);
          if (sub_divide<iseg || sub_divide<jseg) {
            printf(", clamping to (%d,%d)\n",
                   MIN(sub_divide,iseg),MIN(sub_divide,jseg));
          } else {
            printf("\n");
          }
        }
        iseg = MIN(sub_divide,iseg);
        jseg = MIN(sub_divide,jseg);
      } else {
        new_diameter             = min_dist/min_sep;
        iseg                     = (int)(eff_diameter/new_diameter+0.5);
        hemicube->imax = MAX(hemicube->imax,iseg);
        hemicube->jmax = MAX(hemicube->jmax,iseg);
        iseg                     = MIN(sub_divide,iseg);
        jseg                     = iseg;
      }
    } else {
      iseg = 1;
      jseg = 1;
    }
    VF_SetView(&view, poly_i);
    for (VFptr=vf, k=0; k<enclosure->npatches_g; k++) *VFptr++ = 0.0;
    if ((sub_divide>1) && (iseg>1 || jseg>1)) {
      VF_HemicubeSub(facet_i, subfacet_i, poly_i, 
                     &view, VF, vf, iseg, jseg);
    } else {
      VF_HemicubeProjectRow(facet_i, subfacet_i, poly_i,
                            &view, 1, area_i, VF, vf);
    }
  }
  return (area_I);
}
