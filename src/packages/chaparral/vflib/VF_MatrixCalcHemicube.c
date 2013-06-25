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
@(#)    $RCSfile: VF_MatrixCalcHemicube.c,v $
@(#)    $Revision: 1.7 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_MatrixCalcHemicube.c,v $
@(#)
@(#)    DESCRIPTION:  Use hemicube method to calculate the viewfactor matrix.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void VF_MatrixCalcHemicube()
{
  int         i, k, n;
  int         iseg, jseg, sub_divide;
  int         nfacets, next_gid;
  int         facet_i, subfacet_i, nsub_polys_i;
  float       *FF, *FFptr, rowsum;
  double      *VF, *vf, *VFptr;
  double      parea, area_I, area_i, patch_area, min_sep;
  double      eff_diameter, sub_area, new_diameter;
  double      min_dist, aspect, new_area;
  ViewPort    view;
  Poly        *poly_i;
  Facet       *facet;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();
  POLYstack   PolyStack_I;

  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncStart();
    printf("   Processor %d ..........\n",VFLIB_Rank);
  }
  VF_InitPolyStack(&PolyStack_I);
  nfacets    = topology->nfacets_base;
  sub_divide = enclosure->hemicube.sub_divide;
  min_sep    = enclosure->hemicube.min_separation;
  VF_GetSPbuffer0_ptr(&FF);
  VF_GetDPbuffer0_ptr(&VF);
  VF_GetDPbuffer1_ptr(&vf);
  VF_AllocHemicube(&(enclosure->hemicube));
  area_I     = 0.0;
  patch_area = 0.0;
  for (VFptr=VF, FFptr=FF, k=0; k<enclosure->npatches_g; k++) {
    *VFptr++ = 0.0;
    *FFptr++ = 0.0;
  }
  for (n=0, i=0; i<nfacets;) {
    facet = &(topology->facets[i]);
    if (facet->patch_proc!=VFLIB_Rank) {
      i++; continue;
    }
    n++;
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
      printf("\n   Processing facet %d of %d (patch gid %d)..............................\n",
             n,topology->nfacets_l,facet->patch_gid);
      fflush(stdout);
    }
    next_gid = -1;
    for (i++; i<nfacets; i++) {
      if (topology->facets[i].patch_proc!=VFLIB_Rank) continue;
      next_gid = topology->facets[i].patch_gid;
      break;
    }
    facet_i     = facet->index;
    patch_area += VF_FacetArea(facet);
    enclosure->hemicube.min_equiv_radius = MIN(enclosure->hemicube.min_equiv_radius, 
                                               sqrt(VF_FacetArea(facet)/M_PI));
    VF_FacetToPolyStack(facet, &PolyStack_I, &nsub_polys_i);
    for (subfacet_i=0; subfacet_i<nsub_polys_i; subfacet_i++) {
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1 && nsub_polys_i>1) {
        printf("       Sub poly_i %d of %d\n",subfacet_i,nsub_polys_i);
        fflush(stdout);
      }
      VF_PopPolyStackPtr(&PolyStack_I, poly_i);
      area_i       = VF_PolyArea(poly_i);
      area_I      += area_i;
      if (poly_i->np==2) {
        eff_diameter = area_i;
      } else {
        eff_diameter = 2.0*sqrt(area_i/M_PI);
      }
      min_dist = VF_FindMinSeperationDist(poly_i, facet_i, subfacet_i);
      enclosure->hemicube.min_distance = MIN(min_dist,
                                             enclosure->hemicube.min_distance);
      if (min_dist<min_sep*eff_diameter) {
        switch (poly_i->np) {
        case 2:
          new_diameter             = min_dist/min_sep;
          iseg                     = (int)(eff_diameter/new_diameter+0.5);
          enclosure->hemicube.imax = MAX(enclosure->hemicube.imax,iseg);
          iseg                     = MIN(sub_divide,iseg);
          jseg                     = 1;
          break;
        case 3:
          new_diameter             = min_dist/min_sep;
          iseg                     = (int)(eff_diameter/new_diameter+0.5);
          enclosure->hemicube.imax = MAX(enclosure->hemicube.imax,iseg);
          enclosure->hemicube.jmax = MAX(enclosure->hemicube.jmax,iseg);
          iseg                     = MIN(sub_divide,iseg);
          jseg                     = iseg;
          break;
        case 4:
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
          iseg                     = MAX(1,iseg);
          jseg                     = MAX(1,jseg);
          enclosure->hemicube.imax = MAX(enclosure->hemicube.imax,iseg);
          enclosure->hemicube.jmax = MAX(enclosure->hemicube.jmax,jseg);
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
          break;
        }
      } else {
        iseg = 1;
        jseg = 1;
      }
      VF_SetView(&view, poly_i);
      if ((sub_divide>1) && (iseg>1 || jseg>1)) {
        VF_HemicubeSub(facet_i, subfacet_i, poly_i, 
                       &view, VF, vf, iseg, jseg);
      } else {
        VF_HemicubeProjectRow(facet_i, subfacet_i, poly_i, 
                              &view, 1, area_i, VF, vf);
      }
    }
    if (facet->patch_gid!=next_gid) {
      parea = 1.0/area_I;
      for (VFptr=VF, FFptr=FF, k=0; k<enclosure->npatches_g; k++) {
        if (*VFptr>0.0) {
          *FFptr = parea*(*VFptr);
        }
        VFptr++;
        FFptr++;
      }
      /*=======================================================*/
      /* IF A PARTIAL ENCLOSURE, CALCULATE THE LAST VIEWFACTOR */
      /*=======================================================*/
      if (enclosure->partial) {
        for (rowsum=0.0, FFptr=FF, k=0; k<enclosure->npatches_g; k++) {
          rowsum += *FFptr++;
        }
        FF[enclosure->npatches_g-1] = MAX((float)0.0,(float)1.0-rowsum);
      }
      /*==============================================*/
      /* STORE THE ROW OF VIEWFACTORS IN THE DATABASE */
      /*==============================================*/
      VF_LoadMatrixRow(facet->patch_local_index, FF, patch_area);
      VF_SetRawRowsum(facet->patch_local_index);
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
        printf("     Loading matrix row  %2d  lid = %2d  gid = %2d  hid =% 2d\n",
               facet->patch_local_index,
               enclosure->row[facet->patch_local_index].local_index,
               enclosure->row[facet->patch_local_index].global_index,
               enclosure->row[facet->patch_local_index].host_gid);
        if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1) {
          for (k=0; k<enclosure->npatches_g; k++) {
            printf("     F(%d,%d) = %.10g\n",facet->patch_global_index,k,FF[k]);
          }
        }
        printf("     rowsum = %g (%g)\n",
               enclosure->row[facet->patch_local_index].raw_rowsum,
               enclosure->row[facet->patch_local_index].raw_rowsum-1.0);
        printf("     area   = %g\n",patch_area);
        fflush(stdout);
      }
      area_I     = 0.0;
      patch_area = 0.0;
      for (VFptr=VF, FFptr=FF, k=0; k<enclosure->npatches_g; k++) {
        *VFptr++ = 0.0;
        *FFptr++ = 0.0;
      }
    }
  }
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncEnd();
    VF_Sync();
  }
  VF_FreeHemicube(&(enclosure->hemicube));
  /*======================================*/
  /* WARN THE USER IF HOLES WERE DETECTED */
  /*======================================*/
  VF_MaxInt(&(enclosure->hemicube.nholes),0);
  if (VFLIB_Rank==0 && !enclosure->partial && enclosure->hemicube.nholes>0) {
    fflush(stdout);
    printf("\n     WARNING: ");
    if (enclosure->hemicube.nholes==1) {
      printf("holes were detected during main pass\n");
    } else if (enclosure->hemicube.nholes==2) {
      printf("holes were detected during subdivision pass\n");
    } else if (enclosure->hemicube.nholes==3) {
      printf("holes were detected during main and subdivision pass\n");
    } else {
      printf("holes were detected\n");
    }
    printf("\n");
    printf("         - this can be due to an ill-defined enclosure\n");
    printf("           or may just be small errors in the polygon\n");
    printf("           clipping and scan conversion routines\n");
    printf("\n");
  }
  /*===============================================================*/
  /* IF A PARTIAL ENCLOSURE, CALCULATE THE LAST ROW OF VIEWFACTORS */
  /*===============================================================*/
  if (enclosure->partial) {
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
      printf("   Completing partial enclosure\n");
    }
    if (VFLIB_Rank==VFLIB_Size-1) {
      for (FFptr=FF, k=0; k<enclosure->npatches_g; k++) *FFptr++ = 0.0;
      VF_LoadMatrixRow(enclosure->npatches_l-1,FF,enclosure->asink);
    }
    VF_PartialEnclosure();
  }
}
