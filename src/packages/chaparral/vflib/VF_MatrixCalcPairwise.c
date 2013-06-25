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
@(#)    $RCSfile: VF_MatrixCalcPairwise.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_MatrixCalcPairwise.c,v $
@(#)
@(#)    DESCRIPTION:  Top-level routine to control the pairwise method
@(#)    to calculate the viewfactor matrix.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include "vf.h"

void VF_MatrixCalcPairwise(void) 
{
  int    i, j, n, next_gid;
  int    subfacet_i, subfacet_j, nsub_polys_i, nsub_polys_j;
  float  *VFptr, *vf;
  double rowsum, area_i, patch_area, parea;
  double t0, t1, t2, t, *Htime;
  Facet  *facet_i, *facet_j;
  Poly   *poly_i, *poly_j;
  POLYstack   PolyStack_I, PolyStack_J;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncStart();
    printf("   Processor %d ..........\n",VFLIB_Rank);
  }
  t0 = VF_Clock();
  enclosure->adaptive.num_total         = 0;
  enclosure->adaptive.num_zero          = 0;
  enclosure->adaptive.num_analytic      = 0;
  enclosure->adaptive.num_numerical     = 0;
  enclosure->adaptive.num_total_row     = 0;
  enclosure->adaptive.num_zero_row      = 0;
  enclosure->adaptive.num_analytic_row  = 0;
  enclosure->adaptive.num_numerical_row = 0;
  VF_InitPolyStack(&PolyStack_I);
  VF_InitPolyStack(&PolyStack_J);
  VF_GetSPbuffer0_ptr(&vf);
  
  /*
  if (!(topology->nonblocking)) {
    enclosure->adaptive.hemicube.resolution     = 500;
    enclosure->adaptive.hemicube.sub_divide     = 5;
    enclosure->adaptive.hemicube.min_separation = 5.0;
    VF_AllocHemicube(&(enclosure->adaptive.hemicube));
  }
  */

  area_i     = 0.0;
  patch_area = 0.0;
  for (VFptr=vf, j=0; j<enclosure->npatches_g; j++) *VFptr++ = 0.0;
  for (n=0, i=0; i<topology->nfacets_base;) {
    facet_i = &(topology->facets[i]);
    if (facet_i->patch_proc!=VFLIB_Rank) {
      i++; continue;
    }
    n++;
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
      printf("\n   Processing base facet %d of %d (patch gid %d)..............................\n",
             n,topology->nfacets_l,facet_i->patch_gid);
      fflush(stdout);
    }
    next_gid = -1;
    for (i++; i<topology->nfacets_base; i++) {
      if (topology->facets[i].patch_proc!=VFLIB_Rank) continue;
      next_gid = topology->facets[i].patch_gid;
      break;
    }
    patch_area += VF_FacetArea(facet_i);
    if (!(topology->nonblocking)) {
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1) {
        printf("     Initializing candidate list\n");
      }
      enclosure->adaptive.montecarlo.sampling.index = i+11;
      VF_SetupSampling(&(enclosure->adaptive.visibility.sampling),
                       enclosure->adaptive.visibility.sampleuv_i);
    }
    VF_InitCandidates(facet_i);
    VF_FacetToPolyStack(facet_i, &PolyStack_I, &nsub_polys_i);
    for (subfacet_i=0; subfacet_i<nsub_polys_i; subfacet_i++) {
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2 && nsub_polys_i>1) {
        printf("     Sub poly_i %d of %d\n",subfacet_i,nsub_polys_i);
        fflush(stdout);
      }
      VF_PopPolyStackPtr(&PolyStack_I, poly_i);
      area_i += VF_PolyArea(poly_i);
      for (j=0; j<topology->nfacets_g; j++) {
        facet_j = &(topology->facets[j]);
        if (facet_j->patch_global_index>facet_i->patch_global_index) continue;
        if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1) {
          printf("\n       Processing target facet %d of %d (patch gid %d)\n",
                 j,topology->nfacets_g,facet_j->patch_gid);
          fflush(stdout);
        }
        VF_FacetToPolyStack(facet_j, &PolyStack_J, &nsub_polys_j);
        for (subfacet_j=0; subfacet_j<nsub_polys_j; subfacet_j++) {
          if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2 && nsub_polys_j>1) {
            printf("         Sub poly_j %d of %d\n",subfacet_j,nsub_polys_j);
            fflush(stdout);
          }
          if (!(topology->nonblocking)) {
            enclosure->adaptive.montecarlo.sampling.index = j+11;
            enclosure->adaptive.visibility.sampling.index = j+11;
            VF_SetupSampling(&(enclosure->adaptive.visibility.sampling), 
                             enclosure->adaptive.visibility.sampleuv_j);
          }
          VF_PopPolyStackPtr(&PolyStack_J, poly_j);
          vf[facet_j->patch_global_index] += (float)(VF_PolyArea(poly_i)*VF_CalcVF_ComputePair(poly_i, poly_j));
        }
      }
    }
    if (facet_i->patch_gid!=next_gid) {
      parea = 1.0/area_i;
      for (rowsum=0.0, VFptr=vf, j=0; j<enclosure->npatches_g; j++) {
        if (*VFptr>0.0) {
          *VFptr *= parea;
          rowsum += *VFptr;
        }
        if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
          printf("       vf(%d, %d) = %g\n",facet_i->patch_global_index,j,*VFptr);
          fflush(stdout);
        }
        VFptr++;
      }
      /*==============================================*/
      /* STORE THE ROW OF VIEWFACTORS IN THE DATABASE */
      /*==============================================*/
      VF_LoadMatrixRow(facet_i->patch_local_index, vf, patch_area);
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
        printf("     Loading matrix row  %2d  lid = %2d  gid = %2d  hid =% 2d\n",
               facet_i->patch_local_index,
               enclosure->row[facet_i->patch_local_index].local_index,
               enclosure->row[facet_i->patch_local_index].global_index,
               enclosure->row[facet_i->patch_local_index].host_gid);
        printf("     N total       = %d\n",enclosure->adaptive.num_total_row);
        printf("       N zero      = %d\n",enclosure->adaptive.num_zero_row);
        printf("       N analytic  = %d\n",enclosure->adaptive.num_analytic_row);
        printf("       N numerical = %d\n",enclosure->adaptive.num_numerical_row);
        printf("     rowsum        = %g\n",rowsum);
      }
      area_i     = 0.0;
      patch_area = 0.0;
      for (VFptr=vf, j=0; j<enclosure->npatches_g; j++) *VFptr++ = 0.0;
      enclosure->adaptive.num_total_row     = 0;
      enclosure->adaptive.num_zero_row      = 0;
      enclosure->adaptive.num_analytic_row  = 0;
      enclosure->adaptive.num_numerical_row = 0;
    }
  }
  
  
  
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncEnd();
    VF_Sync();
  }
  if (enclosure->partial && (VFLIB_Rank==VFLIB_Size-1)) {
    for (VFptr=vf, j=0; j<enclosure->npatches_g; j++) *VFptr++ = 0.0;
    VF_LoadMatrixRow(enclosure->npatches_l-1,vf,enclosure->asink);
  }
  t1 = VF_Clock();
  t  = t0-t1;
  VF_Sync();
  /*=====================================================================*/
  /* FILL OUT THE UPPER DIAGONAL PORTION OF THE MATRIX USING RECIPROCITY */
  /*=====================================================================*/
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
      printf("   Filling out upper tiangle of matrix\n");
      fflush(stdout);
  }
  t1 = VF_Clock();
  VF_FillUpperDiagonal();
  t2 = VF_Clock();
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
    Htime = VF_Newd(VFLIB_Size);
    VF_AllgatherDouble(&t,Htime,1);
    if (VFLIB_Rank==0) {
      printf("     ===== HISTOGRAM INFO =====\n");
      for (i=0; i<VFLIB_Size; i++) {
        printf("     %d  %e\n",i,Htime[i]);
      }
      printf("     ==========================\n\n");
      printf("     Time to fill upper diagonal = %g\n",t2-t1);
      fflush(stdout);
    }
    VF_Free(Htime);
  }
  /*===============================================================*/
  /* IF A PARTIAL ENCLOSURE, CALCULATE THE LAST ROW OF VIEWFACTORS */
  /*===============================================================*/
  if (enclosure->partial) {
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_0) {
      printf("   Completing partial enclosure\n");
    }
    VF_PartialEnclosure();
  }
  if (enclosure->debug_level>=VF_OUTPUT_SUMMARY) {
    VF_Sync();
    VF_GlobalSumInt(&(enclosure->adaptive.num_total));
    VF_GlobalSumInt(&(enclosure->adaptive.num_zero));
    VF_GlobalSumInt(&(enclosure->adaptive.num_analytic));
    VF_GlobalSumInt(&(enclosure->adaptive.num_numerical));
  }
}
