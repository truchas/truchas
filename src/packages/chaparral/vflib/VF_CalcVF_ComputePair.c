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
@(#)    $RCSfile: VF_CalcVF_ComputePair.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_ComputePair.c,v $
@(#)
@(#)    DESCRIPTION:  Top-level routine to control the adaptive method
@(#)    to calculate the viewfactor matrix.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>

#include "vf.h"

double VF_CalcVF_ComputePair(Poly *poly_i, Poly *poly_j) 
{
  char          *indent="           ";
  int           i, j, n, ncandidates;
  int           npoly_i, npoly_j;
  double        ff, area;
  double        visibility;
  Poly          poly_ii, poly_jj, poly_I, poly_J;
  CandidateList *candidates;
  POLYstack     PolyStack_I, PolyStack_J;
  VFenclosure   *enclosure=VF_CurrentEnclosure();
  
  enclosure->adaptive.hc_mode = 0;
  VF_InitPolyStack(&PolyStack_I);
  VF_InitPolyStack(&PolyStack_J);
  /*========================================*/
  /* CHECK IF POLYGON J IS BEHIND POLYGON I */
  /*========================================*/
  if (VF_BehindPoly(poly_i, poly_j) || VF_BehindPoly(poly_j, poly_i)) {
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
      printf("%spoly j is behind poly i\n",indent);
    }
    ff = 0.0;
    enclosure->adaptive.num_zero++;
    enclosure->adaptive.num_zero_row++;
  } else {
    VF_ClipToPolyPlane(poly_i, &poly_ii, poly_j);
    VF_ClipToPolyPlane(poly_j, &poly_jj, poly_i);
    if (VF_BackFaceCullPolys(&poly_ii, &poly_jj)) {
      /*==================================================*/
      /* CHECK IF POLYGON I AND POLYGON J FACE EACH OTHER */
      /*==================================================*/
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
        printf("%spoly j is not facing poly i\n",indent);
      }
      ff = 0.0;
      enclosure->adaptive.num_zero++;
      enclosure->adaptive.num_zero_row++;
    } else if (enclosure->nonblocking) {
      /*======================================*/
      /* IF NONBLOCKING MESH (NO OCCLUSIONS), */
      /* DON'T CHECK FOR BLOCKING, JUST       */
      /* CALCULATE THE VIEWFACTOR             */
      /*======================================*/
      ff = VF_CalcVF_Unoccluded(&poly_ii, &poly_jj);
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
        printf("%ssimple mesh, ",indent);
        printf("calling VF_Unoccluded() = %g\n",ff);
      }
      enclosure->adaptive.num_analytic++;
      enclosure->adaptive.num_analytic_row++;
    } else {
      /*=================================================*/
      /* IF COMPLEX MESH (POSSIBLE OCCLUSIONS), EVALUATE */
      /* VIEW FACTOR EITHER ANALYTICALLY OR NUMERICALLY  */
      /* DEPENDING UPON PATCH-TO-PATCH VISIBLITY         */
      /*=================================================*/

      /*=============================================*/
      /* FIND THE LIST OF POSSIBLE OCCLUDING PATCHES */
      /*=============================================*/
      ncandidates = VF_FindCandidates(&poly_ii, &poly_jj,&candidates);
      if (ncandidates==0) {
        /*===================================*/
        /* NO POSSIBLE OCCLUDING PATCHES,    */
        /* EVALUATE VIEW FACTOR ANALYTICALLY */
        /*===================================*/
        ff = VF_CalcVF_Unoccluded(&poly_ii, &poly_jj);
        if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
          printf("%sno candidates, ",indent);
          printf("calling VF_Unoccluded() = %g\n",ff);
        }
        enclosure->adaptive.num_analytic++;
        enclosure->adaptive.num_analytic_row++;
      } else {
        /*==========================================*/
        /* POSSIBLE OCCLUDING PATCHES, EVALUATE THE */
        /* VIEW FACTORS ANALYTICALLY OR NUMERICALLY */
        /* DEPENDING UPON THEIR VISIBILITY          */
        /*==========================================*/
        if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
          printf("%snum candidates = %d\n",indent,ncandidates);
          printf("%s  ->  %d",indent,candidates[0].facet_num);
          for (n=1; n<ncandidates; n++) {
            if (n%10 == 0) 
              printf("\n%s      %d",indent,candidates[n].facet_num);
            else 
              printf(", %d",candidates[n].facet_num);
          }
          printf("\n");
        }
        if (poly_ii.np<=4) {
          npoly_i = 1;
          VF_PushPolyStack(&PolyStack_I, poly_ii);
        } else if (poly_ii.np==5) {
          if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
            printf("       setting npoly_i = 2\n");
          }
          npoly_i = 2;
          VF_PushSplitPolyStack(&PolyStack_I, &poly_ii);
        } else {
          printf("Whoa!!!  poly_ii->np = %d\n",poly_ii.np);
        }
        for (area=0.0, ff=0.0, i=0; i<npoly_i; i++) {
          VF_PopPolyStack(&PolyStack_I, poly_I);
          if (poly_jj.np<=4) {
            npoly_j = 1;
            VF_PushPolyStack(&PolyStack_J, poly_jj);
          } else if (poly_jj.np==5) {
            if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
              printf("%ssetting npoly_j = 2\n",indent);
            }
            npoly_j = 2;
            VF_PushSplitPolyStack(&PolyStack_J, &poly_jj);
          } else {
            printf("Whoa!!!  poly_jj->np = %d\n",poly_jj.np);
          }
          area += VF_PolyArea(&poly_I);
          for (j=0; j<npoly_j; j++) {
            VF_PopPolyStack(&PolyStack_J, poly_J);
            if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {

            }
#define NOT_NEW_WAY
#ifdef NEW_WAY
            ff += VF_PolyArea(&poly_I)*VF_CalcVF_Occluded(&poly_I, 
                                                          &poly_J, 
                                                          candidates, 
                                                          ncandidates,
                                                          visibility);
            enclosure->adaptive.num_numerical++;
            enclosure->adaptive.num_numerical_row++;
#else
            visibility = VF_Visibility(&poly_I, &poly_J,
                                       candidates, ncandidates);
            if (enclosure->debug_level>=VF_OUTPUT_DEBUG_3) {
              printf("           visibility = %g\n",visibility);
            }
            /*==========================================================*/
            /* CHECK IF POLYGON J IS COMPLETELY OCCLUDED FROM POLYGON I */
            /*==========================================================*/
            if (visibility==0.0) {
              ff += 0.0;
              enclosure->adaptive.num_zero++;
              enclosure->adaptive.num_zero_row++;
            /*===============================================*/
            /* CHECK IF POLYGON J IS COMPLETELY VISIBLE FROM */
            /* POLYGON I EVALUATE VIEW FACTOR ANALYTICALLY   */
            /*===============================================*/
            } else if (visibility==1.0) {
              if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
                printf("%s  no occlusions, calling VF_Unoccluded(), ",indent);
              }
              ff += VF_PolyArea(&poly_I)*VF_CalcVF_Unoccluded(&poly_I,&poly_J);
              enclosure->adaptive.num_analytic++;
              enclosure->adaptive.num_analytic_row++;
            /*==============================================*/
            /* CHECK IF POLYGON j IS PARTIALLY VISIBLE FROM */
            /* POLYGON I EVALUATE VIEW FACTOR ADAPTIVELY    */
            /*==============================================*/
            } else if (visibility>0.0) {
              if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {

              }
              ff += VF_PolyArea(&poly_I)*VF_CalcVF_Occluded(&poly_I, 
                                                            &poly_J, 
                                                            candidates, 
                                                            ncandidates,
                                                            visibility);
              enclosure->adaptive.num_numerical++;
              enclosure->adaptive.num_numerical_row++;
            }
#endif
          }
        }
        ff /= area;
      }
    }
  }
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
    printf("%sff = %g\n",indent,ff);
  }
  return ff;
}
