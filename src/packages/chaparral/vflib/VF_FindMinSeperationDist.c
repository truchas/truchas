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
@(#)    $RCSfile: VF_FindMinSeperationDist.c,v $
@(#)    $Revision: 1.6 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_FindMinSeperationDist.c,v $
@(#)
@(#)    DESCRIPTION:  Use hemicube method to calculate one row of
@(#)    the viewfactor matrix.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

#define DO_CLIP_TO_FRUSTUM  1
#define DO_FRONT_TO_BACK    0

double VF_FindMinSeperationDist(Poly *poly_i, int facet_i, int subfacet_i)
{
  int         clip, facet_j, subfacet_j;
  int         inview, nsub_polys_j;
  float       dist;
  double      min_dist=MAXDOUBLE;
  Plane       view_plane;
  Point       center;
  Poly        *poly_j, pp0;
  POLYstack   PolyStack;
  BinNodePtr  currentNode;
  BSPstackPtr stack;
  Facet       *facet;
  ViewPort    view;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  stack = topology->BSP_stack;
  VF_BSP_InitStack(stack);
  VF_InitPolyStack(&PolyStack);
  /*==============================================*/
  /* SET THE VIEWPORT AND DEFINE THE VIEW PLANE   */
  /*==============================================*/
  VF_SetView(&view, poly_i);
  VF_SetViewPort(&view, 0);
  view_plane.normal = view.view_normal;
  view_plane.d      = -V3_Dot(&(view_plane.normal), &view.view_point);
  /*===================================================*/
  /* TAG ALL THE BSP NODES THAT ARE INFRONT/BEHIND THE */
  /* PATCH AND COMPUTE THE TRANSFORMED BOUNDING BOX OF */
  /* THE BSP NODES THAT ARE IN FRONT OF THE PATCH      */
  /*===================================================*/
  currentNode = topology->BSP_tree.root;
  while (currentNode != NULL) {
    currentNode->behind = VF_VoxelBehindPlane(&(currentNode->bounds), 
                                              &view_plane);
    if (!currentNode->behind) {
      if (currentNode->child[0]!=NULL) {
        VF_BSP_PushStack(stack, currentNode->child[0]);
      }
      if (currentNode->child[1]!=NULL) {
        VF_BSP_PushStack(stack, currentNode->child[1]);
      }
    }
    VF_BSP_PopStack(stack, currentNode);
  }
  /*=======================*/
  /* TRAVERSE THE BSP TREE */
  /*=======================*/
  VF_BSP_InitStack(stack);
  currentNode = topology->BSP_tree.root;
  currentNode->clipped = VF_FRUSTUM_CLIP_UNKNOWN;
  while (currentNode != NULL) {
    inview = VF_FRUSTUM_CLIP_OUT;
    if (!currentNode->behind) {
      if (DO_CLIP_TO_FRUSTUM) {
  	if (currentNode->clipped==VF_FRUSTUM_CLIP_IN) {
  	  inview = VF_FRUSTUM_CLIP_IN;
  	} else {
  	  inview = VF_ClipBoxToFrustum(&view,&(currentNode->bounds),0);
  	}
      } else {
  	inview = VF_FRUSTUM_CLIP_PARTIAL;
      }
    }
    if (inview!=VF_FRUSTUM_CLIP_OUT) {
      if (currentNode->child[0]==NULL && currentNode->child[1]==NULL) {
  	/*======================================*/
  	/* NO MORE BSP NODES BELOW THIS ONE SO  */
  	/* PROJECT ALL THE PATCHES IN THIS NODE */
  	/*======================================*/
  	facet = VF_FirstFacetOfLinkList(&currentNode->members);
  	while (facet) {
  	  facet_j = facet->index;
  	  VF_FacetToPolyStack(facet, &PolyStack, &nsub_polys_j);
  	  for (subfacet_j=0; subfacet_j<nsub_polys_j; subfacet_j++) {
  	    VF_PopPolyStackPtr(&PolyStack, poly_j);
  	    if (facet_i==facet_j && subfacet_i==subfacet_j) continue;
            /*=====================================================*/
            /* CHECK IF POLYGON J IS BEHIND OR NOT FACING POYGON I */
            /*=====================================================*/
            if (!VF_BehindAndBackFaceCullViewPoint(&view, poly_j)) {
              /*================================================*/
              /* PERFORM COORDINATE TRANSFORMATION OF POLYGON J */
              /*================================================*/
              VF_TransformPoly(poly_j, &pp0,  &(view.xform));
              if (topology->geom==VF_2Dplanar) {
                VF_PolyNormal_Aux(&pp0, &(pp0.normal), 0);
                pp0.d = -V3_Dot(&(pp0.normal), &(pp0.p[0]));
              }
              /*========================================================*/
              /* DO POLYGON CLIPPING DEPENDING UPON WHETHER OR NOT THE  */
              /* CURRENT BSP NODE IS EITHER PARTIALLY OR COMPLETELY	*/
              /* INSIDE THE VIEW FRUSTUM AND WHETHER OR NOT THE POLYGON */
              /* IS EITHER PARTIALLY OR COMPLETELY INSIDE THE BSP NODE  */
              /*========================================================*/
              switch (inview) {
              case VF_FRUSTUM_CLIP_PARTIAL:
            	clip = VF_ClipToFrustum(&pp0, 0);
            	break;
              case VF_FRUSTUM_CLIP_IN:
            	if (poly_j->mask&VF_FACET_MASK_SPLIT) {
          	  clip = VF_ClipToFrustum(&pp0, 0);
            	} else {
          	  clip = VF_POLY_CLIP_IN;
            	}
            	break;
              default:
            	clip = VF_POLY_CLIP_OUT;
            	break;
              }
              /*========================================*/
              /* SCAN CONVERT THE POLYGON IF IT WAS NOT */
              /* COMPLETELY OUTSIDE THE VIEW FRUSTUM	*/
              /*========================================*/
              if (clip!=VF_POLY_CLIP_OUT) {
            	/*====================================*/
            	/* COMPUTE DISTANCE BETWEEN CENTROIDS */
            	/*====================================*/
            	VF_PolyCenter(&pp0,&center);
            	min_dist = MIN(min_dist, V3_Length(&center));
              }
	    }
	  }
  	  facet = VF_NextFacetOfLinkList(&currentNode->members);
  	}
      } else {
  	if (currentNode->clipped==VF_FRUSTUM_CLIP_IN) {
  	  if (currentNode->child[0]!=NULL) {
  	    currentNode->child[0]->clipped = VF_FRUSTUM_CLIP_IN;
  	  }
  	  if (currentNode->child[1]!=NULL) {
  	    currentNode->child[1]->clipped = VF_FRUSTUM_CLIP_IN;
  	  }
  	} else {
  	  if (currentNode->child[0]!=NULL) {
  	    currentNode->child[0]->clipped = VF_FRUSTUM_CLIP_UNKNOWN;
  	  }
  	  if (currentNode->child[1]!=NULL) {
  	    currentNode->child[1]->clipped = VF_FRUSTUM_CLIP_UNKNOWN;
  	  }
  	}
  	switch (currentNode->axis) {
  	case VF_BSP_XAXIS:
  	  dist = SGN(view.view_normal.x);
  	  break;
  	case VF_BSP_YAXIS:
  	  dist = SGN(view.view_normal.y);
  	  break;
  	case VF_BSP_ZAXIS:
  	  dist = SGN(view.view_normal.z);
  	  break;
  	}
  	if (dist>0) {
  	  if (currentNode->child[1]!=NULL) {
  	    VF_BSP_PushStack(stack, currentNode->child[1]); /* far  */
  	  }
  	  if (currentNode->child[0]!=NULL) {
  	    VF_BSP_PushStack(stack, currentNode->child[0]); /* near */
  	  }
  	} else {
  	  if (currentNode->child[0]!=NULL) {
  	    VF_BSP_PushStack(stack, currentNode->child[0]); /* far  */
  	  }
  	  if (currentNode->child[1]!=NULL) {
  	    VF_BSP_PushStack(stack, currentNode->child[1]); /* near */
  	  }
  	}
      }
    }
    VF_BSP_PopStack(stack, currentNode);
  }
  return(min_dist);
}
