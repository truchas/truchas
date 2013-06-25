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
@(#)    $RCSfile: VF_ProjectOntoHemicube.c,v $
@(#)    $Revision: 1.5 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_ProjectOntoHemicube.c,v $
@(#)
@(#)    DESCRIPTION:  Project a surface onto a hemicube.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void VF_ProjectOntoHemicube(Poly *poly_i, Poly *poly_j, int patch_j,
                            ViewPort *view, int inview, int face, 
                            Hemicube *hc)
{  
  int	 i, clip;
  double scale, scale_x, scale_y, tmp;
  Poly   pp, pp0;
  VFtopology  *topology=VF_CurrentTopology();

  /*=====================================================*/
  /* CHECK IF POLYGON J IS BEHIND OR NOT FACING POYGON I */
  /*=====================================================*/
  if (VF_BehindAndBackFaceCullViewPoint(view, poly_j)) return;
  /*================================================*/
  /* PERFORM COORDINATE TRANSFORMATION OF POLYGON J */
  /*================================================*/
  VF_TransformPoly(poly_j, &pp0,  &(view->xform));
  if (topology->geom==VF_2Dplanar) {
    VF_PolyNormal_Aux(&pp0, &(pp0.normal), face);
    pp0.d = -V3_Dot(&(pp0.normal), &(pp0.p[0]));
  }
  /*========================================================*/
  /* DO POLYGON CLIPPING DEPENDING UPON WHETHER OR NOT THE  */
  /* CURRENT BSP NODE IS EITHER PARTIALLY OR COMPLETELY     */
  /* INSIDE THE VIEW FRUSTUM AND WHETHER OR NOT THE POLYGON */
  /* IS EITHER PARTIALLY OR COMPLETELY INSIDE THE BSP NODE  */
  /*========================================================*/
  switch (inview) {
  case VF_FRUSTUM_CLIP_PARTIAL:
    clip = VF_ClipToFrustum(&pp0, face);
    break;
  case VF_FRUSTUM_CLIP_IN:
    if (poly_j->mask&VF_FACET_MASK_SPLIT) {
      clip = VF_ClipToFrustum(&pp0, face);
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
  /* COMPLETELY OUTSIDE THE VIEW FRUSTUM    */
  /*========================================*/
  if (clip!=VF_POLY_CLIP_OUT) {
    /*===============================================*/
    /* PROJECT POLYGON (X,Y) COORDINATES TO THE VIEW */
    /* PLANE AND RESCALE TO THE HEMICUBE RESOLUTION  */
    /*===============================================*/
    if (face==0) {
      /*==========================*/
      /* x = [-1,+1] -> [0,win_x] */
      /* y = [-1,+1] -> [0,win_y] */
      /* with win_x = win_y	  */
      /*==========================*/
      scale = 0.5*(view->window.x1+1);
      for (pp.np=pp0.np, i=0; i<pp.np; i++) {
        tmp	  = -scale/pp0.p[i].z;
        pp.p[i].x = pp0.p[i].x*tmp+scale;
        pp.p[i].y = pp0.p[i].y*tmp+scale;
        pp.p[i].z = pp0.p[i].z;
      }
    } else {
      /*==========================*/
      /* x = [-1,+1] -> [0,win_x] */
      /* y = [ 0,+1] -> [0,win_y] */
      /* with win_x = 2*win_y	  */
      /*==========================*/
      scale_x = 0.5*(view->window.x1+1);
      scale_y = view->window.y1+1;
      for (pp.np=pp0.np, i=0; i<pp.np; i++) {
        tmp	  = -1.0/pp0.p[i].z;
        pp.p[i].x = scale_x*(pp0.p[i].x*tmp+1.0);
        pp.p[i].y = scale_y*(pp0.p[i].y*tmp    );
        pp.p[i].z = pp0.p[i].z;
      }
    }
    /*==========================*/
    /* SCAN CONVERT THE POLYGON */
    /*==========================*/
    VF_PolyScan(&pp0,&pp,&(view->window),patch_j,face,hc);
  }
}

