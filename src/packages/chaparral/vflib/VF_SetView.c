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
@(#)    $RCSfile: VF_SetView.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_SetView.c,v $
@(#)
@(#)    DESCRIPTION:  Set the view for the hemicube faces.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void VF_SetView(ViewPort *view, Poly *poly)
{
  Vector       offset, u, v, c;
  VFtopology*  topology=VF_CurrentTopology();
  VFenclosure* enclosure=VF_CurrentEnclosure();

  offset          = poly->normal;
  view->window.x0 = 0;
  view->window.y0 = 0;
  view->window.x1 = enclosure->hemicube.resolution-1;
  view->window.y1 = enclosure->hemicube.resolution-1;
    
  VF_PolyCenter(poly,&c);
  V3_Scale(&offset, topology->spatial_tol);
  V3_Add(&offset, &c, &(view->view_point));
    
  VF_JitterHemicube(&(poly->normal), &u, &v, poly);
  view->u[0] = v;
  view->u[1] = poly->normal;
  view->u[2] = poly->normal;
  view->u[3] = poly->normal;
  view->u[4] = poly->normal;
  view->n[0] = poly->normal;
  view->n[1] = u;
  view->n[2] = u;
  view->n[3] = v;
  view->n[4] = v;
  V3_Negate(&(view->n[2]));
  V3_Negate(&(view->n[4]));
  view->frustum_planes0[0][0].normal.x =  0.0;
  view->frustum_planes0[0][0].normal.y =  0.0;
  view->frustum_planes0[0][0].normal.z = -1.0;
  view->frustum_planes0[0][0].d        =  0.0;
  view->frustum_planes0[0][1].normal.x = -SQ2;
  view->frustum_planes0[0][1].normal.y =  0.0;
  view->frustum_planes0[0][1].normal.z = -SQ2;
  view->frustum_planes0[0][1].d        =  0.0;
  view->frustum_planes0[0][2].normal.x =  SQ2;
  view->frustum_planes0[0][2].normal.y =  0.0;
  view->frustum_planes0[0][2].normal.z = -SQ2;
  view->frustum_planes0[0][2].d        =  0.0;
  view->frustum_planes0[0][3].normal.x =  0.0;
  view->frustum_planes0[0][3].normal.y = -SQ2;
  view->frustum_planes0[0][3].normal.z = -SQ2;
  view->frustum_planes0[0][3].d        =  0.0;
  view->frustum_planes0[0][4].normal.x =  0.0;
  view->frustum_planes0[0][4].normal.y =  SQ2;
  view->frustum_planes0[0][4].normal.z = -SQ2;
  view->frustum_planes0[0][4].d        =  0.0;
    
  view->frustum_planes0[1][0].normal.x =  0.0;
  view->frustum_planes0[1][0].normal.y =  0.0;
  view->frustum_planes0[1][0].normal.z = -1.0;
  view->frustum_planes0[1][0].d        =  0.0;
  view->frustum_planes0[1][1].normal.x = -SQ2;
  view->frustum_planes0[1][1].normal.y =  0.0;
  view->frustum_planes0[1][1].normal.z = -SQ2;
  view->frustum_planes0[1][1].d        =  0.0;
  view->frustum_planes0[1][2].normal.x =  SQ2;
  view->frustum_planes0[1][2].normal.y =  0.0;
  view->frustum_planes0[1][2].normal.z = -SQ2;
  view->frustum_planes0[1][2].d        =  0.0;
  view->frustum_planes0[1][3].normal.x =  0.0;
  view->frustum_planes0[1][3].normal.y = -SQ2;
  view->frustum_planes0[1][3].normal.z = -SQ2;
  view->frustum_planes0[1][3].d        =  0.0;
  view->frustum_planes0[1][4].normal.x =  0.0;
  view->frustum_planes0[1][4].normal.y =  1.0;
  view->frustum_planes0[1][4].normal.z =  0.0;
  view->frustum_planes0[1][4].d        =  0.0;
}
