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
@(#)    $RCSfile: VF_SetViewPort.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_SetViewPort.c,v $
@(#)
@(#)    DESCRIPTION:  Set the ViewPort xform matrix for the hemicube method. 
@(#)
@(#)    Reference: "Computer Graphics, Systems and Concepts",
@(#)               Slater, R. and Salmon, M., Addison-Wesley, 1987.
@(#)               Chapter 13: Basics of 3D Graphics 
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <math.h>

#include "vf.h"

void VF_SetViewPort(ViewPort *view, int face)
{
  int    i, n=0;
  double tmp;
  Vector rx, ry, rz;
  Vector *normal0, *normal1;
  Matrix *m;

  view->up_vector   = view->u[face];
  view->view_normal = view->n[face];
  view->d           = -V3_Dot(&(view->view_normal),
                              &(view->view_point));
  m  = &(view->xform);
  rz = view->view_normal;
  ry = view->up_vector;
  V3_Negate(&rz);
  V3_Cross(&ry,&rz,&rx);
  V3_Normalize(&rx,tmp);
  V3_Normalize(&ry,tmp);
  V3_Normalize(&rz,tmp);

  /*========================================*/
  /* SETUP COORDINATE TRANSFORMATION MATRIX */
  /*========================================*/
  m->element[0][0] = rx.x;
  m->element[0][1] = rx.y;
  m->element[0][2] = rx.z;
  m->element[0][3] = -V3_Dot(&(view->view_point),&rx);

  m->element[1][0] = ry.x;
  m->element[1][1] = ry.y;
  m->element[1][2] = ry.z;
  m->element[1][3] = -V3_Dot(&(view->view_point),&ry);

  m->element[2][0] = rz.x;
  m->element[2][1] = rz.y;
  m->element[2][2] = rz.z;
  m->element[2][3] = -V3_Dot(&(view->view_point),&rz);
    
  m->element[3][0] = 0.0;
  m->element[3][1] = 0.0;
  m->element[3][2] = 0.0;
  m->element[3][3] = 1.0;
    
  if (face>0) n = 1;
  for (i=0; i<5; i++) {
    normal0    = &(view->frustum_planes0[n][i].normal);
    normal1    = &(view->frustum_planes1[i].normal);
    normal1->x = normal0->x * m->element[0][0] +
                 normal0->y * m->element[1][0] +
                 normal0->z * m->element[2][0];
    normal1->y = normal0->x * m->element[0][1] +
                 normal0->y * m->element[1][1] +
                 normal0->z * m->element[2][1];
    normal1->z = normal0->x * m->element[0][2] +
                 normal0->y * m->element[1][2] +
                 normal0->z * m->element[2][2];
    view->frustum_planes1[i].d = -V3_Dot(normal1, &(view->view_point));
  }
}
