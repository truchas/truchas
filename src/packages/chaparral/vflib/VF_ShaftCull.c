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
@(#)    $RCSfile: VF_ShaftCull.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_ShaftCull.c,v $
@(#)
@(#)    DESCRIPTION:  Shaft cull test.
@(#)
@(#)    REFERENCE:
@(#)    "Shaft Culling for Efficient Ray-Traced Radiosity," 
@(#)    Eric A. Haines and John R. Wallace, Photorealistic Rendering in 
@(#)    Computer Graphics (Proceedings of the Second Eurographics Workshop
@(#)    on Rendering), Springer-Verlag, New York, 1994, p.122-138. Also in 
@(#)    SIGGRAPH '91 Frontiers in Rendering course notes.
@(#)
@(#)    paper available online at: 
@(#)    http://wuarchive.wustl.edu/graphics/graphics/mirrors/ftp.funet.fi/pub/sci/papers/graphics/hain91.ps.Z
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vf.h"

void VF_CreateShaft(Shaft *shaft, Poly *poly_i, Poly *poly_j)
{
  int    i, j, dir1, dir2, dir3;
  double u1, u2, v1, v2, du, dv;
  Box    minmax;
  Plane  *plane;
    
  VF_ComputePolyBoundingBox(&shaft->box_i, poly_i);
  VF_ComputePolyBoundingBox(&shaft->box_j, poly_j);
  VF_ExtentsBoundingBox(&shaft->box_i, &shaft->box_j, &shaft->extent, &minmax);
  shaft->nplanes = 0;
  for (i=0; i<6; i++) {
    dir1 = (i%3);
    switch (i) {
    case 0:
      if (minmax.xmin>0.0) continue;
      u1 = shaft->box_i.xmin;
      u2 = shaft->box_j.xmin;
      break;
    case 1:
      if (minmax.ymin>0.0) continue;
      u1 = shaft->box_i.ymin;
      u2 = shaft->box_j.ymin;
      break;
    case 2:
      if (minmax.zmin>0.0) continue;
      u1 = shaft->box_i.zmin;
      u2 = shaft->box_j.zmin;
      break;
    case 3:
      if (minmax.xmax>0.0) continue;
      u1 = shaft->box_i.xmax;
      u2 = shaft->box_j.xmax;
      break;
    case 4:
      if (minmax.ymax>0.0) continue;
      u1 = shaft->box_i.ymax;
      u2 = shaft->box_j.ymax;
      break;
    case 5:
      if (minmax.zmax>0.0) continue;
      u1 = shaft->box_i.zmax;
      u2 = shaft->box_j.zmax;
      break;
    }
    for (j=0; j<6; j++) {
      dir2 = (j%3);
      dir3 = 3-dir1-dir2;
      if (dir1==dir2) continue;
      switch (j) {
      case 0:
        if (minmax.xmin<0.0) continue;
        v1 = shaft->box_i.xmin;
        v2 = shaft->box_j.xmin;
        break;
      case 1:
        if (minmax.ymin<0.0) continue;
        v1 = shaft->box_i.ymin;
        v2 = shaft->box_j.ymin;
        break;
      case 2:
        if (minmax.zmin<0.0) continue;
        v1 = shaft->box_i.zmin;
        v2 = shaft->box_j.zmin;
        break;
      case 3:
        if (minmax.xmax<0.0) continue;
        v1 = shaft->box_i.xmax;
        v2 = shaft->box_j.xmax;
        break;
      case 4:
        if (minmax.ymax<0.0) continue;
        v1 = shaft->box_i.ymax;
        v2 = shaft->box_j.ymax;
        break;
      case 5:
        if (minmax.zmax<0.0) continue;
        v1 = shaft->box_i.zmax;
        v2 = shaft->box_j.zmax;
        break;
      }
      if ((i<=2 && j<=2) || (i>=3 && j>=3)) {
        du = v2-v1;
        dv = u1-u2;
      } else { /* normal must point outwards shaft */
        du = v1-v2;
        dv = u2-u1;
      }
      plane = &shaft->plane[shaft->nplanes];
      switch (dir1) {
      case 0:
        plane->normal.x = du;
        break;
      case 1:
        plane->normal.y = du;
        break;
      case 2:
        plane->normal.z = du;
        break;
      default:
        printf("VF_CreateShaft(): bad value for dir1 (%d)\n", dir1);
        break;
      }
      switch (dir2) {
      case 0:
        plane->normal.x = dv;
        break;
      case 1:
        plane->normal.y = dv;
      break;
      case 2:
        plane->normal.z = dv;
        break;
      default:
        printf("VF_CreateShaft(): bad value for dir2 (%d)\n", dir2);
        break;
      }
      switch (dir3) {
      case 0:
        plane->normal.x = 0.0;
        break;
      case 1:
        plane->normal.y = 0.0;
        break;
      case 2:
        plane->normal.z = 0.0;
        break;
      default:
        printf("VF_CreateShaft(): bad value for dir3 (%d)\n", dir3);
        break;
      }
      plane->d = -(du*u1 + dv*v1);
      shaft->nplanes++;
    }
  }
}

/*================================================================*/
/* "returns" nonzero if the two given bounding boxes are disjunct */
/*================================================================*/
#define DisjunctBoundsOld(b1, b2) ((b1->xmin >= b2->xmax) || \
                                   (b2->xmin >= b1->xmax) || \
                                   (b1->ymin >= b2->ymax) || \
                                   (b2->ymin >= b1->ymax) || \
                                   (b1->zmin >= b2->zmax) || \
                                   (b2->zmin >= b1->zmax))
                                   
#define DisjunctBounds3D(b1, b2, t) ((b1->xmin+t >= b2->xmax) || \
                                     (b2->xmin+t >= b1->xmax) || \
                                     (b1->ymin+t >= b2->ymax) || \
                                     (b2->ymin+t >= b1->ymax) || \
                                     (b1->zmin+t >= b2->zmax) || \
                                     (b2->zmin+t >= b1->zmax))
                                   
#define DisjunctBounds2D(b1, b2, t) ((b1->xmin+t >= b2->xmax) || \
                                     (b2->xmin+t >= b1->xmax) || \
                                     (b1->ymin+t >= b2->ymax) || \
                                     (b2->ymin+t >= b1->ymax))

/*===========================================*/
/* RETURNS XXX                               */
/*         |||                               */
/*         ||+--- IF INTERSECTING WITH BOX_I */
/*         ||                                */
/*         |+---- IF INTERSECTING WITH BOX_J */
/*         |                                 */
/*         +----- IF INTERSECTING WITH SHAFT */
/*===========================================*/
int VF_ShaftCull(Box *bounds, Shaft *shaft, double tol)
{
  int   i=0, outside=0, status=0;
  Plane *plane;
  Point close;
  VFtopology *topology=VF_CurrentTopology();
    
  /*======================================*/
  /* TEST AGAINST THE EXTENT BOUNDING BOX */
  /*======================================*/
  if (topology->geom==VF_2Dplanar) {
    if (DisjunctBounds2D(bounds, (&shaft->extent), tol)) return 0;
  } else {
    if (DisjunctBounds3D(bounds, (&shaft->extent), tol)) return 0;
  }

  /*===============================================*/
  /* TEST AGAINST THE TWO REFERENCE BOUNDING BOXES */
  /*===============================================*/
  if (topology->geom==VF_2Dplanar) {
    if (!DisjunctBounds2D(bounds, (&shaft->box_i), tol)) status |= 0x1;
  } else {
    if (!DisjunctBounds3D(bounds, (&shaft->box_i), tol)) status |= 0x1;
  }
  if (topology->geom==VF_2Dplanar) {
    if (!DisjunctBounds2D(bounds, (&shaft->box_j), tol)) status |= 0x2;
  } else {
    if (!DisjunctBounds3D(bounds, (&shaft->box_j), tol)) status |= 0x2;
  }
  if (status>0) return status;

  /*==================================================*/
  /* TEST AGAINST THE SHAFT PLANE SET:                */
  /* IF CLOSEST CORNER OF THE BOUNDING BOX IS OUTSIDE */ 
  /* ANY SHAFT-PLANE, THE OBJECT IS OUTSIDE THE SHAFT */
  /*==================================================*/
  if (shaft->nplanes==0) return status;
  for (i=0; i<shaft->nplanes; i++) {
    plane = &shaft->plane[i];
    if (topology->geom==VF_2Dplanar && plane->normal.z != 0.0) continue;
    close.x = plane->normal.x < 0.0 ? bounds->xmax : bounds->xmin;
    close.y = plane->normal.y < 0.0 ? bounds->ymax : bounds->ymin;
    close.z = plane->normal.z < 0.0 ? bounds->zmax : bounds->zmin;
    if ((V3_Dot(&(plane->normal),&close)+plane->d)>=0.0) outside++;
  }
  if (outside==shaft->nplanes) return status;
  status |= 0x4;
  return status;
}
