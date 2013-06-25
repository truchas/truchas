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
@(#)    $RCSfile: VF_BackFaceCull.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_BackFaceCull.c,v $
@(#)
@(#)    DESCRIPTION:  Several culling utilities to check if surfaces 
@(#)    are behind or not facing each other.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vf.h"

/*===========================================*/
/* CHECK IF TWO POLYGONS CAN SEE EACH OTHER  */
/* RETURN: 0 IF VISIBLE AND 1 IF NOT VISIBLE */
/*===========================================*/
int VF_BackFaceCullPolys(Poly *poly_i, Poly *poly_j)
{
  int    n, m;
  double tmp;
  Vector *normal_i, *normal_j, r;

  normal_i = &(poly_i->normal);
  normal_j = &(poly_j->normal);
  for (m=0; m<poly_i->np; m++) {
    for (n=0; n<poly_j->np; n++) {
      V3_Sub(&(poly_j->p[n]), &(poly_i->p[m]), &r);
      V3_Normalize(&r,tmp);
      if ((V3_Dot(&r, normal_i)>0.0) && (-V3_Dot(&r, normal_j)>0.0)) return 0;
    }
  }
  return 1;
}

/*===========================================*/
/* CHECK IF POLYGON CAN SEE A VIEWPOINT      */
/* RETURN: 0 IF VISIBLE AND 1 IF NOT VISIBLE */
/*===========================================*/
int VF_BackFaceCullViewPoint(ViewPort *view, Poly *poly_j)
{
  int    n;
  Vector *normal_i, *normal_j, r;

  normal_i = &(view->view_normal);
  normal_j = &(poly_j->normal);
  for (n=0; n<poly_j->np; n++) {
    V3_Sub(&(poly_j->p[n]), &(view->view_point), &r);
    if ((V3_Dot(&r, normal_i)>0.0) && (-V3_Dot(&r, normal_j)>0.0)) return 0;
  }
  return 1;
}

/*=========================================*/
/* CHECK IF A POLY J IS BEHIND POLY I      */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_BehindPoly(Poly *poly_i, Poly *poly_j)
{
  int    n;
  double value, spatial_tol;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  for (n=0; n<poly_j->np; n++) {
    value = V3_Dot(&(poly_j->p[n]), &(poly_i->normal)) + poly_i->d;
    if (value > spatial_tol) return 0;
  }
  return 1;
}

/*=========================================*/
/* CHECK IF A POLYGON IS BEHIND A PLANE    */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_BehindPlane(Plane *plane, Poly *poly)
{
  int    n;
  double value, spatial_tol;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  for (n=0; n<poly->np; n++) {
    value = V3_Dot(&(poly->p[n]), &(plane->normal)) + plane->d;
    if (value > spatial_tol) return 0;
  }
  return 1;
}

/*=========================================*/
/* CHECK IF A PATCH IS BEHIND A PLANE      */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_BehindViewPoint(ViewPort *view, Poly *poly)
{
  int    n;
  double value, spatial_tol;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  for (n=0; n<poly->np; n++) {
    value = V3_Dot(&(poly->p[n]), &(view->view_normal)) + view->d;
    if (value > spatial_tol) return 0;
  }
  return 1;
}

/*=========================================*/
/* CHECK IF A POLY IS BEHIND A PLANE       */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_BehindAndBackFaceCullViewPoint(ViewPort *view, Poly *poly)
{
  int    n;
  double value, spatial_tol;
  Vector *normal_i, *normal_j, r;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  for (n=0; n<poly->np; n++) {
    value = V3_Dot(&(poly->p[n]), &(view->view_normal)) + view->d;
    if (value > spatial_tol) {
      /*===========================================*/
      /* THIS VERTEX IS IN FRONT OF THE VIEWPOINT, */
      /* NOW CHECK IF IT IS ALSO FACING IT         */
      /*===========================================*/
      normal_i = &(view->view_normal);
      normal_j = &(poly->normal);
      V3_Sub(&(poly->p[n]), &(view->view_point), &r);
      if ((V3_Dot(&r, normal_i)>0.0) && (-V3_Dot(&r, normal_j)>0.0)) return 0;
      break;
    }
  }
  return 1;
}

/*=========================================*/
/* CHECK IF A VOXEL IS BEHIND A PLANE      */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_VoxelBehindPlane(Box *bounds, Plane *plane)
{
  double spatial_tol, d;
  double xMax, yMax, zMax;
  Vector *normal;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  normal      = &(plane->normal);
  d           = plane->d;
    
  xMax = MAX( bounds->xmin * normal->x, bounds->xmax * normal->x );
  yMax = MAX( bounds->ymin * normal->y, bounds->ymax * normal->y );
  zMax = MAX( bounds->zmin * normal->z, bounds->zmax * normal->z );
  if (xMax + yMax + zMax + d > spatial_tol) return FALSE;
  return TRUE;
}

/*=========================================*/
/* CHECK IF A VOXEL IS BEHIND A PATCH      */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_VoxelBehindFacet(Box *bounds, Facet *facet)
{
  double spatial_tol;
  double xMax, yMax, zMax;
  Vector *normal;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  normal      = &(facet->normal);
    
  xMax = MAX( bounds->xmin * normal->x, bounds->xmax * normal->x );
  yMax = MAX( bounds->ymin * normal->y, bounds->ymax * normal->y );
  zMax = MAX( bounds->zmin * normal->z, bounds->zmax * normal->z );
  if (xMax + yMax + zMax + facet->d > spatial_tol) return FALSE;
  return TRUE;
}


/*=========================================*/
/* CHECK IF A VOXEL IS BEHIND A POLYGON    */
/* RETURN: 0 IF NOT BEHIND AND 1 IF BEHIND */
/*=========================================*/
int VF_VoxelBehindPoly(Box *bounds, Poly *poly)
{
  double spatial_tol;
  double xMax, yMax, zMax;
  Vector *normal;
  VFtopology *topology=VF_CurrentTopology();

  spatial_tol = topology->spatial_tol;
  normal      = &(poly->normal);
    
  xMax = MAX( bounds->xmin * normal->x, bounds->xmax * normal->x );
  yMax = MAX( bounds->ymin * normal->y, bounds->ymax * normal->y );
  zMax = MAX( bounds->zmin * normal->z, bounds->zmax * normal->z );
  if (xMax + yMax + zMax + poly->d > spatial_tol) return FALSE;
  return TRUE;
}
