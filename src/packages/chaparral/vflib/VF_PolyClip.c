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
@(#)    $RCSfile: VF_PolyClip.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_PolyClip.c,v $
@(#)
@(#)    DESCRIPTION:  Polygon clipping routines based on the routines from
@(#)    the section on "Generic Convex Polygon Scan Conversion and Clipping"
@(#)    by Paul Heckbert from "Graphics Gems", Academic Press, 1990.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "vf.h"

#define CLIP_AND_SWAP(p, q, r, d) { \
  VF_ClipToFrustumPlane(p, q, d); \
  if (q->np==0) {p1->np = 0; return VF_POLY_CLIP_OUT;} \
  r = p; p = q; q = r;}

int VF_ClipToFrustum(Poly *p1, int face)
{
  int     i, in=0;
  int     x0out=0, x1out=0;
  int     y0out=0, y1out=0;
  double  d[VF_POLY_NMAX], ztest, n_dot_d;
  Poly    p2, *p, *q, *r;
  
  p         = p1;
  q         = &p2;
  q->d      = p->d;
  q->normal = p->normal;
  /*====================================*/
  /* FIRST CHECK TO SEE IF THE PLANE OF */
  /* THE POLYGON CONTAINS THE VIEWPOINT */
  /*====================================*/
  n_dot_d = V3_Dot(&(p->normal), &(p->p[0]));
  if (fabs(n_dot_d)<VF_RAY_TEST) return VF_POLY_CLIP_OUT;
  /*========================================*/
  /* NEXT, CHECK TO SEE IF ALL THE VERTICES */
  /* ARE INSIDE OR OUTSIDE THE VIEW FRUSTUM */
  /*========================================*/
  if (face==0) {
    for (i=0; i<p->np; i++) {
      ztest = p->p[i].z;
      /* TEST AGAINST X = Z PLANE */
      if (p->p[i].x<ztest) {
        x0out++;
      } else {
        in++;
      }
      /* TEST AGAINST X = -Z PLANE */
      if (p->p[i].x>-ztest) {
        x1out++;
      } else {
        in++;
      }
      /* TEST AGAINST Y = Z PLANE */
      if (p->p[i].y<ztest) {
        y0out++;
      } else {
        in++;
      }
      /* TEST AGAINST Y = -Z PLANE */
      if (p->p[i].y>-ztest) {
        y1out++;
      } else {
        in++;
      }
    }
  } else {
    for (i=0; i<p->np; i++) {
      ztest = p->p[i].z;
      /* TEST AGAINST X = Z PLANE */
      if (p->p[i].x<ztest) {
        x0out++;
      } else {
        in++;
      }
      /* TEST AGAINST X = -Z PLANE */
      if (p->p[i].x>-ztest) {
        x1out++;
      } else {
        in++;
      }
      /* TEST AGAINST Y = 0 PLANE */
      if (p->p[i].y<0.0) {
        y0out++;
      } else {
        in++;
      }
      /* TEST AGAINST Y = -Z PLANE */
      if (p->p[i].y>-ztest) {
        y1out++;
      } else {
        in++;
      }
    }
  }
  if (x0out==p->np || x1out==p->np ||
      y0out==p->np || y1out==p->np) return VF_POLY_CLIP_OUT;
  if (in==4*p->np) return VF_POLY_CLIP_IN;
    
  /*===========================================*/
  /* CLIP THE POLYGONS THAT ARE NOT COMPLETELY */
  /* INSIDE OR OUTSIDE THE VIEW FRUSTUM        */
  /*===========================================*/
    
  /*=============*/
  /* X = Z PLANE */
  /*=============*/
  for (i=0; i<p->np; i++) {
    d[i] = SQ2*(-p->p[i].z - -p->p[i].x);
    /*if (fabs(d[i])<mesh->spatial_tol) d[i] = 0.0;*/
  }
  CLIP_AND_SWAP(p, q, r, d);
    
  /*==============*/
  /* X = -Z PLANE */
  /*==============*/
  for (i=0; i<p->np; i++) {
    d[i] = SQ2*(-p->p[i].z - p->p[i].x);
    /*if (fabs(d[i])<mesh->spatial_tol) d[i] = 0.0;*/
  }
  CLIP_AND_SWAP(p, q, r, d);
    
  if (face==0) {
    /*==============*/
    /* Y = Z PLANE */
    /*==============*/
    for (i=0; i<p->np; i++) {
      d[i] = SQ2*(-p->p[i].z - -p->p[i].y);
      /*if (fabs(d[i])<mesh->spatial_tol) d[i] = 0.0;*/
    }
  } else {
    /*=============*/
    /* Y = 0 PLANE */
    /*=============*/
    for (i=0; i<p->np; i++) {
      d[i] = p->p[i].y;
      /*if (fabs(d[i])<mesh->spatial_tol) d[i] = 0.0;*/
    }
  }
  CLIP_AND_SWAP(p, q, r, d);
    
  /*=============*/
  /* Y = -Z PLANE */
  /*=============*/
  for (i=0; i<p->np; i++) {
    d[i] = SQ2*(-p->p[i].z - p->p[i].y);
    /*if (fabs(d[i])<mesh->spatial_tol) d[i] = 0.0;*/
  }
  CLIP_AND_SWAP(p, q, r, d);
    
  if (p==&p2) {
    memcpy(p1, &p2, sizeof(Poly)-(VF_POLY_NMAX-p2.np)*sizeof(Point));
  }
  return VF_POLY_CLIP_PARTIAL;
}

void VF_ClipToFrustumPlane(Poly *p, Poly *q, double d[])
{
  int    i, n;
  double t, tu, tv;
  Vector *u, *v, *w;

  q->np     = 0;
  q->d      = p->d;
  q->normal = p->normal;
  if (p->np==2) {
    u  = &(p->p[1]);
    tu = d[1];
    v  = &(p->p[0]);
    tv = d[0];
    /*==============================*/
    /* VERTEX V IS IN, COPY IT TO Q */
    /*==============================*/
    if (tv>=0.0) {
      w    = &(q->p[q->np]);
      w->x = v->x;
      w->y = v->y;
      w->z = v->z;
      q->np++;
    }
    if ((tu>0.0 && tv<0.0) ||
        (tu<0.0 && tv>0.0)) {
      /*=================================================*/
      /* EDGE CROSSES PLANE; ADD INTERSECTION POINT TO Q */
      /*=================================================*/
      t      = tu/(tu-tv);
      w      = &(q->p[q->np]);
      w->x   = u->x+t*(v->x-u->x);
      w->y   = u->y+t*(v->y-u->y);
      w->z   = u->z+t*(v->z-u->z);
      q->np++;
    }
    /*==============================*/
    /* VERTEX V IS IN, COPY IT TO Q */
    /*==============================*/
    if (tu>=0.0) {
      w    = &(q->p[q->np]);
      w->x = u->x;
      w->y = u->y;
      w->z = u->z;
      q->np++;
    }
    if (q->np<2) q->np=0;
  } else {
    /*========================*/
    /* START WITH U=VERT[N-1] */
    /*========================*/
    u  = &(p->p[p->np-1]);
    tu = d[p->np-1];
    for (n=0, i=p->np; i>0; i--, u=v, tu=tv, n++) {
      /*=================================================*/
      /* ON OLD POLYGON (P), U IS PREVIOUS VERTEX, V IS  */
      /* CURRENT VERTEX TV IS NEGATIVE IF VERTEX V IS IN */
      /*=================================================*/
      v  = &(p->p[n]);
      tv = d[n];
      if ((tu>0.0 && tv<0.0) ||
          (tu<0.0 && tv>0.0)) {
        /*=================================================*/
        /* EDGE CROSSES PLANE; ADD INTERSECTION POINT TO Q */
        /*=================================================*/
        t      = tu/(tu-tv);
        w      = &(q->p[q->np]);
        w->x   = u->x+t*(v->x-u->x);
        w->y   = u->y+t*(v->y-u->y);
        w->z   = u->z+t*(v->z-u->z);
        q->np++;
      }
      /*==============================*/
      /* VERTEX V IS IN, COPY IT TO Q */
      /*==============================*/
      if (tv>=0.0) {
        w    = &(q->p[q->np]);
        w->x = v->x;
        w->y = v->y;
        w->z = v->z;
        q->np++;
      }
    }
    if (q->np < 3) q->np = 0;
  }
}

void VF_ClipToPolyPlane(Poly *p, Poly *q, Poly *poly)
{
  int    i, n;
  double t, tu, tv;
  Vector *u, *v, *w;
  VFtopology* topology=VF_CurrentTopology();

  q->np     = 0;
  q->d      = p->d;
  q->normal = p->normal;
  if (topology->geom==VF_2Dplanar) {
    /*===========*/
    /* U=VERT[0] */
    /* V=VERT[1] */
    /*===========*/
    u  = &(p->p[0]);
    tu = V3_Dot(u, &poly->normal)+poly->d;
    if (fabs(tu)<topology->spatial_tol) tu=0.0;
    v  = &(p->p[1]);
    tv = V3_Dot(v, &poly->normal)+poly->d;
    if (fabs(tv)<topology->spatial_tol) tv=0.0;
    if (tu>=0.0) {
      /*==============================*/
      /* VERTEX U IS IN: COPY IT TO Q */
      /*==============================*/
      w    = &(q->p[q->np]);
      w->x = u->x;
      w->y = u->y;
      w->z = u->z;
      q->np++;
    }
    if ((tu>0.0 && tv<0.0) ||
        (tu<0.0 && tv>0.0)) {
      /*=================================================*/
      /* EDGE CROSSES PLANE: ADD INTERSECTION POINT TO Q */
      /*=================================================*/
      t    = tu/(tu-tv);
      w    = &(q->p[q->np]);
      w->x = u->x+t*(v->x-u->x);
      w->y = u->y+t*(v->y-u->y);
      w->z = u->z+t*(v->z-u->z);
      q->np++;
    }
    if (tv>=0.0) {
      /*==============================*/
      /* VERTEX V IS IN: COPY IT TO Q */
      /*==============================*/
      w    = &(q->p[q->np]);
      w->x = v->x;
      w->y = v->y;
      w->z = v->z;
      q->np++;
    }
    if (q->np<2) q->np=0;
    if (q->np>2) {
      printf("VF_ClipToPolyPlane():  q->np = %d\n",q->np);
      VF_PrintPoly(p,   "Poly p");
      VF_PrintPoly(q,   "Poly q");
      VF_PrintPoly(poly,"Plane Poly");
      VF_Exit(1);
    }
  } else {
    /*========================*/
    /* START WITH U=VERT[N-1] */
    /*========================*/
    u  = &(p->p[p->np-1]);
    tu = V3_Dot(u, &poly->normal)+poly->d;
    if (fabs(tu)<topology->spatial_tol) tu=0.0;
    for (n=0, i=p->np; i>0; i--, u=v, tu=tv, n++) {
      /*=================================================*/
      /* ON OLD POLYGON (P), U IS PREVIOUS VERTEX, V IS  */
      /* CURRENT VERTEX TV IS POSITIVE IF VERTEX V IS IN */
      /*=================================================*/
      v  = &(p->p[n]);
      tv = V3_Dot(v, &poly->normal)+poly->d;
      if (fabs(tv)<topology->spatial_tol) tv=0.0;
      if ((tu>0.0 && tv<0.0) ||
          (tu<0.0 && tv>0.0)) {
        /*=================================================*/
        /* EDGE CROSSES PLANE: ADD INTERSECTION POINT TO Q */
        /*=================================================*/
        t    = tu/(tu-tv);
        w    = &(q->p[q->np]);
        w->x = u->x+t*(v->x-u->x);
        w->y = u->y+t*(v->y-u->y);
        w->z = u->z+t*(v->z-u->z);
        q->np++;
      }
      if (tv>=0.0) {
        /*==============================*/
        /* VERTEX V IS IN: COPY IT TO Q */
        /*==============================*/
        w    = &(q->p[q->np]);
        w->x = v->x;
        w->y = v->y;
        w->z = v->z;
        q->np++;
      }
    }
    if (q->np<3) q->np=0;
    if (q->np>5) {
      printf("VF_ClipToPolyPlane():  q->np = %d\n",q->np);
      VF_PrintPoly(p,   "Poly p");
      VF_PrintPoly(q,   "Poly q");
      VF_PrintPoly(poly,"Plane Poly");
      VF_Exit(1);
    }
  }
}

int VF_ClipBoxToFrustum(ViewPort *view, Box *bounds, int face)
{
  int    i;
  double d, xMax, yMax, zMax;
  Vector *normal;
  for (i=1; i<5; i++) {
    d      =   view->frustum_planes1[i].d;
    normal = &(view->frustum_planes1[i].normal);
    xMax   = MAX( bounds->xmin * normal->x, bounds->xmax * normal->x );
    yMax   = MAX( bounds->ymin * normal->y, bounds->ymax * normal->y );
    zMax   = MAX( bounds->zmin * normal->z, bounds->zmax * normal->z );
    if ( xMax + yMax + zMax + d < 0.0 ) {
      return VF_FRUSTUM_CLIP_OUT;
    }
  }
  return VF_FRUSTUM_CLIP_PARTIAL;
}
