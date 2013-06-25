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
@(#)    $RCSfile: VF_Visibility.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_Visibility.c,v $
@(#)
@(#)    DESCRIPTION:  Ray-polygon intersection utilities 
@(#)    (used for visibility tests).
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

double VF_Visibility(Poly *poly_i, Poly *poly_j, 
                     CandidateList *candidates, int nc)
{
  int    i, j, n, ns, cnt, nt, test_vertices=0;
  double tmax, tmp, n_dot_d;
  Ray    ray0,*ray;
  Point  samples_i[4096], samples_j[4096];
  Vector q;
  Adaptive    *adaptive;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  adaptive = &(enclosure->adaptive);
  if (nc==0) return(1.0);
  VF_SamplesUVtoXYZ(poly_i, &(adaptive->visibility.sampling), 
                    adaptive->visibility.sampleuv_i, samples_i);
  VF_SamplesUVtoXYZ(poly_j, &(adaptive->visibility.sampling), 
                    adaptive->visibility.sampleuv_j, samples_j);
  ns = adaptive->visibility.sampling.n;
  if (test_vertices) {
    cnt = nt = ns+poly_i->np*poly_j->np;
  } else {
    cnt = nt = ns;
  }
  /*=================================================*/
  /* TEST VISIBILITY OF RAYS BETWEEN SAMPLING POINTS */
  /*=================================================*/
  for (n=0; n<ns; n++) {
    ray      = &(adaptive->visibility.ray[n]);
    ray->O   = samples_i[n];
    ray->vis = 0;
    V3_Sub(&(samples_j[n]), &(ray->O), &(ray->D));
    V3_Normalize(&(ray->D),tmp);
    n_dot_d = V3_Dot(&(poly_j->normal), &(ray->D));
    if (fabs(n_dot_d)<VF_RAY_TEST) {
      cnt--;
    } else {
      V3_Sub(&(poly_j->p[0]), &(ray->O), &q);
      tmax = V3_Dot(&(poly_j->normal), &q)/n_dot_d;
      if (VF_RayPatchListTest(&(adaptive->visibility.ray[n]), 
                              n, tmax, candidates, nc)) cnt--;
    }
  }
  if (test_vertices) {
    /*===============================================================*/
    /* ALSO TEST RAYS BETWEEN ALL VERTEX COMBINATIONS FOR VISIBILITY */
    /*===============================================================*/
    for (i=0; i<poly_i->np; i++) {
      for (j=0; j<poly_j->np; j++) {
        ray0.vis = 0;
        ray0.O   = poly_i->p[i];
        V3_Sub(&(poly_j->p[j]), &(ray0.O), &(ray0.D));
        V3_Normalize(&(ray0.D),tmp);
        n_dot_d = V3_Dot(&(poly_j->normal), &(ray0.D));
        if (fabs(n_dot_d)<VF_RAY_TEST) {
          cnt--;
        } else {
          V3_Sub(&(poly_j->p[0]), &(ray0.O), &q);
          tmax = V3_Dot(&(poly_j->normal), &q)/n_dot_d;
          if (VF_RayPatchListTest(&ray0, n, tmax, candidates, nc)) {
            cnt--;
          }
        }
      }
    }
  }
  return ((double)(cnt)/(double)(nt));
}

/*=======================================*/
/* CHECK RAY <-> PATCH LIST INTERSECTION */
/*   RETURN 0 - NO INTERSECTION          */
/*   RETURN 1 - HAVE INTERSECTION        */
/*=======================================*/
int 
VF_RayPatchListTest(Ray *ray, int sample_num, double tmax,
                    CandidateList *candidates, int ncandidates)
{
  int    n, nn, intersect=0;
  static int cache_size=10;
  static int cache_cnt;
  static int cache_loc;
  static int cache[10];

  if (sample_num==0) {
    /*====================================*/
    /* IF THIS IS THE FIRST SAMPLE POINT, */
    /* INITIALIZE THE SHADOW CACHE        */
    /*====================================*/
    cache_cnt = 0;
    cache_loc = 0;
    for (n=0; n<ncandidates; n++) {
      candidates[n].cached = 0;
    }
  } else {
    /*==========================================*/
    /* IF NOT THE FIRST SAMPLE POINT, CHECK ALL */
    /* THE CANDIDATES IN THE SHADOW CACHE FIRST */
    /*==========================================*/
    nn = cache_loc;
    for (n=0; n<cache_cnt; n++) {
      intersect = VF_RayIntersect(ray, tmax, &(candidates[cache[nn]].poly));
      if (intersect) {
        cache_loc = nn;
        break;
      }
      nn--;
      if (nn<0) nn = cache_cnt-1;
    }
  }
  if (!intersect) {
    /*=======================================*/
    /* IF NO INTERSECTIONS FROM THE SHADOW   */
    /* CACHE, CHECK ALL THE OTHER CANDIDATES */
    /*=======================================*/
    for (n=0; n<ncandidates; n++) {
      if (!candidates[n].cached) {
        intersect = VF_RayIntersect(ray, tmax, &(candidates[n].poly));
        if (intersect) {
          /*================================*/
          /* IF INTERSECTION FOUND, ADD THE */
          /* CANDIDATE TO THE SHADOW CACHE  */
          /*================================*/
          if (cache_cnt<cache_size) {
            /*========================================*/
            /* CACHE NOT FULL, JUST ADD THE CANDIDATE */
            /*========================================*/
            cache_loc            = cache_cnt;
            cache[cache_loc]     = n;
            candidates[n].cached = 1;
            cache_cnt++;
          } else {
            /*============================================*/
            /* CACHE IS FULL, RECYCLE THE OLDEST POSITION */
            /*============================================*/
            cache_loc++;
            if (cache_loc>=cache_size) cache_loc = 0;
            candidates[cache[cache_loc]].cached = 0;
            cache[cache_loc]     = n;
            candidates[n].cached = 1;
          }
          break;
        }
      }
    }
  }
  return (intersect);
}

/*====================================*/
/* CHECK RAY <-> POLYGON INTERSECTION */
/*   RETURN 0 - NO INTERSECTION       */
/*   RETURN 1 - HAVE INTERSECTION     */
/*====================================*/
int VF_RayIntersect(Ray *ray, double tmax, Poly *poly)
{
  double t, n_dot_d;
  Vector q;   
    
  n_dot_d = V3_Dot(&(poly->normal), &ray->D);
  if (fabs(n_dot_d)<VF_RAY_TEST) return 0;
  V3_Sub(&(poly->p[0]), &ray->O, &q);
  t = V3_Dot(&(poly->normal), &q)/n_dot_d;
  if (t<=0.0 || t>=tmax) return 0;
  if (poly->np!=2) {
    return (VF_InsidePoly(ray, t, poly));
  } else {
    return (VF_InsideLine(ray, t, poly));
  }
}

/*======================================================*/
/* CHECK RAY - LINE INTERSECTION IS INSIDE LINE SEGMENT */
/*   RETURN 0 - NO                                      */
/*   RETURN 1 - YES                                     */
/*======================================================*/
int OldInsideLine(Ray *ray, double t, Poly *poly)
{
    double d, d0, d1;
    Vector normal;
    
    normal.x =  ray->D.y;
    normal.y = -ray->D.x;
    normal.z = 0.0;
    d        = -V2_Dot(&(ray->O)    , &normal);
    d0       =  V2_Dot(&(poly->p[0]), &normal) + d;
    d1       =  V2_Dot(&(poly->p[1]), &normal) + d;
    /*=============================================================*/
    /* NO INTERSECTION IF THE RAY PASSES THRU ONE OF THE ENDPOINTS */
    /*=============================================================*/
    if (d0==0.0 || d1==0.0) return 0;
    /*===============================================================*/
    /* NO INTERSECTION IF THE RAY IS COLLINEAR WITH THE LINE SEGMENT */
    /*===============================================================*/
    if (d0==d1) return 0;
    /*==============================================*/
    /* NO INTERSECTION IF BOTH ENDPOINTS OF THE     */
    /* LINE SEGMENT ARE ON THE SAME SIDE OF THE RAY */
    /*==============================================*/
    if (ZSGN(d0)==ZSGN(d1)) return 0;
    /*================================================*/
    /* OTHERWISE, THE RAY INTERSECTS THE LINE SEGMENT */
    /*================================================*/
    return 1;
}

/*======================================================*/
/* CHECK RAY - LINE INTERSECTION IS INSIDE LINE SEGMENT */
/*   RETURN 0 - NO                                      */
/*   RETURN 1 - YES                                     */
/*======================================================*/
int VF_InsideLine(Ray *ray, double t, Poly *poly)
{
  double d, d0, d1, tol;
  Vector normal;
  VFtopology *topology=VF_CurrentTopology();
    
  tol      = topology->spatial_tol;
  normal.x =  ray->D.y;
  normal.y = -ray->D.x;
  normal.z = 0.0;
  d        = -V2_Dot(&(ray->O)    , &normal);
  d0       =  V2_Dot(&(poly->p[0]), &normal) + d;
  d1       =  V2_Dot(&(poly->p[1]), &normal) + d;
  /*=============================================================*/
  /* NO INTERSECTION IF THE RAY PASSES THRU ONE OF THE ENDPOINTS */
  /*=============================================================*/
  if (fabs(d0)<tol || fabs(d1)<tol) return 0;
  /*===============================================================*/
  /* NO INTERSECTION IF THE RAY IS COLLINEAR WITH THE LINE SEGMENT */
  /*===============================================================*/
  if (fabs(d0-d1)<2.0*tol) return 0;
  /*==============================================*/
  /* NO INTERSECTION IF BOTH ENDPOINTS OF THE     */
  /* LINE SEGMENT ARE ON THE SAME SIDE OF THE RAY */
  /*==============================================*/
  if (ZSGN(d0)==ZSGN(d1)) return 0;
  /*================================================*/
  /* OTHERWISE, THE RAY INTERSECTS THE LINE SEGMENT */
  /*================================================*/
  return 1;
}

/*==================================================*/
/* CHECK RAY - PLANE INTERSECTION IS INSIDE POLYGON */
/*   RETURN 0 - NO                                  */
/*   RETURN 1 - YES                                 */
/*==================================================*/
int VF_InsidePoly(Ray *ray, double t, Poly *poly)
{
  int    i0, i=2, inter=0;
  double u0, u1, u2, v0, v1, v2, alpha, beta;
  Point  P, mag;
    
  P.x = ray->O.x + ray->D.x*t;
  P.y = ray->O.y + ray->D.y*t;
  P.z = ray->O.z + ray->D.z*t;
  V3_Mul(&(poly->normal),&(poly->normal),&mag);
  if (mag.x>=mag.y && mag.x>=mag.z) {
    i0 = 0;
    u0 = P.y - poly->p[0].y;
    v0 = P.z - poly->p[0].z;
  }
  else if (mag.y>=mag.x && mag.y>=mag.z) { 
    i0 = 1;
    u0 = P.z - poly->p[0].z;
    v0 = P.x - poly->p[0].x;
  }
  else if (mag.z>=mag.x && mag.z>=mag.y) {
    i0 = 2;
    u0 = P.x - poly->p[0].x;
    v0 = P.y - poly->p[0].y;
  }
  /*==========================================*/
  /* THE POLYGON IS VIEWED AS (N-2) TRIANGLES */
  /*==========================================*/
  do {
    /*==========================================*/
    /* BEGIN BARYCENTRIC INTERSECTION ALGORITHM */
    /*==========================================*/
    switch (i0) {
    case 0:
      u1 = poly->p[i-1].y - poly->p[0].y;
      u2 = poly->p[i  ].y - poly->p[0].y;
      v1 = poly->p[i-1].z - poly->p[0].z;
      v2 = poly->p[i  ].z - poly->p[0].z;
      break;
    case 1:
      u1 = poly->p[i-1].z - poly->p[0].z;
      u2 = poly->p[i  ].z - poly->p[0].z;
      v1 = poly->p[i-1].x - poly->p[0].x;
      v2 = poly->p[i  ].x - poly->p[0].x;
      break;
    case 2:
      u1 = poly->p[i-1].x - poly->p[0].x;
      u2 = poly->p[i  ].x - poly->p[0].x;
      v1 = poly->p[i-1].y - poly->p[0].y;
      v2 = poly->p[i  ].y - poly->p[0].y;
      break;
    }
    /*===============================================*/
    /* CALCULATE AND COMPARE BARYCENTRIC COORDINATES */
    /*===============================================*/
    if (u1!=0.0)    {  /* common case */
      beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
      if ((beta>=0.0) && (beta<=1.0)) {
        alpha = (u0 - beta*u2)/u1;
        inter = ((alpha>=0.0) && ((alpha+beta)<=1.0));
      }
    } else {           /* uncommon case */
      beta = u0/u2;
      if ((beta>=0.0) && (beta<=1.0)) {
        alpha = (v0 - beta*v2)/v1;
        inter = ((alpha>=0.0) && ((alpha+beta)<=1.0));
      }
    }
  } while ((!inter) && (++i<poly->np));
  if (inter) ray->vis = 1;
  return (inter);
}

