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
@(#)    $RCSfile: VF_FacetAndPolyUtils.c,v $
@(#)    $Revision: 1.4 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_FacetAndPolyUtils.c,v $
@(#)
@(#)    DESCRIPTION:  Various utilities to work on surface geometries.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "vf.h"

/*=============================*/
/* COMPUTE THE AREA OF A PATCH */ 
/*=============================*/
double VF_FacetArea(Facet *facet)
{
  Poly   poly;
  VFtopology* topology=VF_CurrentTopology();

  VF_FacetToPoly(facet, &poly);
  return VF_SurfArea(&poly,topology->geom);
}

double VF_SurfArea(Poly *poly, int dim)
{
  int    i, j;
  double area, jacobian, d0, d1;
  Point  d;
  static double wt[] = { 1.0, 1.0};
  static double gp[] = {-0.577350269189626, 0.577350269189626};

  switch (dim) {
  case VF_2Daxisym:
    V3_Sub(&(poly->p[1]),&(poly->p[0]),&d);
    d0   = sqrt(poly->p[0].x*poly->p[0].x+poly->p[0].z*poly->p[0].z);
    d1   = sqrt(poly->p[1].x*poly->p[1].x+poly->p[1].z*poly->p[1].z);
    area = M_PI*V3_Length(&d)*(d0+d1);
    break;
  case VF_2Dplanar:
    area = VF_PolyArea(poly);
    break;
  case VF_3D:
    if (poly->np == 4) {
      for (area=0.0, i=0; i<2; i++) {
        for (j=0; j<2; j++) {
          jacobian = VF_GetPolyJacobian(poly, gp[i],gp[j]);
          area += wt[i]*wt[j]*jacobian;
        }
      }
    } else {
      area = VF_PolyArea(poly);
    }
    break;
  }
  return area;
}

/*================================================*/
/* CHECK IF A PATCH IS PLANAR: RETURN 0=YES, 1=NO */
/*================================================*/
int VF_IsFacetPlanar(Facet *facet)
{
  Poly    poly;
    
  VF_FacetToPoly(facet, &poly);
  return (VF_IsPolyPlanar(&poly));
}

/*=============================================*/
/* COMPUTE THE AREA OF A CONVEX PLANAR POLYGON */
/*=============================================*/
double VF_PolyArea(Poly *poly)
{
  int    i, ii;
  double area;
  Point  a, c, d1, d2;

  V3_Zero(&a);
  switch (poly->np) {
  case 2:
    V3_Sub(&(poly->p[1]),&(poly->p[0]),&d1);
    area = V3_Length(&d1);
    break;
  default:
    for (i=1, ii=2; ii<poly->np; i++, ii++) {
      V3_Sub(&(poly->p[ i]),&(poly->p[0]),&d1);
      V3_Sub(&(poly->p[ii]),&(poly->p[0]),&d2);
      V3_Cross(&d1, &d2, &c);
      V3_Add(&a, &c, &a);
    }
    area = 0.5*V3_Length(&a);
    break;
  }
  return area;
}

/*===========================================================*/
/* COMPUTE THE AREA OF A CONVEX, POSSIBLY NONPLANAR, POLYGON */
/*===========================================================*/
double VF_GenericPolyArea(Poly *poly)
{
  int    i, ii, n, nn;
  double area;
  Poly   p;
  Point  a, c, d1, d2;

  V3_Zero(&a);
  switch (poly->np) {
  case 2:
    V3_Sub(&(poly->p[1]),&(poly->p[0]),&d1);
    area = V3_Length(&d1);
    break;
  case 3:
    for (i=1, ii=2; ii<3; i++, ii++) {
      V3_Sub(&(poly->p[ i]),&(poly->p[0]),&d1);
      V3_Sub(&(poly->p[ii]),&(poly->p[0]),&d2);
      V3_Cross(&d1, &d2, &c);
      V3_Add(&a, &c, &a);
    }
    area = 0.5*V3_Length(&a);
    break;
  case 4:
    if (VF_IsPolyPlanar(poly)) {
      for (i=1, ii=2; ii<poly->np; i++, ii++) {
         V3_Sub(&(poly->p[ i]),&(poly->p[0]),&d1);
         V3_Sub(&(poly->p[ii]),&(poly->p[0]),&d2);
         V3_Cross(&d1, &d2, &c);
         V3_Add(&a, &c, &a);
      }
      area = 0.5*V3_Length(&a);
    } else {
      area = 0.0;
      VF_PolyCenter(poly, &(p.p[0]));
      for (n=0; nn=1, n<4; n++, nn=(n+1)%4) {
        p.p[1] = poly->p[n];
        p.p[2] = poly->p[nn];
        for (i=1, ii=2; ii<3; i++, ii++) {
          V3_Sub(&(poly->p[ i]),&(poly->p[0]),&d1);
          V3_Sub(&(poly->p[ii]),&(poly->p[0]),&d2);
          V3_Cross(&d1, &d2, &c);
          V3_Add(&a, &c, &a);
        }
        area = 0.5*V3_Length(&a);
      }
    }
    break;
  }
  return area;
}

/*===============================================*/
/* COMPUTE THE NORMAL OF A CONVEX PLANAR POLYGON */
/*===============================================*/
void VF_PolyNormal(Poly *poly, Vector *normal)
{
  double tmp;
  Vector d1, d2;
    
  switch (poly->np) {
  case 2:
    V3_Sub(&(poly->p[1]), &(poly->p[0]), &d1);
    normal->x =  d1.y;
    normal->y = -d1.x;
    normal->z = 0.0;
    break;
  case 3:
    V3_Sub(&(poly->p[2]), &(poly->p[1]), &d1);
    V3_Sub(&(poly->p[0]), &(poly->p[1]), &d2);
    V3_Cross(&d1, &d2, normal);
    break;
  default:
    V3_Sub(&(poly->p[2]), &(poly->p[0]), &d1);
    V3_Sub(&(poly->p[3]), &(poly->p[1]), &d2);
    V3_Cross(&d1, &d2, normal);
    break;
  }
  V3_Normalize(normal,tmp);
}

void VF_PolyNormal_Aux(Poly *poly, Vector *normal, int face)
{
  double tmp;
  Vector d1;
    
  if (poly->np==2) {
    V3_Sub(&(poly->p[1]), &(poly->p[0]), &d1);
    switch (face) {
    case 0:
      normal->x =  d1.z;
      normal->y =  0.0;
      normal->z = -d1.x;
      break;
    case 1:
    case 2:
      normal->x =  0.0;
      normal->y = -d1.z;
      normal->z =  d1.y;
      break;
    case 3:
    case 4:
      normal->x =  d1.y;
      normal->y = -d1.x;
      normal->z =  0.0;
      break;
    }
    V3_Normalize(normal,tmp);
  }
}

/*=============================================================*/
/* COMPUTE THE NORMAL OF A CONVEX, POSSIBLY NONPLANAR, POLYGON */
/*=============================================================*/
void VF_GenericPolyNormal(Poly *poly, Vector *normal)
{
  int    i, ii;
  double tmp;
  Vector d1, d2, nn, center;
    
  if (poly->np==2) {
    V3_Sub(&(poly->p[1]), &(poly->p[0]), &d1);
    normal->x = d1.y;
    normal->y = -d1.x;
    normal->z = 0.0;
  } else {
    if (poly->np==3 || VF_IsPolyPlanar(poly)) {
      V3_Sub(&(poly->p[2]), &(poly->p[1]), &d1);
      V3_Sub(&(poly->p[0]), &(poly->p[1]), &d2);
      V3_Cross(&d1, &d2, normal);
    } else {
      normal->x = 0.0;
      normal->y = 0.0;
      normal->z = 0.0;
      VF_PolyCenter(poly,&center);
      for (i=0, ii=1; i<poly->np; i++, ii=(i+1)%4) {
        V3_Sub(&(poly->p[ i]), &center, &d1);
        V3_Sub(&(poly->p[ii]), &center, &d2);
        V3_Cross(&d1, &d2, &nn);
        V3_Add(&nn, normal, normal);
      }
    }
  }
  V3_Normalize(normal,tmp);
}

/*==================================================*/
/* CHECK IF A POLYGON IS PLANAR: RETURN 0=YES, 1=NO */
/*==================================================*/
int VF_IsPolyPlanar(Poly *poly)
{
  int    n;
  double d, value, tol=1.0e-8;
  Vector d1, d2, normal;
    
  if (poly->np!=4) return(1);
  V3_Sub(&(poly->p[2]), &(poly->p[1]), &d1);
  V3_Sub(&(poly->p[0]), &(poly->p[1]), &d2);
  V3_Cross(&d1, &d2, &normal);
  V3_Normalize(&normal, d);
  d = -V3_Dot(&normal, &(poly->p[0]));
  for (n=0; n<poly->np; n++) {
    value = d+V3_Dot(&(poly->p[n]), &normal);
    if (fabs(value)>tol) return(0);
  }
  return(1);
}

/*=================================*/
/* COMPUTE THE CENTER OF A POLYGON */
/*=================================*/
void VF_PolyCenter(Poly *poly, Vector *center)
{
  int    i;
  double x=0.0, y=0.0, z=0.0, tmp;

  tmp = 1.0/(double)(poly->np);
  for (i=0; i<poly->np; i++) {
    x += poly->p[i].x;
    y += poly->p[i].y;
    z += poly->p[i].z;
  }
  center->x = x*tmp;
  center->y = y*tmp;
  center->z = z*tmp;
}

void VF_FacetToPoly(Facet *facet, Poly *poly)
{
  int    n, nn, mirror, rot=0;
  int    vertex, vertex0, vertex1;
  int    quad_map[] = {0,1,1,0};
  int    rot_map[] = {0,0,1,1};
  double xm=1.0, ym=1.0, zm=1.0;
  double theta[2], costheta, sintheta;
  VFtopology *topology=VF_CurrentTopology();
    
  poly->np        = facet->num_vertices;
  poly->mask      = facet->mask;
  poly->normal    = facet->normal;
  poly->d         = facet->d;
  poly->facet     = facet->index;
  poly->sub_facet = -1;
  if (facet->sector) {
    if ((facet->sector&0x000F0000)>>16) {
      xm = -1.0;
    } else {
      xm = 1.0;
    }
    if ((facet->sector&0x00F00000)>>20) {
      ym = -1.0;
    } else {
      ym = 1.0;
    }
    if ((facet->sector&0x0F000000)>>24) {
      zm = -1.0;
    } else {
      zm = 1.0;
    }
    rot = facet->sector&0xFFFF;
  }
  switch (topology->geom) {
  case VF_2Daxisym:
    theta[0] = topology->theta*(double)(rot);
    theta[1] = topology->theta*(double)(rot+1);
    if (rot+1 == topology->nrotations) {
      theta[1] = 0.0;
    }
    for (n=0; n<facet->num_vertices; n++) {
      if (facet->num_vertices==4) {
        vertex       = facet->vertex_list[quad_map[n]];
        poly->p[n].x = xm*cos(theta[rot_map[n]])*topology->x[vertex];
        poly->p[n].y = ym*topology->y[vertex];
        poly->p[n].z = zm*sin(theta[rot_map[n]])*topology->x[vertex];
      } else {
        if (n<2) {
          vertex       = facet->vertex_list[n];
          poly->p[n].x = xm*cos(theta[0])*topology->x[vertex];
          poly->p[n].y = ym*topology->y[vertex];
          poly->p[n].z = zm*sin(theta[0])*topology->x[vertex];
        } else {
          vertex0 = facet->vertex_list[0];
          vertex1 = facet->vertex_list[1];
          if (fabs((topology->x[vertex0])) < topology->spatial_tol) {
            vertex = vertex1;
          }
          if (fabs((topology->x[vertex1])) < topology->spatial_tol) {
            vertex = vertex0;
          }
          poly->p[n].x = xm*cos(theta[1])*topology->x[vertex];
          poly->p[n].y = ym*topology->y[vertex];
          poly->p[n].z = zm*sin(theta[1])*topology->x[vertex];
        }
      }
    }
    break;
  case VF_2Dplanar:
    if (facet->sector) {
      for (n=0; n<facet->num_vertices; n++) {
        vertex       = facet->vertex_list[n];
        poly->p[n].x = xm*topology->x[vertex];
        poly->p[n].y = ym*topology->y[vertex];
        poly->p[n].z = 0.0;
      }
    } else {
      for (n=0; n<facet->num_vertices; n++) {
        vertex       = facet->vertex_list[n];
        poly->p[n].x = topology->x[vertex];
        poly->p[n].y = topology->y[vertex];
        poly->p[n].z = 0.0;
      }
    }
    break;
  case VF_3D:
    if (facet->sector) {
      costheta = 1.0;
      sintheta = 0.0;
      if (rot) {
        theta[0] = topology->theta*(double)(rot);
        costheta = cos(theta[0]);
        sintheta = sin(theta[0]);
      }
      mirror = (xm*ym*zm)>=0.0?1:0;
      for (n=0; n<facet->num_vertices; n++) {
        if (mirror) {
          nn = n;
        } else {
          nn = facet->num_vertices-n-1;
        }
        vertex = facet->vertex_list[n];
        if (rot) {
          poly->p[nn].x = xm*(costheta*topology->x[vertex]+
                              sintheta*topology->z[vertex]);
          poly->p[nn].y = ym*topology->y[vertex];
          poly->p[nn].z = zm*(-sintheta*topology->x[vertex]+
                               costheta*topology->z[vertex]);
        } else {
          poly->p[nn].x = xm*topology->x[vertex];
          poly->p[nn].y = ym*topology->y[vertex];
          poly->p[nn].z = zm*topology->z[vertex];
        }
      }
    } else {
      for (n=0; n<facet->num_vertices; n++) {
        vertex       = facet->vertex_list[n];
        poly->p[n].x = topology->x[vertex];
        poly->p[n].y = topology->y[vertex];
        poly->p[n].z = topology->z[vertex];
      }
    }
    break;
  }
}

void VF_FacetToPolyStack(Facet *facet, POLYstack *PolyStack, int *cnt)
{
  int   i, ii;
  Point center;
  Poly  *poly, p;
  VFtopology *topology=VF_CurrentTopology();

  switch (topology->geom) {
  case VF_2Daxisym:
    poly = &(PolyStack->poly[PolyStack->cnt]);
    VF_FacetToPoly(facet,poly);
    PolyStack->cnt++;
    break;
  case VF_2Dplanar:
    poly = &(PolyStack->poly[PolyStack->cnt]);
    VF_FacetToPoly(facet,poly);
    PolyStack->cnt++;
    break;
  case VF_3D:
    if (facet->mask&VF_FACET_MASK_PLANAR) {
      poly = &(PolyStack->poly[PolyStack->cnt]);
      VF_FacetToPoly(facet,poly);
      PolyStack->cnt++;
    } else {
      VF_FacetToPoly(facet,&p);
      VF_PolyCenter(&p,&center);
      for (i=0, ii=1; i<4; i++, ii=(i+1)%4) {
        poly       = &(PolyStack->poly[PolyStack->cnt]);
        poly->np   = 3;
        poly->sub_facet = i;
        poly->p[0] = center;
        poly->p[1] = p.p[i];
        poly->p[2] = p.p[ii];
        poly->mask = facet->mask;
        VF_PolyNormal(poly, &(poly->normal));
        poly->d = -V3_Dot(&(poly->normal), &(poly->p[0]));
        PolyStack->cnt++;
      }
    }
    break;
  }
  *cnt = PolyStack->cnt;
}

void VF_SubPoly(int i, int j, int ni, int nj, Poly *poly, Poly *subpoly)
{
    int    k;
    double L1, L2, u0, u1, v0, v1;
    double Ni = (double)ni;
    double Nj = (double)nj;
    
    switch (poly->np) {
    case 2:
        u0            = MAX(0.0,(double)i/Ni);
        u1            = MIN(1.0,(double)(i+1)/Nj);
        subpoly->np   = poly->np;
        subpoly->p[0] = VF_UVtoXYZ(u0, 0.0, poly);
        subpoly->p[1] = VF_UVtoXYZ(u1, 0.0, poly);
        break;
    case 3:
        k  = j/2;
        L1 = (double)i/Ni;
        L2 = 0.0;
        u0 = MAX(0.0,L1-(double)k/Ni);
        v0 = MIN(1.0,L2+(double)k/Nj);
        if (j%2) {
            u1            = MAX(0.0,u0-1.0/Ni);
            v1            = MIN(1.0,v0+1.0/Nj);
            subpoly->p[0] = VF_UVtoXYZ(u0, v0, poly);
            subpoly->p[1] = VF_UVtoXYZ(u0, v1, poly);
            subpoly->p[2] = VF_UVtoXYZ(u1, v1, poly);
        } else {
            u1            = MAX(0.0,u0+1.0/Ni);
            v1            = MIN(1.0,v0+1.0/Nj);
            subpoly->p[0] = VF_UVtoXYZ(u0, v0, poly);
            subpoly->p[1] = VF_UVtoXYZ(u1, v0, poly);
            subpoly->p[2] = VF_UVtoXYZ(u0, v1, poly);
        }
        subpoly->np = poly->np;
        break;
    case 4:
        u0            = MAX(0.0,(double)i/Ni);
        v0            = MIN(1.0,(double)j/Nj);
        u1            = MAX(0.0,(double)(i+1)/Ni);
        v1            = MIN(1.0,(double)(j+1)/Nj);
        subpoly->np   = poly->np;
        subpoly->p[0] = VF_UVtoXYZ(u0, v0, poly);
        subpoly->p[1] = VF_UVtoXYZ(u1, v0, poly);
        subpoly->p[2] = VF_UVtoXYZ(u1, v1, poly);
        subpoly->p[3] = VF_UVtoXYZ(u0, v1, poly);
        break;
    default:
        printf("SubPoly(): poly->np = %d not supported\n",poly->np);
        VF_Exit(1);
        break;
    }
}

void VF_PrintPoly(Poly *poly, char tag[])
{
  int n;
    
  printf("%s:  Npoints = %d\n",tag,poly->np);
  for (n=0; n<poly->np; n++) {
    printf("   Vertex %d = (%g, %g, %g)\n",
           n,poly->p[n].x,poly->p[n].y,poly->p[n].z);
  }
  printf("   normal  = (%g, %g, %g)\n",
         poly->normal.x,poly->normal.y,poly->normal.z);
  printf("   d       = %g\n",poly->d);
}

void VF_PrintVector(Vector *v, char tag[])
{
  printf("%s:\n",tag);
  printf("   x = %g\n",v->x);
  printf("   y = %g\n",v->y);
  printf("   z = %g\n",v->z);
}

int VF_SharedEdge(Facet *facet_i, Facet *facet_j, Poly *poly_i, Poly *poly_j)
{
  int    i, ii, j, jj, shared=0;
  int    ni, nj, vi1, vi2, vj1, vj2;
  double d;
  VFtopology *topology;
    
  topology = VF_CurrentTopology();
  ni   = poly_i->np;
  nj   = poly_j->np;
  for (ii=1, i=0; i<ni; i++, ii=(i+1)%ni) {
    vi1 = facet_i->vertex_list[i];
    vi2 = facet_i->vertex_list[(i+1)%ni];
    for (jj=1, j=0; j<nj; j++, jj=(j+1)%nj) {
      vj1 = facet_j->vertex_list[j];
      vj2 = facet_j->vertex_list[(j+1)%nj];
      if (vi1==vj2 && vi2==vj1) {
        shared = 1;
        break;
      } else {
        d = V3_DistanceBetween2Points(&(poly_i->p[i]), &(poly_j->p[jj]));
        if (d<topology->spatial_tol) {
          d = V3_DistanceBetween2Points(&(poly_i->p[ii]), &(poly_j->p[j]));
          if (d<topology->spatial_tol) {
            shared = 1;
            break;
          }
        }
      }
    }
    if (shared) break;
  }
  return(shared);
}

int VF_SharedPolyEdge(Poly *poly_i, Poly *poly_j)
{
  int    i1, i2, j1, j2, shared=0;
  double d, tmp, dot, value1, value2;
  Vector edge_i, edge_j, dir;
  Point  c;
  VFtopology *topology;
    
  topology = VF_CurrentTopology();
  /*======================================*/
  /* FIRST, CHECK IF ANY VERTEX OF POLY_J */
  /* IS THE SAME AS A VERTEX OF POLY_I    */
  /*======================================*/
  for (j1=0; j1<poly_j->np; j1++) {
    for (i1=0; i1<poly_i->np; i1++) {
      d = V3_DistanceBetween2Points(&(poly_i->p[i1]),&(poly_j->p[j1]));
      if (d < topology->spatial_tol) {
        shared = 1;
        break;
      }
    }
  }
  if (shared) return(shared);
  /*====================================================================*/
  /* SECOND, CHECK IF ANY EDGE OF POLY_J IS IN THE SAME PLANE AS POLY_I */
  /*====================================================================*/
  VF_PolyCenter(poly_i,&c);
  for (j2=1, j1=0; j1<poly_j->np; j1++, j2=(j1+1)%poly_j->np) {
    value1 = V3_Dot(&(poly_j->p[j1]), &(poly_i->normal)) + poly_i->d;
    value2 = V3_Dot(&(poly_j->p[j2]), &(poly_i->normal)) + poly_i->d;
    if ((value1 < topology->spatial_tol && value1 < topology->spatial_tol) &&
        (value2 < topology->spatial_tol && value2 < topology->spatial_tol)) {
      /*===============================================*/
      /* THIS EDGE IS IN THE SAME PLANE AS POLY_I, NOW */
      /* CHECK IF IT IS SHARED WITH ANY EDGE ON POLY_I */
      /*===============================================*/
      V3_Sub(&(poly_j->p[j1]),&(poly_j->p[j2]),&edge_j);
      V3_Normalize(&edge_j,tmp);
      for (i2=1, i1=0; i1<poly_i->np; i1++, i2=(i1+1)%poly_i->np) {
        V3_Sub(&(poly_i->p[i2]),&(poly_i->p[i1]),&edge_i);
        V3_Normalize(&edge_i,tmp);
        dot = V3_Dot(&edge_i, &edge_j);
        if (dot+topology->spatial_tol>=1.0) {
          /*====================*/
          /* EDGES ARE PARALLEL */
          /*====================*/
          V3_Sub(&(poly_i->p[i1]),&(poly_j->p[j1]),&dir);
          V3_Normalize(&dir,tmp);
          dot = V3_Dot(&edge_i, &dir);
          if ((dot > 0.0 && dot+topology->spatial_tol >=  1.0) || 
              (dot < 0.0 && dot-topology->spatial_tol <= -1.0)) {
            shared = 1;
            break;
          }
        }
      }
      if (shared) break;
    }
  }
  return(shared);
}

/*=====================================================*/
/* COMPUTE XYZ COORDINATES OF A POINT ON QUADRILATERAL */
/* GIVEN IT'S U,V COORDINATES USING BI-LINEAR MAPPING  */ 
/*=====================================================*/
void VF_UV_to_XYZ(Poly *q, double u, double v, Vector *p)
{
  int    i;
  double phi[4];
  
  p->x = 0.0;
  p->y = 0.0;
  p->z = 0.0;
  switch (q->np) {
  case 2:
    phi[0] = 1.0-u;
    phi[1] = u;
    break;
  case 3:
    phi[0] = u;
    phi[1] = v;
    phi[2] = 1.0-u-v;
    break;
  case 4:
    phi[0] = (1.0-u)*(1.0-v);
    phi[1] = (    u)*(1.0-v);
    phi[2] = (    u)*(    v);
    phi[3] = (1.0-u)*(    v);
    break;
  }
  for (i=0; i<q->np; i++) {
    p->x += phi[i] * q->p[i].x;
    p->y += phi[i] * q->p[i].y;
    p->z += phi[i] * q->p[i].z;
  }
}
/* U,V = [0,1] */

Point VF_UVtoXYZ(double u, double v, Poly *poly)
{
  int    i;
  double phi[4];
  Point  p;
    
  p.x = 0.0;
  p.y = 0.0;
  p.z = 0.0;
  switch (poly->np) {
  case 2:
    if (u<0.0 || u>1.0) {
      printf("L:(u,v) out of range (%g, %g)\n",u,v);
      VF_PrintPoly(poly,"Sampled Poly");
      VF_Exit(0);
    }
    phi[0] = 1.0-u;
    phi[1] = u;
    for (i=0; i<2; i++) {
      p.x += phi[i]*poly->p[i].x;
      p.y += phi[i]*poly->p[i].y;
    }
    break;
  case 3:
    if (u<0.0 || u>1.0 || v<0.0 || v>1.0) {
      printf("T:(u,v) out of range (%g, %g)\n",u,v);
      VF_PrintPoly(poly,"Sampled Poly");
      VF_Exit(0);
    }
    phi[0] = u;
    phi[1] = v;
    phi[2] = 1.0-u-v;
    for (i=0; i<3; i++) {
      p.x += phi[i]*poly->p[i].x;
      p.y += phi[i]*poly->p[i].y;
      p.z += phi[i]*poly->p[i].z;
    }
    break;
  case 4:
    if (u<0.0 || u>1.0 || v<0.0 || v>1.0) {
      printf("Q:(u,v) out of range (%g, %g)\n",u,v);
      VF_PrintPoly(poly,"Sampled Poly");
      VF_Exit(0);
    }
    phi[0] = (1.0-u)*(1.0-v);
    phi[1] = (    u)*(1.0-v);
    phi[2] = (    u)*(    v);
    phi[3] = (1.0-u)*(    v);
    for (i=0; i<4; i++) {
      p.x += phi[i]*poly->p[i].x;
      p.y += phi[i]*poly->p[i].y;
      p.z += phi[i]*poly->p[i].z;
    }
    break;
  default:
    printf("VF_UVtoXYZ():  poly->np = %d\n",poly->np);
    VF_Exit(1);
    break;
  }
  return (p);
}

double VF_GetPolyJacobian(Poly *poly, double u, double v)
{
    int    n;
    double dpdu[4], dpdv[4], phi[4];
    double dxdu=0.0, dydu=0.0, dzdu=0.0;
    double dxdv=0.0, dydv=0.0, dzdv=0.0;
    double e1_dot_e1, e2_dot_e2, e1_dot_e2;

    /*==========================*/
    /* EVALUATE SHAPE FUNCTIONS */
    /*==========================*/
    phi[0] = 0.25*(1.0-u)*(1.0-v);
    phi[1] = 0.25*(1.0+u)*(1.0-v);
    phi[2] = 0.25*(1.0+u)*(1.0+v);
    phi[3] = 0.25*(1.0-u)*(1.0+v);
    /*===============================================*/
    /* DERIVATIVES OF SHAPE FUNCTIONS W.R.T. U AND V */
    /*===============================================*/
    dpdu[0] = -0.25*(1.0-v);
    dpdu[1] = -dpdu[0];
    dpdu[2] =  0.25*(1.0+v);
    dpdu[3] = -dpdu[2];
    dpdv[0] = -0.25*(1.0-u);
    dpdv[1] = -0.25*(1.0+u);
    dpdv[2] = -dpdv[1];
    dpdv[3] = -dpdv[0];
    /*=====================================*/
    /* TREAT TRIANGLE AS DEGENERATIVE QUAD */
    /*=====================================*/
    if (poly->np==3) {
        phi[2]  += phi[3];
        dpdu[2] += dpdu[3];
        dpdv[2] += dpdv[3];
        phi[3]   = 0.0;
        dpdu[3]  = 0.0;
        dpdv[3]  = 0.0;
    }
    /*===========================================*/
    /* DERIVATIVES OF X, Y, AND Z W.R.T. U AND V */
    /*===========================================*/
    for (n=0; n<poly->np; n++) {
        dxdu += dpdu[n]*poly->p[n].x;
        dydu += dpdu[n]*poly->p[n].y;
        dzdu += dpdu[n]*poly->p[n].z;
        dxdv += dpdv[n]*poly->p[n].x;
        dydv += dpdv[n]*poly->p[n].y;
        dzdv += dpdv[n]*poly->p[n].z;
    }
    /*=====================================*/
    /* COMPUTE DETERMINANT OF THE JACOBIAN */
    /*=====================================*/
    e1_dot_e1 = dxdu*dxdu+dydu*dydu+dzdu*dzdu;
    e2_dot_e2 = dxdv*dxdv+dydv*dydv+dzdv*dzdv;
    e1_dot_e2 = dxdu*dxdv+dydu*dydv+dzdu*dzdv;
    return (sqrt(e1_dot_e1*e2_dot_e2-e1_dot_e2*e1_dot_e2));
}

void VF_GetShapeFunction(Poly *poly, double u, double v, double phi[], double *jacobian)
{
    int    n;
    double dpdu[4], dpdv[4];
    double dxdu=0.0, dydu=0.0, dzdu=0.0;
    double dxdv=0.0, dydv=0.0, dzdv=0.0;
    double e1_dot_e1, e2_dot_e2, e1_dot_e2;

    /*==========================*/
    /* EVALUATE SHAPE FUNCTIONS */
    /*==========================*/
    phi[0] = 0.25*(1.0-u)*(1.0-v);
    phi[1] = 0.25*(1.0+u)*(1.0-v);
    phi[2] = 0.25*(1.0+u)*(1.0+v);
    phi[3] = 0.25*(1.0-u)*(1.0+v);
    /*===============================================*/
    /* DERIVATIVES OF SHAPE FUNCTIONS W.R.T. U AND V */
    /*===============================================*/
    dpdu[0] = -0.25*(1.0-v);
    dpdu[1] = -dpdu[0];
    dpdu[2] =  0.25*(1.0+v);
    dpdu[3] = -dpdu[2];
    dpdv[0] = -0.25*(1.0-u);
    dpdv[1] = -0.25*(1.0+u);
    dpdv[2] = -dpdv[1];
    dpdv[3] = -dpdv[0];
    /*=====================================*/
    /* TREAT TRIANGLE AS DEGENERATIVE QUAD */
    /*=====================================*/
    if (poly->np==3) {
        phi[2]  += phi[3];
        dpdu[2] += dpdu[3];
        dpdv[2] += dpdv[3];
        phi[3]   = 0.0;
        dpdu[3]  = 0.0;
        dpdv[3]  = 0.0;
    }
    /*===========================================*/
    /* DERIVATIVES OF X, Y, AND Z W.R.T. U AND V */
    /*===========================================*/
    for (n=0; n<poly->np; n++) {
        dxdu += dpdu[n]*poly->p[n].x;
        dydu += dpdu[n]*poly->p[n].y;
        dzdu += dpdu[n]*poly->p[n].z;
        dxdv += dpdv[n]*poly->p[n].x;
        dydv += dpdv[n]*poly->p[n].y;
        dzdv += dpdv[n]*poly->p[n].z;
    }
    /*=====================================*/
    /* COMPUTE DETERMINANT OF THE JACOBIAN */
    /*=====================================*/
    e1_dot_e1 = dxdu*dxdu+dydu*dydu+dzdu*dzdu;
    e2_dot_e2 = dxdv*dxdv+dydv*dydv+dzdv*dzdv;
    e1_dot_e2 = dxdu*dxdv+dydu*dydv+dzdu*dzdv;
    *jacobian = sqrt(e1_dot_e1*e2_dot_e2-e1_dot_e2*e1_dot_e2);
}

void VF_TransformPoly(Poly *poly, Poly *p, Matrix *xform)
{
  int   i;
  Point *vertex;
    
  p->np = poly->np;
  for (i=0; i<p->np; i++) {
    vertex    = &(poly->p[i]);
    p->p[i].x = vertex->x * xform->element[0][0] +
                vertex->y * xform->element[0][1] +
                vertex->z * xform->element[0][2] +
                            xform->element[0][3];
    p->p[i].y = vertex->x * xform->element[1][0] +
                vertex->y * xform->element[1][1] +
                vertex->z * xform->element[1][2] +
                            xform->element[1][3];
    p->p[i].z = vertex->x * xform->element[2][0] +
                vertex->y * xform->element[2][1] +
                vertex->z * xform->element[2][2] +
                            xform->element[2][3];
  }
  VF_PolyNormal(p, &(p->normal));
  p->d = -V3_Dot(&(p->normal), &(p->p[0]));
}

double VF_QuadAspectRatio(Poly *poly)
{
    double aspect, u_len, v_len;
    
    u_len  = V3_DistanceBetween2Points(&(poly->p[0]), &(poly->p[1])) +
             V3_DistanceBetween2Points(&(poly->p[3]), &(poly->p[2]));
    v_len  = V3_DistanceBetween2Points(&(poly->p[1]), &(poly->p[2])) +
             V3_DistanceBetween2Points(&(poly->p[0]), &(poly->p[3]));
    aspect = u_len/v_len;
    return(aspect);
}

double VF_TriAspectRatio(Poly *poly, double *len1, double *len2, double *len3)
{
    double aspect, srms, sq_len1, sq_len2, sq_len3;
    
    *len1   = V3_DistanceBetween2Points(&(poly->p[0]), &(poly->p[1]));
    *len2   = V3_DistanceBetween2Points(&(poly->p[1]), &(poly->p[2]));
    *len3   = V3_DistanceBetween2Points(&(poly->p[2]), &(poly->p[0]));
    sq_len1 = (*len1)*(*len1);
    sq_len2 = (*len2)*(*len2);
    sq_len3 = (*len3)*(*len3);
    srms    = (sq_len1+sq_len2+sq_len3)/3.0;
    aspect  = srms/(2.30940108*VF_PolyArea(poly));
    return(aspect);
}
