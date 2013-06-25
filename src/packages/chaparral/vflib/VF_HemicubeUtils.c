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
@(#)    $RCSfile: VF_HemicubeUtils.c,v $
@(#)    $Revision: 1.5 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_HemicubeUtils.c,v $
@(#)
@(#)    DESCRIPTION:  Initialize the hemicube delta-viewfactors.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "vf.h"

static int VF_HemicubeJitter=1;

void 
VF_AllocHemicube(Hemicube *hemicube)
{
  int     i,j,res,res2,size_top,size_side;
  float   *ptr;
  double  x,y,z,xsq,ysq,zsq;
  double  vf,vf_sum=0.0,sum=0.0;
  double  invrsq,areaPI,zareaPI,ds,ds2;
  Vector  *dptr;
  VFtopology* topology=VF_CurrentTopology();

  /*=========================================*/
  /* CALCULATE SOME PARAMETERS FOR LATER USE */
  /*=========================================*/
  res       = hemicube->resolution;
  res2      = hemicube->resolution/2;
  ds        = 2.0/res;
  ds2       = 1.0/res;
  areaPI    = ds*ds/M_PI;
  if (topology->geom==VF_2Dplanar) {
    size_top  = res;
    size_side = res2;
  } else {
    size_top  = res*res;
    size_side = res*res2;
  }
  /*================================*/
  /* ALLOCATE SPACE FOR THE ZBUFFER */
  /*================================*/
  hemicube->zbuffer = VF_Newd(size_top);
  /*====================================*/
  /* ALLOCATE SPACE FOR THE ITEM BUFFER */
  /*====================================*/
  hemicube->ibuffer = VF_Newi(size_top);
  /*====================================*/
  /*==========================================*/
  /* ALLOCATE SPACE FOR THE DELTA VIEWFACTORS */
  /*==========================================*/
  hemicube->top_deltaVF  = VF_Newf(size_top);
  hemicube->side_deltaVF = VF_Newf(size_side);
  /*==========================================*/
  /* ALLOCATE SPACE FOR THE DELTA VIEWFACTORS */
  /*==========================================*/
  hemicube->top_dir  = (Vector*)VF_Newv(size_top*sizeof(Vector));
  hemicube->side_dir = (Vector*)VF_Newv(size_side*sizeof(Vector));
  /*==============================================*/
  /* CALCULATE THE DELTA VIEWFACTORS FOR TOP FACE */
  /*==============================================*/
  dptr = hemicube->top_dir;
  ptr  = hemicube->top_deltaVF;
  if (topology->geom==VF_2Dplanar) {
    x = -1.0+ds2;
    for (j=0; j<res; j++) {
      *ptr = 0.0;
      xsq  = x*x;
      y    = -1.0+ds2;
      for (i=0; i<res; i++) {
        ysq     = y*y;
        invrsq  = 1.0/(xsq+ysq+1.0);
        vf      = areaPI*invrsq*invrsq;
        *ptr   += (float)vf;
        vf_sum += vf;
        y      += ds;
      }   
      invrsq  = 1.0/(xsq+1.0); 
      dptr->x = x*invrsq;
      dptr->y = 0.0;
      dptr->z = -invrsq;
      x      += ds;
      dptr++;
      ptr++;
    }
  } else {
    dptr = hemicube->top_dir;
    ptr  = hemicube->top_deltaVF;
    y    = -1.0+ds2;
    for (i=0; i<res; i++) {
      ysq = y*y;
      x   = -1.0+ds2;
      for (j=0; j<res; j++) {
        xsq      = x*x;
        invrsq   = 1.0/(xsq+ysq+1.0);
        vf       = areaPI*invrsq*invrsq;
        *ptr++   = vf;
        vf_sum  += vf;
        dptr->x  = x*invrsq;
        dptr->y  = y*invrsq;
        dptr->z  = -invrsq;
        x       += ds;
        dptr++;
      }
      y += ds;
    }
  }
  /*===============================================*/
  /* CALCULATE THE DELTA VIEWFACTORS FOR SIDE FACE */
  /*===============================================*/
  dptr = hemicube->side_dir;
  ptr  = hemicube->side_deltaVF;
  if (topology->geom==VF_2Dplanar) {
    z = ds2;
    for (i=0; i<res2; i++) {
      *ptr    = 0.0;
      zsq     = z*z; 
      zareaPI = z*areaPI;
      y       = -1.0+ds2;
      for (j=0; j<res; j++) {
        ysq     = y*y;
        invrsq  = 1.0/(1.0+ysq+zsq);
        vf      = zareaPI*invrsq*invrsq;
        *ptr   += (float)vf;
        sum    += vf;
        y      += ds;
      }   
      invrsq  = 1.0/(1.0+zsq);   
      dptr->x = 0.0;
      dptr->y = z*invrsq;
      dptr->z = -invrsq;
       z     += ds;
      dptr++;
      ptr++;
    }
    vf_sum += 2.0*sum;
  } else {
    z    = ds2;
    for (i=0; i<res2; i++) {
      zsq     = z*z; 
      zareaPI = z*areaPI;
      y       = -1.0+ds2;
      for (j=0; j<res; j++) {
        ysq     = y*y;
        invrsq  = 1.0/(1.0+ysq+zsq);
        vf      = zareaPI*invrsq*invrsq;
        *ptr++  = (float)vf;
        sum    += vf;
        dptr->x = y*invrsq;
        dptr->y = z*invrsq;
        dptr->z = -invrsq;
        y      += ds;
        dptr++;
      }
      z += ds;
    }
    vf_sum += 4.0*sum;
  }
  /*========================================================*/
  /* NORMALIZE THE DELTA-VIEWFACTORS TO THE HEMICUBE VF SUM */
  /*========================================================*/
  ptr = hemicube->top_deltaVF;
  for (i=0; i<size_top; i++) {
    *ptr++ /= vf_sum;
  }
  ptr = hemicube->side_deltaVF;
  for (i=0; i<size_side; i++) {
    *ptr++ /= vf_sum;
  }
  hemicube->nholes           = 0;
  hemicube->min_distance     = MAXDOUBLE;
  hemicube->min_equiv_radius = MAXDOUBLE;
}

void 
VF_AllocHemicubeAux(Hemicube *hemicube)
{
  int     i,j,res,res2,size_top,size_side;
  float   *ptr;
  double  x,y,z,xsq,ysq,zsq;
  double  vf,vf_sum=0.0,sum=0.0;
  double  invrsq,areaPI,zareaPI,ds,ds2;
  Vector  *dptr;
  VFtopology* topology=VF_CurrentTopology();

  /*=========================================*/
  /* CALCULATE SOME PARAMETERS FOR LATER USE */
  /*=========================================*/
  res       = hemicube->resolution;
  res2      = hemicube->resolution/2;
  ds        = 2.0/res;
  ds2       = 1.0/res;
  areaPI    = ds*ds/M_PI;
  if (topology->geom==VF_2Dplanar) {
    size_top  = res;
    size_side = res2;
  } else {
    size_top  = res*res;
    size_side = res*res2;
  }
  /*================================*/
  /* ALLOCATE SPACE FOR THE ZBUFFER */
  /*================================*/
  hemicube->zbuffer = VF_Newd(size_top+4*size_side);
  /*====================================*/
  /* ALLOCATE SPACE FOR THE ITEM BUFFER */
  /*====================================*/
  hemicube->ibuffer = VF_Newi(size_top+4*size_side);
  /*====================================*/
  /*==========================================*/
  /* ALLOCATE SPACE FOR THE DELTA VIEWFACTORS */
  /*==========================================*/
  hemicube->top_deltaVF  = VF_Newf(size_top);
  hemicube->side_deltaVF = VF_Newf(size_side);
  /*==========================================*/
  /* ALLOCATE SPACE FOR THE DELTA VIEWFACTORS */
  /*==========================================*/
  hemicube->top_dir  = (Vector*)VF_Newv(size_top*sizeof(Vector));
  hemicube->side_dir = (Vector*)VF_Newv(size_side*sizeof(Vector));
  /*==============================================*/
  /* CALCULATE THE DELTA VIEWFACTORS FOR TOP FACE */
  /*==============================================*/
  dptr = hemicube->top_dir;
  ptr  = hemicube->top_deltaVF;
  y    = -1.0+ds2;
  for (i=0; i<res; i++) {
    ysq = y*y;
    x   = -1.0+ds2;
    for (j=0; j<res; j++) {
      xsq     = x*x;
      invrsq  = 1.0/(xsq+ysq+1.0);
      vf      = areaPI*invrsq*invrsq;
      *ptr++  = (float)vf;
      vf_sum += vf;
            
      dptr->x = x*invrsq;
      dptr->y = y*invrsq;
      dptr->z = -invrsq;
      dptr++;
            
      x += ds;
    }
    y += ds;
  }
  /*===============================================*/
  /* CALCULATE THE DELTA VIEWFACTORS FOR SIDE FACE */
  /*===============================================*/
  dptr = hemicube->side_dir;
  ptr  = hemicube->side_deltaVF;
  z    = ds2;
  for (i=0; i<res2; i++) {
    zsq     = z*z; 
    zareaPI = z*areaPI;
    y       = -1.0+ds2;
    for (j=0; j<res; j++) {
      ysq     = y*y;
      invrsq  = 1.0/(1.0+ysq+zsq);
      vf      = zareaPI*invrsq*invrsq;
      *ptr++  = (float)vf;
      sum    += vf;
            
      dptr->x = y*invrsq;
      dptr->y = z*invrsq;
      dptr->z = -invrsq;
      dptr++;
            
      y += ds;
    }
    z += ds;
  }
  vf_sum += 4.0*sum;
  /*========================================================*/
  /* NORMALIZE THE DELTA-VIEWFACTORS TO THE HEMICUBE VF SUM */
  /*========================================================*/
  ptr = hemicube->top_deltaVF;
  for (i=0; i<size_top; i++) {
    *ptr++ /= vf_sum;
  }
  ptr = hemicube->side_deltaVF;
  for (i=0; i<size_side; i++) {
    *ptr++ /= vf_sum;
  }
  hemicube->nholes           = 0;
  hemicube->min_distance     = MAXDOUBLE;
  hemicube->min_equiv_radius = MAXDOUBLE;
}

void 
VF_FreeHemicube(Hemicube *hemicube)
{
  free(hemicube->zbuffer);
  free(hemicube->ibuffer);
  free(hemicube->top_deltaVF);
  free(hemicube->side_deltaVF);
  free(hemicube->top_dir);
  free(hemicube->side_dir);
}

void 
VF_JitterHemicube(Vector *n, Vector *u, Vector *v, Poly *poly)
{
  int     i, nn, nv, loop_cnt=0, max_loop=50;
  double  tmp;
  Vector  random_vector;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  VFtopology *topology=VF_CurrentTopology();

  if (VF_HemicubeJitter && topology->geom!=VF_2Dplanar) {
    do {
      random_vector.x = ran2(&enclosure->seed);
      random_vector.y = ran2(&enclosure->seed);
      random_vector.z = ran2(&enclosure->seed);
      V3_Cross(n,&random_vector,u);
      V3_Normalize(u,tmp);
      if (u->x==0.0 && u->y==0.0 && u->z==0.0) {
        loop_cnt++;
        if (loop_cnt>max_loop) {
          fprintf(stderr,"Zero length U vector detected too many times\n");
          fprintf(stderr,"    Proc = %d\n",VFLIB_Rank);
          fprintf(stderr,"    Poly Facet    = %d\n",poly->facet);
          fprintf(stderr,"    Poly SubFacet = %d\n",poly->sub_facet);
          fprintf(stderr,"    Poly Npoints  = %d\n",poly->np);
          for (nn=0; nn<poly->np; nn++) {
            fprintf(stderr,"    V[%d] = (%g, %g, %g)\n",
                    nn,poly->p[nn].x,poly->p[nn].y,poly->p[nn].z);
          }
          fprintf(stderr,"    Poly Normal      = (%g, %g, %g)\n",
                  n->x, n->y, n->z);
          fprintf(stderr,"    Last Rand Vector = (%g, %g, %g)\n",
                  random_vector.x,random_vector.y,random_vector.z);
          fprintf(stderr,"    Last U Vector    = (%g, %g, %g)\n",
                  u->x, u->y, u->z);
          nv = (topology->geom==VF_2Daxisym)?2:topology->facets[poly->facet].num_vertices;
          fprintf(stderr,"Facet %d:\n",poly->facet);
          fprintf(stderr,"   index              = %d\n",topology->facets[poly->facet].index);
          fprintf(stderr,"   patch_proc         = %d\n",topology->facets[poly->facet].patch_proc);
          fprintf(stderr,"   patch_gid          = %d\n",topology->facets[poly->facet].patch_gid);
          fprintf(stderr,"   patch_local_index  = %d\n",topology->facets[poly->facet].patch_local_index);
          fprintf(stderr,"   patch_global_index = %d\n",topology->facets[poly->facet].patch_global_index);
          fprintf(stderr,"   num_vertices       = %d (%d)\n",topology->facets[poly->facet].num_vertices,nv);
          for (i=0; i<nv; i++){
            nn = topology->facets[poly->facet].vertex_list[i];
            fprintf(stderr,"      %d)  %d -> (%f, %f, %f)\n",i,nn,
                    topology->x[nn],topology->y[nn],topology->z[nn]);
          }
          sleep(5);
          VF_Exit(0);
        }
      }
    } while (V3_Length(u)==0.0);
    V3_Cross(n,u,v);
    V3_Normalize(v,tmp);
    if (v->x==0.0 && v->y==0.0 && v->z==0.0) {
      fprintf(stderr,"Zero length V vector detected\n");
      VF_Exit(0);
    }
  } else {
    V3_Sub(&(poly->p[1]),&(poly->p[0]),u);
    V3_Normalize(u,tmp);
    V3_Cross(n,u,v);
    V3_Normalize(v,tmp);
  }
}

void VF_JitterOn(void)
{
  VF_HemicubeJitter = 1;
}

void VF_JitterOff(void)
{
  VF_HemicubeJitter = 0;
}
