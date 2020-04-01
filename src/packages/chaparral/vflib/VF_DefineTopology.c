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
@(#)    $RCSfile: VF_DefineTopology.c,v $
@(#)    $Revision: 1.5 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_DefineTopology.c,v $
@(#)
@(#)    DESCRIPTION:  Setup the topology parameters to define the enclosure
@(#)    and control the viewfactor computation method.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vf.h"

void VF_DefineTopology(int enclosure, int geom_type, 
                       int nfacets, int nnodes, 
                       double *x, double *y, double *z, int *c, 
                       int vertex_offset, int *f2p_map, int nrotations, 
                       int x_mirror, int y_mirror, int z_mirror,
                       int bsp_depth, int bsp_length,
                       double spatial_tol, 
                       int debug_level)
{
  static int FacetListSize=0;
  int    ns, blocks, sector;
  int    i, j, n, nv, *vlist, facet_num, more, xm, ym, zm, r;
  int    vertex, vertex0, vertex1, num_vertices, nfacets_l, nfacets_g;
  int    *conn=NULL, *patch_map, *gid_map;
  float  *buffer;
  double clock0, clock1;
  double xmin, ymin, zmin, xmax, ymax, zmax, radius, value;
  double *xcoord=NULL, *ycoord=NULL, *zcoord=NULL;
  Facet  *facet, *base_facet;
  Poly   poly;
  Point  center;
  VFrow       *row;
  VFenclosure *encl=VF_GetEnclosure(enclosure);
  VFtopology  *topology=VF_GetTopology();

  /* NB: This modification alters the behavior of the function. It now
   * assumes a new topology is being defined with each call, instead of
   * adding to a topology defined by previous calls. */
  FacetListSize = 0;

  clock0 = VF_Clock();
  VF_GetSPbuffer0_ptr(&buffer);
  gid_map        = (int*)buffer;
  encl->topology = topology;
  nfacets_l      = nfacets;
  nfacets_g      = nfacets;
  VF_GlobalSumInt(&nfacets_g);
  
  if (debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncStart();
    printf("Processor %d has %d facets\n\n",VFLIB_Rank,nfacets_l);
    VF_PrintSyncEnd();
    
    VF_PrintSyncStart();
    printf("Processor %d --------------------------------\n",VFLIB_Rank);
    for (n=0; n<nfacets; n++) {
      j = VF_LocateIndexInSortedArray(f2p_map[n], gid_map, encl->npatches_g);
      printf("facet %d:\n",n);
      printf("   patch_gid          = %d\n",f2p_map[n]);
      printf("   patch_index        = %d\n",j);
      for (i=0; i<4; i++){
        if (c[4*n+i] == -1) break;
        vertex = c[4*n+i];
        printf("      %d)  (%f, %f, %f)\n",i,
               x[vertex],y[vertex],z[vertex]);
      }
    }
    VF_PrintSyncEnd();
  }
  
  if (VFLIB_Size>1) {
    VF_TopologyConcat(nfacets_g, nfacets_l, vertex_offset, x, y, z, c, f2p_map, 
                      &xcoord, &ycoord, &zcoord, &conn, &patch_map, &nnodes);
  } else {
    xcoord = VF_Newd(nnodes);
    ycoord = VF_Newd(nnodes);
    zcoord = VF_Newd(nnodes);
    memcpy(xcoord,x,nnodes*sizeof(double));
    memcpy(ycoord,y,nnodes*sizeof(double));
    memcpy(zcoord,z,nnodes*sizeof(double));
    conn = VF_Newi(4*nfacets+1); /* not sure of why the "+1" */
    memcpy(conn,c,4*nfacets*sizeof(int));
    patch_map = f2p_map;
    if (vertex_offset) {
      xcoord--;
      ycoord--;
      zcoord--;
    }
  }
  if (VFLIB_Size>1) {
    blocks    = nfacets_g/2;
    nfacets_l = 2*(int)(blocks/VFLIB_Size);
    more      = blocks%VFLIB_Size;
    if (VFLIB_Rank<more) nfacets_l += 2;
  } else {
    nfacets_l = nfacets_g;
  }
  ns = 1;
  if (x_mirror!=0) ns *= 2;
  if (y_mirror!=0) ns *= 2;
  if (z_mirror!=0) ns *= 2;
  switch (geom_type) {
  case VF_2Daxisym:
    if (nrotations%2 != 0) nrotations++;
    break;
  case VF_2Dplanar:
    nrotations = 1;
    break;
  }
  topology->nsections     = nrotations*ns;
  topology->xmirror       = x_mirror+1;
  topology->ymirror       = y_mirror+1;
  topology->zmirror       = z_mirror+1;
  topology->nrotations    = nrotations;
  topology->theta         = 2.0*M_PI/(double)(nrotations);
  topology->debug_level   = debug_level;
  topology->nonblocking   = encl->nonblocking;
  topology->geom          = geom_type;
  topology->nfacets_g     = nfacets_g*topology->nsections;
  topology->nfacets_l     = nfacets_l*topology->nsections;
  topology->nfacets_base  = nfacets_g;
  topology->nnodes        = nnodes;
  topology->x             = xcoord;
  topology->y             = ycoord;
  topology->z             = zcoord;
  topology->connect       = conn;
  topology->vertex_offset = vertex_offset;
  topology->bsp_depth     = bsp_depth;
  topology->bsp_length    = bsp_length;
  topology->spatial_tol   = spatial_tol;
  if (FacetListSize==0) {
    FacetListSize = topology->nfacets_g*(int)sizeof(Facet);
    topology->facets = (Facet *)VF_Newv(FacetListSize);
  } else if (topology->nfacets_g*(int)sizeof(Facet)>FacetListSize) {
    FacetListSize = topology->nfacets_g*(int)sizeof(Facet);
    topology->facets = (Facet *)VF_ReNewv(topology->facets,FacetListSize);
  }
  sector    = 0;
  facet_num = 0;
  facet     = topology->facets;
  vlist     = topology->connect;
  for (n=0; n<topology->nfacets_base; n++) {
    switch (topology->geom) {
    case VF_2Dplanar:
      num_vertices = 2;
      break;
    case VF_2Daxisym:
      vertex0 = vlist[0];
      vertex1 = vlist[1];
      if (fabs(topology->x[vertex0])<topology->spatial_tol &&
          fabs(topology->x[vertex1])<topology->spatial_tol) {
        printf("topologyInit() - Bad element, both points at R=0\n");
        VF_Exit(1);
      } else if (fabs(topology->x[vertex0])<topology->spatial_tol) {
        num_vertices = 3;
      } else if (fabs(topology->x[vertex1])<topology->spatial_tol) {
        num_vertices = 3;
      } else {
        num_vertices = 4;
      }
      break;
    case VF_3D:
      if (vlist[3]==-1) {
        num_vertices = 3;
      } else {
        num_vertices = 4;
      }
    }
    row                       = VF_LocateLocalRow(patch_map[n]);
    facet->sector             = sector;
    facet->num_vertices       = num_vertices;
    facet->vertex_list        = vlist;
    facet->mask               = 0;
    facet->index              = facet_num++;
    facet->patch_gid          = patch_map[n];
    facet->patch_proc         = row!=NULL?VFLIB_Rank:-1;
    facet->patch_local_index  = row!=NULL?row->local_index:-1;
    facet->patch_global_index = VF_LocateIndexInSortedArray(patch_map[n], 
                                                            gid_map, 
                                                            encl->npatches_g);
    /*
    VF_GlobalMaxInt(&(facet->patch_proc));
    */
    if (num_vertices!=4 || topology->geom==VF_2Daxisym) {
      facet->mask |= VF_FACET_MASK_PLANAR;
    } else {
      if (VF_IsFacetPlanar(facet)) facet->mask |= VF_FACET_MASK_PLANAR;
    }
    VF_FacetToPoly(facet, &poly);
    if (facet->mask & VF_FACET_MASK_PLANAR) {
      VF_PolyNormal(&poly, &(facet->normal));
      facet->d = -V3_Dot(&(facet->normal), &(poly.p[0]));
    } else {
      VF_GenericPolyNormal(&poly, &(facet->normal));
      VF_PolyCenter(&poly, &center);
      facet->d = -V3_Dot(&(facet->normal), &center);
      value    = V3_Dot(&(poly.p[0]), &facet->normal) + facet->d;
      if (value < 0.0) {
        facet->d = -V3_Dot(&(facet->normal), &(poly.p[0]));
      }
    }
    VF_ComputeFacetBoundingBox(&(facet->bounds), facet);
    vlist += 4;
    facet++;
  }

  for (xm=0; xm<topology->xmirror; xm++) {
    for (ym=0; ym<topology->ymirror; ym++) {
      for (zm=0; zm<topology->zmirror; zm++) {
        for (r=0; r<topology->nrotations; r++) {
          sector = r+(xm<<16)+(ym<<20)+(zm<<24);
          if (sector==0) continue;
          base_facet = topology->facets;
          for (n=0; n<nfacets_g; n++) {
            facet->sector             = sector;
            facet->num_vertices       = base_facet->num_vertices;
            facet->vertex_list        = base_facet->vertex_list;
            facet->mask               = base_facet->mask;
            facet->index              = facet_num++;
            facet->patch_proc         = base_facet->patch_proc;
            facet->patch_gid          = base_facet->patch_gid;
            facet->patch_local_index  = base_facet->patch_local_index;
            facet->patch_global_index = base_facet->patch_global_index;
            VF_FacetToPoly(facet, &poly);
            if (facet->mask & VF_FACET_MASK_PLANAR) {
              VF_PolyNormal(&poly, &(facet->normal));
              facet->d = -V3_Dot(&(facet->normal), &(poly.p[0]));
            } else {
              VF_GenericPolyNormal(&poly, &(facet->normal));
              VF_PolyCenter(&poly, &center);
              facet->d = -V3_Dot(&(facet->normal), &center);
              value    = V3_Dot(&(poly.p[0]), &facet->normal) + facet->d;
              if (value < 0.0) {
                facet->d = -V3_Dot(&(facet->normal), &(poly.p[0]));
              }
            }
            VF_ComputeFacetBoundingBox(&(facet->bounds), facet);
            base_facet++;
            facet++;
          }
        }
      }
    }
  }
  VF_SortFacets();
  for (facet=topology->facets, n=0; n<topology->nfacets_g; n++, facet++) {
    facet->index = n;
  }
  topology->nfacets_l = 0;
  for (facet=topology->facets, n=0; n<topology->nfacets_base; n++, facet++) {
    if (VFLIB_Rank==facet->patch_proc) {
      topology->nfacets_l++;
    }
  }

  radius = MAXDOUBLE;
  for (facet=topology->facets, n=0; n<topology->nfacets_g; n++, facet++) {
    radius = MIN(radius,sqrt(VF_FacetArea(facet)/M_PI));
    if (topology->geom==VF_2Daxisym) {
      nv = 2;
    } else {
      nv = facet->num_vertices;
    }
    for (i=0; i<nv; i++){
      vertex = facet->vertex_list[i];
      if (n==0 && i==0) {
        xmin = topology->x[vertex];
        ymin = topology->y[vertex];
        xmax = topology->x[vertex];
        ymax = topology->y[vertex];
        switch (topology->geom) {
        case VF_2Daxisym:
          zmin = xmin;
          zmax = xmax;
          break;
        case VF_2Dplanar:
          zmin = 0.0;
          zmax = 0.0;
          break;
        case VF_3D:
          zmin = topology->z[vertex];
          zmax = topology->z[vertex];
          break;
        }
      } else {
        xmin = MIN(xmin,(double)(topology->x[vertex]));
        ymin = MIN(ymin,(double)(topology->y[vertex]));
        xmax = MAX(xmax,(double)(topology->x[vertex]));
        ymax = MAX(ymax,(double)(topology->y[vertex]));
        if (topology->geom==VF_3D) {
          zmin = MIN(zmin,(double)(topology->z[vertex]));
          zmax = MAX(zmax,(double)(topology->z[vertex]));
        } else if (topology->geom==VF_2Daxisym) {
          xmin = MIN(xmin,-(double)(topology->x[vertex]));
          ymin = MIN(ymin,-(double)(topology->y[vertex]));
          xmax = MAX(xmax,-(double)(topology->x[vertex]));
          ymax = MAX(ymax,-(double)(topology->y[vertex]));
          zmin = xmin;
          zmax = xmax;
        }
      }
    }
    topology->bounds.xmin = xmin;
    topology->bounds.ymin = ymin;
    topology->bounds.zmin = zmin;
    topology->bounds.xmax = xmax;
    topology->bounds.ymax = ymax;
    topology->bounds.zmax = zmax;
  }
  if (debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncStart();
    printf("Processor %d --------------------------------\n",VFLIB_Rank);
    for (n=0; n<topology->nfacets_g; n++) {
      facet = &(topology->facets[n]);
      nv    = (topology->geom==VF_2Daxisym)?2:facet->num_vertices;
      printf("facet %d:\n",n);
      printf("   index              = %d\n",facet->index);
      printf("   patch_proc         = %d\n",facet->patch_proc);
      printf("   patch_gid          = %d\n",facet->patch_gid);
      printf("   patch_local_index  = %d\n",facet->patch_local_index);
      printf("   patch_global_index = %d\n",facet->patch_global_index);
      printf("   mask               = %d\n",facet->mask);
      printf("   sector             = %d\n",facet->sector);
      printf("   d                  = %f\n",facet->d);
      printf("   normal             = (%f, %f, %f)\n",facet->normal.x,
                                                      facet->normal.y,
                                                      facet->normal.z);
      printf("   num_vertices       = %d (%d)\n",facet->num_vertices,nv);
      for (i=0; i<nv; i++){
        vertex = facet->vertex_list[i];
        /*printf("      %d)  %d -> (%f, %f, %f)\n",i,vertex,*/
        printf("      %d)  (%f, %f, %f)\n",i,
               topology->x[vertex],topology->y[vertex],topology->z[vertex]);
      }
    }
    VF_PrintSyncEnd();
  }
  VF_BSP_CreateTree(topology);
  clock1 = VF_Clock();
  encl->time_vf_init += clock1-clock0;
}
