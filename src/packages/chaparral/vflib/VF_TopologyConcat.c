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
@(#)    $RCSfile: VF_TopologyConcat.c,v $
@(#)    $Revision: 1.4 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_TopologyConcat.c,v $
@(#)
@(#)    DESCRIPTION:  Concatenate all the distributed surfaces 
@(#)    that reside separately on each processor onto each processor.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>

#include "vf.h"

/*==================================================================*/
/* NEED TO CONCATENATE ALL THE SURFACES ONTO EACH OF THE PROCESSORS */
/*==================================================================*/
void VF_TopologyConcat(int nsurf_g, int nsurf_l,int vertex_offset,
                       double x[], double y[], double z[], int c[],
                       int *f2p_map,
                       double **xx, double **yy, double **zz, int **cc, 
                       int **pmap, int* nnodes)
{
  int    i, j, n;
  int    num_nodes, node_offset, surf_offset, map_offset;
  int    num_unique_nodes, num_surf_nodes, node_num;
  int    *node_list=NULL, *conn_g=NULL, *map; 
  static int *Csize=NULL, *Coffset=NULL;
  static int *Msize=NULL, *Moffset=NULL;
  static int *Nsize=NULL, *Noffset=NULL;
  double *surf_x=NULL, *surf_y=NULL, *surf_z=NULL;
  
  VF_GetINTbuffer_ptr(&map);
  if (Csize==NULL)   Csize   = VF_Newi(VFLIB_Size);
  if (Coffset==NULL) Coffset = VF_Newi(VFLIB_Size);
  if (Msize==NULL)   Msize   = VF_Newi(VFLIB_Size);
  if (Moffset==NULL) Moffset = VF_Newi(VFLIB_Size);
  if (Nsize==NULL)   Nsize   = VF_Newi(VFLIB_Size);
  if (Noffset==NULL) Noffset = VF_Newi(VFLIB_Size);
  VF_AllgatherInt(&nsurf_l,Nsize,1);
  for (n=0, i=0; i<VFLIB_Size; i++) {
    n += Nsize[i];
  }
  if (n!=nsurf_g) {
    fprintf(stderr,"nsurf_g != SUM(nsurf_l)\n");
    VF_Exit(1);
  }
  VF_AllgatherInt(&nsurf_l,Csize,1);
  for (Coffset[0]=0, i=1; i<VFLIB_Size; i++) {
    Coffset[i] = Coffset[i-1]+Csize[i-1];
  }
  surf_offset = Coffset[VFLIB_Rank];
  VF_AllgatherInt(&nsurf_l,Msize,1);
  for (Moffset[0]=0, i=1; i<VFLIB_Size; i++) {
    Moffset[i] = Moffset[i-1]+Msize[i-1];
  }
  map_offset = Moffset[VFLIB_Rank];
  /*=================================================================*/
  /* COLLAPSE THE LOCAL CONNECTIVITY ARRAY TO ONLY ITS VALID ENTRIES */
  /*=================================================================*/
  node_list = VF_Newi(4*(nsurf_l+5));
  for (num_nodes=0, i=0; i<4*nsurf_l; i++) {
    if (c[i]==VF_INVALID) continue;
    node_list[num_nodes++] = c[i]-vertex_offset;
  }
  /*=======================================*/
  /* FIND THE NUMBER OF UNIQUE LOCAL NODES */
  /*=======================================*/
  if (num_nodes>0) {
    VF_SortIntegerArray(node_list,num_nodes);
    for (num_unique_nodes=1, i=1; i<num_nodes; i++) {
      if (node_list[i]!=node_list[num_unique_nodes-1]) {
        node_list[num_unique_nodes] = node_list[i];
        num_unique_nodes++;
      }
    }
  } else {
    num_unique_nodes = 0;
  }
  VF_AllgatherInt(&num_unique_nodes,Nsize,1);
  for (num_surf_nodes=0, i=0; i<VFLIB_Size; i++) {
    num_surf_nodes += Nsize[i];
  }
  for (Noffset[0]=0, i=1; i<VFLIB_Size; i++) {
    Noffset[i] = Noffset[i-1]+Nsize[i-1];
  }
  node_offset = Noffset[VFLIB_Rank];
  *nnodes     = num_surf_nodes;
  /*========================================================*/
  /* ALLOCATE GLOBAL COORDINATE ARRAYS AND OFFSET THE LOCAL */
  /* COORDINATE ARRAYS INTO THE GLOBAL COORDINATE ARRAYS    */
  /*========================================================*/
  surf_x = VF_Newd(num_surf_nodes+5);
  surf_y = VF_Newd(num_surf_nodes+5);
  surf_z = VF_Newd(num_surf_nodes+5);
  for (i=0; i<num_unique_nodes; i++) {
    surf_x[node_offset+i] = x[node_list[i]];
    surf_y[node_offset+i] = y[node_list[i]];
    surf_z[node_offset+i] = z[node_list[i]];
  }
  /*=========================================================*/
  /* ALLOCATE GLOBAL CONNECTIVITY ARRAY AND OFFSET THE LOCAL */
  /* CONNECTIVITY ARRAY INTO THE GLOBAL CONNECTIVITY ARRAY   */
  /* AND MAP THE LOCAL NODE NUMBERS TO GLOBAL NODE NUMBERS   */
  /*=========================================================*/
  conn_g = VF_Newi(4*(nsurf_g+5));
  for (i=0; i<nsurf_l; i++) {
    for (n=0; n<4; n++) {
      node_num = c[4*i+n];
      if (node_num!=VF_INVALID) {
        j = VF_LocateIndexInArray(node_num-vertex_offset,node_list,num_unique_nodes);
        if (j==VF_INVALID) {
          fprintf(stdout,"Can't find cross reference for node %d\n",node_num);
        } else {
          conn_g[4*(surf_offset+i)+n] = node_offset+j;
        }
      } else {
        conn_g[4*(surf_offset+i)+n] = -1;
      }
    }
  }
  VF_Free(node_list);
  /*========================================*/
  /* FILLIN THE COORDINATE AND CONNECTIVITY */
  /* ARRAYS FROM THE OTHER PROCESSORS       */
  /*========================================*/
  n = 4*nsurf_l;
  VF_AllgatherInt(&n,Csize,1);
  for (Coffset[0]=0, i=1; i<VFLIB_Size; i++) {
    Coffset[i] = Coffset[i-1]+Csize[i-1];
  }
  /*=================================================================*/
  /* HAVE TO DO THIS BECAUSE MPI_Allgatherv() SEEMS TO HAVE PROBLEMS */
  /* IF ANY OF THE SIZES ARE ZERO, PROBABLY BECAUSE THIS CAUSES MORE */
  /* THAN ONE PROCESSOR TO REFERENCE THE SAME MEMORY LOCATION - AT   */
  /* LEAST THIS IS THE CASE ON THE SP2 AND MPICH/Solaris             */
  /*=================================================================*/
  for (j=0, i=0; i<Msize[VFLIB_Rank]; j++, i++) {
    map[map_offset+j] = f2p_map[i];
  }
  for (j=0, i=0; i<VFLIB_Size; i++) {
    if (Csize[i]>0) j++;
  }
  if (j==VFLIB_Size) {
    VF_AllgathervInt(MPI_IN_PLACE,conn_g,
                     Csize[VFLIB_Rank],Csize,Coffset);
    VF_AllgathervDouble(MPI_IN_PLACE,surf_x,
                        Nsize[VFLIB_Rank],Nsize,Noffset);
    VF_AllgathervDouble(MPI_IN_PLACE,surf_y,
                        Nsize[VFLIB_Rank],Nsize,Noffset);
    VF_AllgathervDouble(MPI_IN_PLACE,surf_z,
                        Nsize[VFLIB_Rank],Nsize,Noffset);
    VF_AllgathervInt(MPI_IN_PLACE,map,
                     Msize[VFLIB_Rank],Msize,Moffset);
  } else {
    for (i=0; i<VFLIB_Size; i++) {
      if (Csize[i]>0) {
        VF_BroadcastInt(&conn_g[Coffset[i]],Csize[i],i);
        VF_BroadcastDouble(&surf_x[Noffset[i]],Nsize[i],i);
        VF_BroadcastDouble(&surf_y[Noffset[i]],Nsize[i],i);
        VF_BroadcastDouble(&surf_z[Noffset[i]],Nsize[i],i);
        VF_BroadcastInt(&map[Moffset[i]],Msize[i],i);
      }
      VF_Sync();
    }
  }
  *xx   = surf_x;
  *yy   = surf_y;
  *zz   = surf_z;
  *cc   = conn_g;
  *pmap = map;
}

