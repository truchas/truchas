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
@(#)    $RCSfile: VF_TopologyUtils.c,v $
@(#)    $Revision: 1.5.4.1 $  $Date: 2006/04/10 22:51:03 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_TopologyUtils.c,v $
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

int VF_CompareFacets(Facet*, Facet*);

static VFtopology VF_DATA_Topology;

VFtopology*
VF_GetTopology()
{
  return &VF_DATA_Topology;
}

void
VF_ResetTopology(int enclosure)
{
  VFtopology  *topology=&VF_DATA_Topology;
  BinTree     *BSPtree=&(topology->BSP_tree);
  VFenclosure* encl = VF_GetEnclosure(enclosure);
  encl->topology = NULL;

  BSPtree = &(topology->BSP_tree);
  VF_BSP_DeleteTree(BSPtree->root);
  VF_Free(BSPtree->facetLinkList);
  VF_Free(topology->BSP_stack);
  /* We also want to free these when VF_LIB=1 but in this case
     the addresses in X,Y,Z are sometimes shifted for vertex offset */
  if (VFLIB_Size>1) {
    VF_Free(topology->x);
    VF_Free(topology->y);
    VF_Free(topology->z);
    VF_Free(topology->connect);
  }
}

void
VF_DeleteTopology()
{
  VFtopology  *topology=&VF_DATA_Topology;
  /*
  BinTree     *BSPtree=&(topology->BSP_tree);
  VF_BSP_DeleteTree(BSPtree->root);
  VF_Free(BSPtree->facetLinkList);
  VF_Free(topology->BSP_stack);
  */
  VF_Free(topology->facets);
}

void VF_SortFacets(void)
{
  int        l,j,ir,i;
  int        nelements;
  Facet      rra;
  VFtopology *topology=VF_CurrentTopology();

  nelements = topology->nfacets_g;
  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        rra = topology->facets[--l-1];
      } else {
        rra                    = topology->facets[ir-1];
        topology->facets[ir-1] = topology->facets[0];
        if (--ir == 1) {
          topology->facets[0] = rra;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (VF_CompareFacets(&(topology->facets[j-1]), &(topology->facets[j]))<0)) ++j;
        if (VF_CompareFacets(&(rra), &(topology->facets[j-1]))<0) {
          topology->facets[i-1] = topology->facets[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      topology->facets[i-1] = rra;
    }
  }
}

int VF_CompareFacets(Facet* f1, Facet* f2)
{
  if (f1->sector<f2->sector)  return -1;
  if (f1->sector>f2->sector)  return  1;

  if (f1->patch_global_index<f2->patch_global_index)  return -1;
  if (f1->patch_global_index>f2->patch_global_index)  return  1;

  if (f1->index<f2->index)  return -1;
  if (f1->index>f2->index)  return  1;

  return(0);
}



