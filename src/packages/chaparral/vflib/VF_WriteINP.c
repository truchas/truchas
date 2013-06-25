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
@(#)    $RCSfile: VF_WriteINP.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_WriteINP.c,v $
@(#)
@(#)    DESCRIPTION:  Write the current enclosure (as defined by MeshData)
@(#)    to a .inp (facet like) formatted file.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "vf.h"

void VF_WriteINP(char *filename)
{
  FILE        *fp;
  int         nnodes, i, j, k, nfacets, *nlist;
  Facet       *facet;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (VFLIB_Rank==0) {
    printf("\n   VF_WriteINP(%s)\n\n",filename);
    nfacets = topology->nfacets_g;
    if ((fp=fopen(filename,"w")) == NULL) {
      fprintf(stderr,"cannot open .inp file %s",filename);
      VF_Exit(-1);
    }
    nlist = VF_Newi(topology->nnodes);
    for (i=0; i<topology->nnodes; i++) {
      nlist[i] = -1;
     }
    for (facet=topology->facets, i=0; i<nfacets; i++, facet++) {
      for (j=0; j<facet->num_vertices; j++) {
        k        = facet->vertex_list[j] - topology->vertex_offset;
        nlist[k] = 0;
      }
    }
    for (nnodes=0, i=0; i<topology->nnodes; i++) {
      if (nlist[i]==0) nlist[i] = nnodes++ + topology->vertex_offset;
    }
    fprintf(fp,"test problem\n");
    fprintf(fp," %d %d 1 %d %d 0 %d 1 %d 1 %d\n",
            enclosure->vf_method, topology->geom, nnodes, nfacets,
            topology->nrotations,enclosure->hemicube.sub_divide,topology->debug_level);
    for (j=1, i=0; i<topology->nnodes; i++) {
      if (nlist[i]>=0) {
        if (topology->geom<3) {
          fprintf(fp," %6d\t%f\t%f\n",j,topology->x[i],topology->y[i]);
        } else {
          fprintf(fp," %6d\t%f\t%f\t%f\n",
                  j,topology->x[i],topology->y[i],topology->z[i]);
        }
        j++;
      }
    }
    for (facet=topology->facets, i=0; i<nfacets; i++, facet++) {
      fprintf(fp," %6d",i+1);
      for (j=0; j<facet->num_vertices; j++) {
        k = facet->vertex_list[j] - topology->vertex_offset;
        fprintf(fp,"\t%6d",nlist[k]);
      }
      fprintf(fp,"\n");
    }
    VF_Free(nlist);
    fclose(fp);
  } 
}
