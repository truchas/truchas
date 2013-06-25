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
@(#)    $RCSfile: VF_WriteGenesis.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_WriteGenesis.c,v $
@(#)
@(#)    DESCRIPTION:  Write the current enclosure (as defined by MeshData)
@(#)    to an ExodusII formatted file.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>

#ifndef VF_NO_EXODUS_OUTPUT

#include "vf.h"

#undef TRUE
#include "exodusII.h"

void VF_WriteGenesis(char *filename)
{
  char       exo_title[81], *cnames[3];
  int        i, j, k, m, n, dim, nblocks, exoid, exoerr, ws1, ws2;
  int        nbar, ntri, nquad, nnodes, nfacets, *elist, *slist;
  int        *conn, npe[3], nelems[3], elemID[] = {100, 200, 300};
  double     *xx=NULL, *yy=NULL, *zz=NULL;
  Poly       poly;
  Facet      *facet;
  VFtopology *topology=VF_CurrentTopology();

  if (VFLIB_Rank==0) {
    nbar    = 0;
    ntri    = 0;
    nquad   = 0;
    nblocks = 0;
    nfacets = topology->nfacets_g;
    nnodes  = 0;
    for (facet=topology->facets, i=0; i<nfacets; i++, facet++) {
      switch (facet->num_vertices) {
      case 2:
        dim = 2;
        nbar++;
        break;
      case 3:
        dim = 3;
        ntri++;
        break;
      case 4:
        dim = 3;
        nquad++;
        break;
      }
      nnodes += facet->num_vertices;
    }
    if (nbar>0) {
      npe[nblocks] = 2;
      nelems[nblocks] = nbar;
      nblocks++;
    }
    if (ntri>0) {
      npe[nblocks] = 3;
      nelems[nblocks] = ntri;
      nblocks++;
    }
    if (nquad>0) {
      npe[nblocks] = 4;
      nelems[nblocks] = nquad;
      nblocks++;
    }
    k  = 0;
    xx = VF_Newd(nnodes);
    yy = VF_Newd(nnodes);
    zz = VF_Newd(nnodes);
    for (facet=topology->facets, i=0; i<nfacets; i++, facet++) {
      VF_FacetToPoly(facet, &poly);
      for (j=0; j<poly.np; j++) {
        xx[k] = poly.p[j].x;
        yy[k] = poly.p[j].y;
        zz[k] = poly.p[j].z;
        k++;
      }
    }
    sprintf(exo_title,"enclosure topology");
    ws1       = sizeof(double);
    ws2       = sizeof(float);
    cnames[0] = "x-coord";
    cnames[1] = "y-coord";
    cnames[2] = "z-coord";
    exoid  = ex_create(filename, EX_CLOBBER, &ws1, &ws2);
    exoerr = ex_put_init(exoid,exo_title,dim,nnodes,
                         nfacets,nblocks,0,1);
    exoerr = ex_put_coord(exoid,xx,yy,zz);
    exoerr = ex_put_coord_names(exoid,cnames);
    VF_Free(xx);
    VF_Free(yy);
    VF_Free(zz);
    for (i=0; i<nblocks; i++) {
      n    = 1;
      m    = 0;
      conn = VF_Newi(npe[i]*nelems[i]);
      for (facet=topology->facets, j=0; j<nfacets; j++, facet++) {
        if (facet->num_vertices==npe[i]) {
          for (k=0; k<npe[i]; k++) {
            conn[m++] = n+k;
          }
        }
        n += facet->num_vertices;
      }
      switch (npe[i]) {
      case 2:
        exoerr = ex_put_elem_block(exoid,elemID[i],"BAR",
                                   nelems[i],npe[i],0);
        break;
      case 3:
        exoerr = ex_put_elem_block(exoid,elemID[i],"TRIANGLE",
                                   nelems[i],npe[i],0);
        break;
      case 4:
        exoerr = ex_put_elem_block(exoid,elemID[i],"SHELL",
                                   nelems[i],npe[i],0);
        break;
      }
      exoerr = ex_put_elem_conn(exoid,elemID[i],conn);
      VF_Free(conn);
    }
    elist = VF_Newi(nfacets);
    slist = VF_Newi(nfacets);
    for (i=0; i<nfacets; i++) {
      elist[i] = i+1;
      slist[i] = 1;
    }
    exoerr = ex_put_side_set_param(exoid,10,nfacets,0);
    exoerr = ex_put_side_set(exoid,10,elist,slist);
    VF_Free(elist);
    VF_Free(slist);
    ex_close(exoid);
    if (exoerr) exoerr = 0;
  }
  VF_Sync();
}

#else

void VF_WriteGenesis(char *filename) 
{
  fprintf(stderr,"WARNING: VF_WriteGenesis(%s) not supported\n",filename);
}

#endif
