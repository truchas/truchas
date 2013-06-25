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
@(#)    $RCSfile: VF_WriteExodus.c,v $
@(#)    $Revision: 1.5.2.1 $  $Date: 2006/04/10 22:51:03 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_WriteExodus.c,v $
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

#define VFLIB_NelemVar 4

#define VFLIB_NelemVarAux 10

void VF_WriteExodus(char *filename)
{
  char        exo_title[81], *cnames[3], *vnames[VFLIB_NelemVar];
  int         i, j, k, m, n, dim, nblocks, offset, exoid, exoerr=0, ws1, ws2;
  int         nverts, nbar, ntri, nquad, nfacets, *key=NULL, *index=NULL;
  int         *truth=NULL, *elist=NULL, *slist=NULL, *conn=NULL, *emap=NULL;
  int         nelems[3], npe[3], elemID[] = {100, 200, 300};
  double      *buffer=NULL, *exo_data=NULL;
  double      *xx=NULL, *yy=NULL, *zz=NULL, value, time=0.0;
  Facet       *facet=NULL;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  nfacets = topology->nfacets_base;
  key     = VF_Newi(nfacets);
  if (VFLIB_Rank==0) {
    nbar     = 0;
    ntri     = 0;
    nquad    = 0;
    nblocks  = 0;
    exo_data = VF_Newd(nfacets);
    for (facet=topology->facets, i=0; i<nfacets; i++, facet++) {
      if (topology->geom==VF_2Daxisym) {
        nverts = 2;
      } else {
        nverts = facet->num_vertices;
      }
      switch (nverts) {
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
    sprintf(exo_title,"viewfactor geometry output");
    ws1       = sizeof(double);
    ws2       = sizeof(float);
    time      = 0.0;
    xx        = topology->x;
    yy        = topology->y;
    zz        = topology->z;
    if (VFLIB_Size==1 && topology->vertex_offset) {
      xx++;
      yy++;
      zz++;
    }
    cnames[0] = "x-coord";
    cnames[1] = "y-coord";
    cnames[2] = "z-coord";
    exoid  = ex_create(filename, EX_CLOBBER, &ws1, &ws2);
    exoerr = ex_put_init(exoid,exo_title,dim,topology->nnodes,
                         nfacets,nblocks,0,1);
    exoerr = ex_put_coord(exoid,xx,yy,zz);
    exoerr = ex_put_coord_names(exoid,cnames);
    for (i=0; i<nblocks; i++) {
      n    = 0;
      m    = 0;
      conn = VF_Newi(npe[i]*nelems[i]);
      for (facet=topology->facets, j=0; j<nfacets; j++, facet++) {
        if (topology->geom==VF_2Daxisym) {
          nverts = 2;
        } else {
          nverts = facet->num_vertices;
        }
        if (nverts==npe[i]) {
          for (k=0; k<npe[i]; k++) {
            conn[m++] = facet->vertex_list[k]+1;
          }
          key[j] = n++;
        }
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
    vnames[0] = "Patch GID";
    vnames[1] = "CHAPARRAL Proc";
    vnames[2] = "Raw Rowsum Error";
    vnames[3] = "Smoothed Rowsum Error";
    exoerr = ex_put_var_param(exoid,"e",VFLIB_NelemVar);
    exoerr = ex_put_var_names(exoid,"e",VFLIB_NelemVar,vnames);
    truth  = VF_Newi(VFLIB_NelemVar*nblocks);
    for (i=0; i<VFLIB_NelemVar*nblocks; i++) truth[i] = 1;
    exoerr = ex_put_elem_var_tab(exoid,nblocks,VFLIB_NelemVar,truth);
    exoerr = ex_put_time(exoid, 1, &time);
    VF_Free(truth);
  }
  offset = 0;
  index  = VF_Newi(nfacets);
  buffer = VF_Newd(nfacets);
  emap   = VF_Newi(nfacets);
  VF_Sync();
  VF_BroadcastInt(&nblocks,1,0);
  VF_BroadcastInt(key,nfacets,0);
  VF_BroadcastInt(npe,nblocks,0);
  for (i=0; i<nblocks; i++) {
    for (j=0; j<VFLIB_NelemVar; j++) {
      for (n=0, k=0; k<nfacets; k++) {
        facet = &(topology->facets[k]);
        if (topology->geom==VF_2Daxisym) {
          nverts = 2;
        } else {
          nverts = facet->num_vertices;
        }
        if (nverts==npe[i] && facet->patch_proc==VFLIB_Rank) {
          switch (j) {
          case 0:
            value = (double)facet->patch_gid;
            break;
          case 1:
            value = (double)VFLIB_Rank;
            break;
          case 2:
            value = VF_GetRawRowsum(facet->patch_local_index)-1.0;
            break;
          case 3:
            value = VF_GetSmoothedRowsum(facet->patch_local_index)-1.0;
            break;
          }
          if (VFLIB_Rank==0) {
            exo_data[key[k]] = value;
          } else {
            buffer[n] = value;
            index[n]  = key[k];
            n++;
          }
        }
      }
      /*VF_GatherExodusData(exo_data, buffer, index, n, VFLIB_NelemVar);*/
      VF_GatherExodusData(exo_data, buffer, index, n, 1);
      if (VFLIB_Rank==0) {
        if (j==0) {
          for (n=0; n<nelems[i]; ++n) {
            emap[offset+n] = (int)exo_data[n];
          }
          offset += nelems[i];
        }
        exoerr = ex_put_elem_var(exoid, 1, j+1, elemID[i], 
                                 nelems[i], exo_data);
      }
    }
  }
  if (VFLIB_Rank==0) {
    exoerr = ex_put_elem_num_map(exoid, emap);
    ex_close(exoid);
    VF_Free(exo_data);
  }
  VF_Free(emap);
  VF_Free(index);
  VF_Free(buffer);
  VF_Free(key);
  if (exoerr) exoerr = 0;
}

int VF_WriteExodusAux1(char *filename, int step, double time)
{
  char        exo_title[81], *cnames[3], *vnames[VFLIB_NelemVarAux];
  int         i, j, k, m, n, dim, nblocks, offset, exoid, exoerr=0, ws1, ws2;
  int         nverts, nbar, ntri, nquad, nfacets, *key=NULL, *index=NULL;
  int         *truth=NULL, *elist=NULL, *slist=NULL, *conn=NULL;
  int         nelems[3], npe[3], elemID[] = {100, 200, 300};
  double      *buffer=NULL, *exo_data=NULL;
  double      *xx=NULL, *yy=NULL, *zz=NULL, value;
  Facet       *facet=NULL;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  nfacets = topology->nfacets_base;
  key     = VF_Newi(nfacets);
  if (VFLIB_Rank==0) {
    nbar     = 0;
    ntri     = 0;
    nquad    = 0;
    nblocks  = 0;
    exo_data = VF_Newd(nfacets);
    for (facet=topology->facets, i=0; i<nfacets; i++, facet++) {
      if (topology->geom==VF_2Daxisym) {
        nverts = 2;
      } else {
        nverts = facet->num_vertices;
      }
      switch (nverts) {
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
    sprintf(exo_title,"viewfactor geometry output");
    ws1  = sizeof(double);
    ws2  = sizeof(float);
    time = 0.0;
    xx   = topology->x;
    yy   = topology->y;
    zz   = topology->z;
    if (VFLIB_Size==1 && topology->vertex_offset) {
      xx++;
      yy++;
      zz++;
    }
    cnames[0] = "x-coord";
    cnames[1] = "y-coord";
    cnames[2] = "z-coord";
    exoid  = ex_create(filename, EX_CLOBBER, &ws1, &ws2);
    exoerr = ex_put_init(exoid,exo_title,dim,topology->nnodes,
                         nfacets,nblocks,0,1);
    exoerr = ex_put_coord(exoid,xx,yy,zz);
    exoerr = ex_put_coord_names(exoid,cnames);
    for (i=0; i<nblocks; i++) {
      n    = 0;
      m    = 0;
      conn = VF_Newi(npe[i]*nelems[i]);
      for (facet=topology->facets, j=0; j<nfacets; j++, facet++) {
        if (topology->geom==VF_2Daxisym) {
          nverts = 2;
        } else {
          nverts = facet->num_vertices;
        }
        if (nverts==npe[i]) {
          for (k=0; k<npe[i]; k++) {
            conn[m++] = facet->vertex_list[k]+1;
          }
          key[j] = n++;
        }
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
    vnames[0] = "Tmp Key";
    vnames[1] = "Tmp Blk";
    vnames[2] = "Tmp PID";
    vnames[3] = "Tmp GID";
    vnames[4] = "Patch GID";
    vnames[5] = "CHAPARRAL Proc";
    vnames[6] = "Raw Rowsum Error";
    vnames[7] = "Smoothed Rowsum Error";
    vnames[8] = "Radiosity";
    vnames[9] = "Flux";
    exoerr = ex_put_var_param(exoid,"e",VFLIB_NelemVarAux);
    exoerr = ex_put_var_names(exoid,"e",VFLIB_NelemVarAux,vnames);
    truth  = VF_Newi(VFLIB_NelemVarAux*nblocks);
    for (i=0; i<VFLIB_NelemVarAux*nblocks; i++) truth[i] = 1;
    exoerr = ex_put_elem_var_tab(exoid,nblocks,VFLIB_NelemVarAux,truth);
    exoerr = ex_put_time(exoid, step, &time);
    VF_Free(truth);
    
    for (k=0; k<nfacets; k++) {
      exo_data[k] = key[k];
    }
    for (offset=0, i=0; i<nblocks; i++) {
      exoerr = ex_put_elem_var(exoid, 1, 1, elemID[i], 
                               nelems[i], &exo_data[offset]);
      offset += nelems[i];
    }
    
    for (k=0; k<nfacets; k++) {
      facet = &(topology->facets[k]);
      if (topology->geom==VF_2Daxisym) {
        nverts = 2;
      } else {
        nverts = facet->num_vertices;
      }
      for (i=0; i<nblocks; i++) {
        if (nverts==npe[i]) {
          exo_data[k] = i;
          break;
        }
      }
    }
    for (offset=0, i=0; i<nblocks; i++) {
      exoerr = ex_put_elem_var(exoid, 1, 2, elemID[i], 
                               nelems[i], &exo_data[offset]);
      offset += nelems[i];
    }
    
    for (k=0; k<nfacets; k++) {
      facet = &(topology->facets[k]);
      exo_data[k] = facet->patch_proc;
    }
    for (offset=0, i=0; i<nblocks; i++) {
      exoerr = ex_put_elem_var(exoid, 1, 3, elemID[i], 
                               nelems[i], &exo_data[offset]);
      offset += nelems[i];
    }
      
    for (k=0; k<nfacets; k++) {
      facet = &(topology->facets[k]);
      exo_data[k] = facet->patch_global_index;
    }
    for (offset=0, i=0; i<nblocks; i++) {
      exoerr = ex_put_elem_var(exoid, 1, 4, elemID[i], 
                               nelems[i], &exo_data[offset]);
      offset += nelems[i];
    }
  }
  index  = VF_Newi(nfacets);
  buffer = VF_Newd(nfacets);
  VF_Sync();
  VF_BroadcastInt(&nblocks,1,0);
  VF_BroadcastInt(key,nfacets,0);
  VF_BroadcastInt(npe,nblocks,0);
  for (i=0; i<nblocks; i++) {
    for (j=4; j<VFLIB_NelemVarAux-2; j++) {
      for (n=0, k=0; k<nfacets; k++) {
        facet = &(topology->facets[k]);
        if (topology->geom==VF_2Daxisym) {
          nverts = 2;
        } else {
          nverts = facet->num_vertices;
        }
        if (nverts==npe[i] && facet->patch_proc==VFLIB_Rank) {
          switch (j) {
          case 4:
            value = (double)facet->patch_gid;
            break;
          case 5:
            value = (double)VFLIB_Rank;
            break;
          case 6:
            value = VF_GetRawRowsum(facet->patch_local_index)-1.0;
            break;
          case 7:
            value = VF_GetSmoothedRowsum(facet->patch_local_index)-1.0;
            break;
          }
          if (VFLIB_Rank==0) {
            exo_data[key[k]] = value;
          } else {
            buffer[n] = value;
            index[n]  = key[k];
            n++;
          }
        }
      }
      /*VF_GatherExodusData(exo_data, buffer, index, n, VFLIB_NelemVarAux-4-2);*/
      VF_GatherExodusData(exo_data, buffer, index, n, 1);
      if (VFLIB_Rank==0) {
        exoerr = ex_put_elem_var(exoid, step, j+1, elemID[i], 
                                 nelems[i], exo_data);
      }
    }
  }
  if (VFLIB_Rank==0) {
    ex_update(exoid);
    VF_Free(exo_data);
  }
  VF_Free(index);
  VF_Free(buffer);
  VF_Free(key);
  if (exoerr) exoerr = 0;
  return exoid;
}

void VF_WriteExodusAux2(int encl, int exoid, int step, 
                        double time, double radq[], double radj[])
{
  char        cdum[81], elem_type[81];
  int         max_surfaces;
  int         i, j, k, m, n, dim, offset, exoerr=0;
  int         numnodes, numelems, numblocks, numattr;
  int         *key=NULL, *proc_id=NULL, *g_index=NULL, *block=NULL, *index=NULL;
  int         nelems[3], npe[3], elemID[] = {100, 200, 300};
  float       fdum;
  double      *buffer=NULL, *exo_data=NULL;
  double      *radq_g=NULL, *radj_g=NULL, value;
  VFenclosure *enclosure=VF_GetEnclosure(encl);

  max_surfaces = VF_MaxPatches();
  radj_g       = VF_Newd(max_surfaces);
  radq_g       = VF_Newd(max_surfaces);
 
  VF_Gather_ExodusVals(radj, radq, radj_g, radq_g);

  if (VFLIB_Rank==0) {
    exoerr   = ex_inquire(exoid, EX_INQ_DIM,      &dim,       &fdum, cdum);
    exoerr   = ex_inquire(exoid, EX_INQ_NODES,    &numnodes,  &fdum, cdum);
    exoerr   = ex_inquire(exoid, EX_INQ_ELEM,     &numelems,  &fdum, cdum);
    exoerr   = ex_inquire(exoid, EX_INQ_ELEM_BLK, &numblocks, &fdum, cdum);
    exo_data = VF_Newd(numelems);
    key      = VF_Newi(numelems);
    block    = VF_Newi(numelems);
    proc_id  = VF_Newi(numelems);
    g_index  = VF_Newi(numelems);
    for (offset=0, i=0; i<numblocks; i++) {
      exoerr = ex_get_elem_block(exoid,elemID[i],elem_type,
                                 &nelems[i],&npe[i],&numattr);
      
      exoerr = ex_get_elem_var(exoid, 1, 1, elemID[i], 
                               nelems[i], exo_data);
      for (k=0; k<nelems[i]; k++) {
        key[offset+k] = exo_data[k];
      }
      
      exoerr = ex_get_elem_var(exoid, 1, 2, elemID[i], 
                               nelems[i], exo_data);
      for (k=0; k<nelems[i]; k++) {
        block[offset+k] = exo_data[k];
      }
      
      exoerr = ex_get_elem_var(exoid, 1, 3, elemID[i], 
                               nelems[i], exo_data);
      for (k=0; k<nelems[i]; k++) {
        proc_id[offset+k] = exo_data[k];
      }
      
      exoerr = ex_get_elem_var(exoid, 1, 4, elemID[i], 
                               nelems[i], exo_data);
      for (k=0; k<nelems[i]; k++) {
        g_index[offset+k] = exo_data[k];
      }
      
      offset += nelems[i];
    }
  }
  VF_BroadcastInt(&numelems,1,0);
  VF_BroadcastInt(&numblocks,1,0);
  if (VFLIB_Rank!=0) {
    key     = VF_Newi(numelems);
    block   = VF_Newi(numelems);
    proc_id = VF_Newi(numelems);
    g_index = VF_Newi(numelems);
  }
  index  = VF_Newi(numelems);
  buffer = VF_Newd(numelems);
  VF_BroadcastInt(key,numelems,0);
  VF_BroadcastInt(block,numelems,0);
  VF_BroadcastInt(proc_id,numelems,0);
  VF_BroadcastInt(g_index,numelems,0);
  for (i=0; i<numblocks; i++) {
    for (j=VFLIB_NelemVarAux-2; j<VFLIB_NelemVarAux; j++) {
      for (n=0, k=0; k<numelems; k++) {
        if (block[k]==i && proc_id[k]==VFLIB_Rank) {
          switch (j) {
          case 8:
            value = radj_g[g_index[k]];
            break;
          case 9:
            value = radq_g[g_index[k]];
            break;
          }
          if (VFLIB_Rank==0) {
            exo_data[key[k]] = value;
          } else {
            buffer[n] = value;
            index[n]  = key[k];
            n++;
          }
        }
      }
      /*VF_GatherExodusData(exo_data, buffer, index, n, 2);*/
      VF_GatherExodusData(exo_data, buffer, index, n, 1);
      if (VFLIB_Rank==0) {
        exoerr = ex_put_elem_var(exoid, step, j+1, elemID[i], 
                                 nelems[i], exo_data);
      }
    }
  }
  if (VFLIB_Rank==0) {
    ex_update(exoid);
    VF_Free(exo_data);
  }
  VF_Free(index);
  VF_Free(buffer);
  VF_Free(key);
  VF_Free(block);
  VF_Free(proc_id);
  VF_Free(g_index);
  VF_Free(radq_g);
  VF_Free(radj_g);
  if (exoerr) exoerr = 0;
}

void VF_WriteExodusClose(int exoid)
{
  ex_close(exoid);
}

#else

void VF_WriteExodus(char *filename)  
{
  fprintf(stderr,"WARNING: VF_WriteExodus(%s) not supported\n",filename);
}

int VF_WriteExodusAux1(char *filename)  
{
  fprintf(stderr,"WARNING: VF_WriteExodusAux1(%s) not supported\n",filename);
  return -1;
}

void VF_WriteExodusAux2(int encl, int exoid, int step, double time, double radq[], double radj[])  
{
  fprintf(stderr,"WARNING: VF_WriteExodusAux2() not supported\n");
}

void VF_WriteExodusClose(int exoid)  
{
  fprintf(stderr,"WARNING: VF_WriteExodusClose() not supported\n");
}

#endif
