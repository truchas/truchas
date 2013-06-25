#ifdef VF_MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#undef TRUE
#include "exodusII.h"

int main(int argc,char **argv)
{
  char     file0[132], file1[132], title[132], *cnames[3];
  int      i, j, k, m, n, nss, found;
  int      ndim, code_ws, disk_ws, exoid, exoerr;
  int      nnodes0, nelems0, nblocks0, nns0, nss0;
  int      nnodes1, nelems1, nblocks1, nns1, nss1;
  int      nsides, nfactors, nbar, ntri, nquad;
  int      *c, *conn, *ss_id, *id, *nnodes, *nlist, *tag, *nptr;
  int      nelems[3], npe[3], elemID[] = {100, 200, 300};
  float    version;
  double   *x0, *y0, *z0;
  double   *x1, *y1, *z1;

#ifndef WIN32

  if (argc==1) {
    printf("\n");
    printf("usage:\tcmd sideset_id_list in_file out_file\n");
    printf("\n");
    exit(2);
  } else {
    nss = argc-3;
    id  = (int*)malloc(nss*sizeof(int));
    for (i=0; i<nss; i++) {
      id[i] = atoi(argv[i+1]);
    }
    strcpy(file0,argv[argc-2]);
    strcpy(file1,argv[argc-1]);
  }

#else

  strcpy(file0,"..\\tests\\SphereInSphereOctant_Hemicube\\full_mesh.g");
  strcpy(file1,"..\\tests\\SphereInSphereOctant_Hemicube\\mesh.g");
  nss = 2;
  id = (int*)malloc(nss*sizeof(int));
  id[0] = 1;
  id[1] = 2;

#endif

  code_ws = sizeof(double);
  disk_ws = 0;
  exoid   = ex_open(file0, EX_READ, &code_ws, &disk_ws, &version);
  exoerr  = ex_get_init(exoid,title,&ndim,&nnodes0,
                        &nelems0,&nblocks0,&nns0, &nss0);

  /* READ THE NODES COORDINATES */
  x0 = (double *)malloc(nnodes0*sizeof(double));
  y0 = (double *)malloc(nnodes0*sizeof(double));
  z0 = (double *)malloc(nnodes0*sizeof(double));
  for (i=0; i<nnodes0; i++) {
    x0[i] = y0[i] = z0[i] = 0.0;
  }
  exoerr = ex_get_coord(exoid,x0,y0,z0);

  ss_id = (int *)malloc(nss0*sizeof(int));
  ex_get_side_set_ids(exoid, ss_id);
  for (i=0; i<nss; i++) {
    found = 0;
    for (j=0; j<nss0; j++) {
      if (id[i]==ss_id[j]) {
        found = 1;
        break;
      }
    }
    if (!found) {
      printf("Side Set ID %d was not found in the input mesh\n",id[i]);
      exit(1);
    }
  }
  n = 0;
  for (nelems1=0, i=0; i<nss; i++) {
    ex_get_side_set_param (exoid, id[i], &nsides, &nfactors);
    nelems1 += nsides;
    n = nsides>n?nsides:n;
  }
  c      = (int*)malloc(4*nelems1*sizeof(int));
  nnodes = (int*)malloc(n*sizeof(int));
  nlist  = (int*)malloc(4*n*sizeof(int));
  for (i=0; i<4*nelems1; i++) c[i] = 0;
  for (m=0, i=0; i<nss; i++) {
    ex_get_side_set_param (exoid, id[i], &nsides, &nfactors);
    ex_get_side_set_node_list (exoid, id[i], nnodes, nlist);
    for (nptr=nlist, j=0; j<nsides; j++) {
      for (k=0; k<nnodes[j]; k++) {
        c[4*m+k] = *nptr++;
      }
      m++;
    }
  }
  ex_close(exoid);

  nnodes1 = 0;
  tag     = (int*)malloc(nnodes0*sizeof(int));
  for (i=0; i<nnodes0; i++) tag[i]=0;
  for (i=0; i<nelems1; i++) {
    for (j=0; j<4; j++) {
      if (c[4*i+j]>0) {
        if (tag[c[4*i+j]-1]==0) nnodes1++;
        tag[c[4*i+j]-1] = 1;
      }
    }
  }
  x1 = (double *)malloc(nnodes1*sizeof(double));
  y1 = (double *)malloc(nnodes1*sizeof(double));
  z1 = (double *)malloc(nnodes1*sizeof(double));
  for (n=0, i=0; i<nnodes0; i++) {
    if (tag[i]>0) {
      x1[n]  = x0[i];
      y1[n]  = y0[i];
      z1[n]  = z0[i];
      tag[i] = ++n;
    }
  }
  for (i=0; i<nelems1; i++) {
    for (j=0; j<4; j++) {
      k = 4*i+j;
      n = c[k];
      if (n>0) {
        c[k] = tag[n-1];
      }
    }
  }

  nbar     = 0;
  ntri     = 0;
  nquad    = 0;
  nblocks1 = 0;
  for (i=0; i<nelems1; i++) {
    if (c[4*i+2]==0) 
      nbar++;
    else if (c[4*i+3]==0) 
      ntri++;
    else 
      nquad++;
  }
  if (nbar>0) {
    npe[nblocks1] = 2;
    nelems[nblocks1] = nbar;
    nblocks1++;
  }
  if (ntri>0) {
    npe[nblocks1] = 3;
    nelems[nblocks1] = ntri;
    nblocks1++;
  }
  if (nquad>0) {
    npe[nblocks1] = 4;
    nelems[nblocks1] = nquad;
    nblocks1++;
  }
  nns1      = 0;
  nss1      = 0;
  cnames[0] = "x-coord";
  cnames[1] = "y-coord";
  cnames[2] = "z-coord";
  exoid  = ex_create(file1, EX_CLOBBER, &code_ws, &disk_ws);
  exoerr = ex_put_init(exoid,title,ndim,nnodes1,nelems1,nblocks1,nns1,nss1);
  exoerr = ex_put_coord(exoid,x1,y1,z1);
  exoerr = ex_put_coord_names(exoid,cnames);
  for (i=0; i<nblocks1; i++) {
    conn = (int*)malloc(npe[i]*nelems[i]*sizeof(int));
    for (m=0, j=0; j<nelems1; j++) {
      if (c[4*j+2]==0) 
        n = 2;
      else if (c[4*j+3]==0) 
        n = 3;
      else 
        n = 4;
      if (n==npe[i]) {
        for (k=0; k<npe[i]; k++) {
          conn[m++] = c[4*j+k];
        }
      }
    }
    switch (npe[i]) {
    case 2:
      exoerr = ex_put_elem_block(exoid,elemID[i],"BAR",nelems[i],npe[i],0);
      break;
    case 3:
      exoerr = ex_put_elem_block(exoid,elemID[i],"TRIANGLE",nelems[i],npe[i],0);
      break;
    case 4:
      exoerr = ex_put_elem_block(exoid,elemID[i],"SHELL",nelems[i],npe[i],0);
      break;
    }
    exoerr = ex_put_elem_conn(exoid,elemID[i],conn);
    free(conn);
  }
  exit(0);
}
