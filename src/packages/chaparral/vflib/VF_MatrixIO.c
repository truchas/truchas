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
@(#)    $RCSfile: VF_MatrixIO.c,v $
@(#)    $Revision: 1.4 $  $Date: 2005/09/12 06:51:04 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_MatrixIO.c,v $
@(#)
@(#)    DESCRIPTION:  read/write the viewfactors.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#ifdef HAVE_XDR
# include <unistd.h>
# include <rpc/types.h>
# include <rpc/xdr.h>
#endif

#include "vf.h"

int    GetIntValue  (FILE *fp, char *delimiter);
float  GetFloatValue(FILE *fp, char *delimiter);
double GetDoubleValue(FILE *fp, char *delimiter);
#ifdef HAVE_XDR
void  ReadXDRarray_int   (XDR *xdr, char *ptr, int cnt);
void  ReadXDRarray_float (XDR *xdr, char *ptr, int cnt);
void  WriteXDRarray_int  (XDR *xdr, char *ptr, int cnt);
void  WriteXDRarray_float(XDR *xdr, char *ptr, int cnt);
#endif

#ifdef VF_IO_READ

void
VF_AllInfoRead(char *filename)
{
#ifdef HAVE_XDR
  XDR  xdr;
#endif
  FILE     *fp;
  char     magic_string[9],*magic_ptr=magic_string,line[80];
  unsigned int magic_number;
  int      max_surfaces, num_enclosures, ext_scheme;

  max_surfaces   = 0;
  num_enclosures = 0;
  if (VFLIB_Rank==0) {
    ext_scheme = -1;
    if (ext_scheme<0) {
      /*==========================================*/
      /* CHECK TO SEE IF THE FILE IS ASCII FORMAT */
      /*==========================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(line,"cannot open viewfactor file %s",filename);
        VF_Error(line);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("ASCII   ",magic_string) == 0) {
        ext_scheme = ASCII;
      }
      fclose(fp);
    }
    if (ext_scheme<0) {
      /*==================================================*/
      /* CHECK TO SEE IF THE FILE IS NATIVE BINARY FORMAT */
      /*==================================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(line,"cannot open viewfactor file %s",filename);
        VF_Error(line);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("BINARY  ",magic_string) == 0) {
        fread(&magic_number,1,sizeof(int),fp);
        if (magic_number == 0xFF00FF00) ext_scheme = BINARY;
      }
      fclose(fp);
    }
#ifdef HAVE_XDR
    if (ext_scheme<0) {
      /*========================================*/
      /* CHECK TO SEE IF THE FILE IS XDR FORMAT */
      /*========================================*/
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(line,"cannot open viewfactor file %s",filename);
        VF_Error(line);
      }
      xdrstdio_create(&xdr, fp, XDR_DECODE);
      xdr_string(&xdr, &magic_ptr, 8);
      if (strcmp("XDR     ",magic_string) == 0) ext_scheme = XDRFMTVF_;
      xdr_destroy(&xdr);
      fclose(fp);
    }
#endif
    /*===================================*/
    /* CHECK TO SEE IF THE FILE IS VALID */
    /*===================================*/
    if ((fp=fopen(filename,"r")) == NULL) {
      perror("vf dbase");
      sprintf(line,"cannot open viewfactor file %s",filename);
      VF_Error(line);
    }
    switch (ext_scheme) {
    case ASCII:
      fgets(line, 80, fp);
      num_enclosures = GetIntValue(fp,":");
      max_surfaces   = GetIntValue(fp,":");
      break;
    case BINARY:
      fread(magic_string,1,8,fp);
      fread(&magic_number,1,sizeof(int),fp);
      fread(&num_enclosures,sizeof(int),1,fp);
      fread(&max_surfaces,sizeof(int),1,fp);
      break;
#ifdef HAVE_XDR
    case XDRFMTVF_:
      xdrstdio_create(&xdr, fp, XDR_DECODE);
      xdr_string(&xdr, &magic_ptr, 8);
      xdr_u_int(&xdr, &magic_number);
      xdr_int(&xdr, &num_enclosures);
      xdr_int(&xdr, &max_surfaces);
      break;
#endif
    }
    fclose(fp);
  }
  MP_Broadcast_Int(&num_enclosures,1,0,"read_vfdbase_info(): num_enclosures");
  MP_Broadcast_Int(&max_surfaces,1,0,"read_vfdbase_info(): max_surfaces");
  VF_SetNumEnclosures(nenclosures);
  VF_SetMaxSurfaces(max_patches);
}

void
VF_AllMatrixRead(char *filename, int output)
{
#ifdef HAVE_XDR
  XDR  xdr;
#endif
  FILE     *fp;
  char     magic_string[9],*magic_ptr=magic_string;
  char     s[132],line[132],encl_id[132];
  unsigned int magic_number;
  int      i,j,k,i0,i1,id,ext_scheme;
  int      nonzeros,format;
  int      max_surfaces, num_enclosures, *vf_index;
  long     *toc;
  float    *buffer0, *buffer1;
  double   area;
  VFenclosure *e;

  max_surfaces   = vf_get_nsurfaces_max();
  num_enclosures = vf_get_nenclosures();
  VF_GetSPbuffer0_ptr(&buffer0);
  VF_GetSPbuffer1_ptr(&buffer1);
  VF_GetINTbuffer_ptr(&vf_index);
  if (VFLIB_Rank==0) {
    ext_scheme = -1;
    if (ext_scheme<0) {
      /*==========================================*/
      /* CHECK TO SEE IF THE FILE IS ASCII FORMAT */
      /*==========================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("ASCII   ",magic_string) == 0) {
        ext_scheme = VF_ASCII;
      }
      fclose(fp);
    }
    if (ext_scheme<0) {
      /*==================================================*/
      /* CHECK TO SEE IF THE FILE IS NATIVE BINARY FORMAT */
      /*==================================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("BINARY  ",magic_string) == 0) {
        fread(&magic_number,1,sizeof(int),fp);
        if (magic_number == 0xFF00FF00) ext_scheme = VF_BINARY;
      }
      fclose(fp);
    }
#ifdef HAVE_XDR
    if (ext_scheme<0) {
      /*========================================*/
      /* CHECK TO SEE IF THE FILE IS XDR FORMAT */
      /*========================================*/
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      xdrstdio_create(&xdr, fp, XDR_DECODE);
      xdr_string(&xdr, &magic_ptr, 8);
      if (strcmp("XDR     ",magic_string) == 0) ext_scheme = VF_XDRFMT;
      xdr_destroy(&xdr);
      fclose(fp);
    }
#endif
  }
  VF_Broadcast_Int(&ext_scheme,1,0);
  if (ext_scheme<0) {
    printf("Viewfactor file format not supported\n");
    VF_Exit(1);
  }
  VF_StagedIO_Start();
  /*===================================*/
  /* CHECK TO SEE IF THE FILE IS VALID */
  /*===================================*/
  if ((fp=fopen(filename,"r")) == NULL) {
    perror("vf dbase");
    sprintf(s,"cannot open viewfactor file %s",filename);
    VF_Error(s);
  }
  switch (ext_scheme) {
  case VF_ASCII:
    fgets(line, 80, fp);
    num_enclosures = GetIntValue(fp,":");
    max_surfaces   = GetIntValue(fp,":");
    for (i=0; i<num_enclosures; i++) {
      e = VF_GetEnclosure(i);
      fgets(encl_id, 80, fp);
      e->nonblocking   = GetIntValue(fp,":");
      e->partial       = GetIntValue(fp,":");
      e->asink         = GetDoubleValue(fp,":");
      e->npatches_g    = GetIntValue(fp,":");
      e->npatches_l    = GetIntValue(fp,":");
      e->host_npatches = GetIntValue(fp,":");
      e->vf_method     = GetIntValue(fp,":");
      e->sym_method    = GetIntValue(fp,":");
      e->smoothed      = GetIntValue(fp,":");
      VF_InitializeEnclosure(encl);
      for (j=0; j<e->host_npatches; j++) {
        e->host2vflib_map[j] = GetIntValue(fp,":");
      }
      for (j=0; j<e->npatches_l; j++) {
        VF_InitializeSparseArray(&(e->row[j].array0));
        VF_InitializeSparseArray(&(e->row[j].array1));
        e->row[j].local_index   = GetIntValue(fp,":");
        e->row[j].global_index  = GetIntValue(fp,":");
        e->row[j].host_gid      = GetIntValue(fp,":");
        area                    = GetDoubleValue(fp,":");
        for (k=0; k<e->npatches_g; k++) {
          buffer0[k] = GetFloatValue(fp,":");
        }
        VF_LoadMatrixRow(e->row[j].local_index, buffer0, area);
        e->row[j].raw_rowsum    = GetDoubleValue(fp,":");
        e->row[j].sym_rowsum    = GetDoubleValue(fp,":");
        e->row[j].smooth_rowsum = GetDoubleValue(fp,":");
      }
      if (VFLIB_Size>1) {
        for (j=0; j<VFLIB_Size; j++) {
          e->comm_plan.partners[j].proc         = GetIntValue(fp,":");
          e->comm_plan.partners[j].cnt          = GetIntValue(fp,":");
          e->comm_plan.partners[j].offset       = GetIntValue(fp,":");
          e->comm_plan.partners[j].global_index = VF_Newi(e->comm_plan.partners[j].cnt);
          for (k=0; k<e->comm_plan.partners[j].cnt; k++) {
            e->comm_plan.partners[j].global_index[k] = GetIntValue(fp,":");
          }
        }
      }
    }
    break;
  case VF_BINARY:
    fread(magic_string,1,8,fp);
    fread(&magic_number,1,sizeof(int),fp);
    fread(&num_enclosures,sizeof(int),1,fp);
    fread(&max_surfaces,sizeof(int),1,fp);
    toc = newl(num_enclosures, "file toc");
    fread(toc,sizeof(long),num_enclosures,fp);
    if (encl<0) {
      i0 = 0;
      i1 = num_enclosures;
    } else {
      i0 = encl;
      i1 = encl+1;
      fseek(fp, toc[encl], SEEK_SET);
    }
    VF_Free(toc);
    for (i=i0; i<i1; i++) {
      e = VF_GetEnclosure(i);
      fread(&j,sizeof(int),1,fp);
      fread(encl_id,sizeof(char),j,fp);
      fread(&e->nonblocking,sizeof(int),1,fp);
      fread(&e->partial,sizeof(int),1,fp);
      fread(&e->asink,sizeof(double),1,fp);
      fread(&e->npatches_g,sizeof(int),1,fp);
      fread(&e->npatches_l,sizeof(int),1,fp);
      fread(&e->host_npatches,sizeof(int),1,fp);
      fread(&e->vf_method,sizeof(int),1,fp);
      fread(&e->sym_method,sizeof(int),1,fp);
      fread(&e->smoothed,sizeof(int),1,fp);
      VF_InitializeEnclosure(encl);
      fread(e->host2vflib_map,sizeof(int),e->host_npatches,fp);
      for (j=0; j<e->npatches_l; j++) {
        VF_InitializeSparseArray(&(e->row[j].array0));
        VF_InitializeSparseArray(&(e->row[j].array1));
        fread(&(e->row[j].local_index),sizeof(int),1,fp);
        fread(&(e->row[j].global_index),sizeof(int),1,fp);
        fread(&(e->row[j].host_gid),sizeof(int),1,fp);
        fread(&area,sizeof(double),1,fp);
        fread(&nonzeros,sizeof(int),1,fp);
        fread(buffer1,sizeof(float),nonzeros,fp);
        fread(vf_index,sizeof(int),nonzeros,fp);
        for (k=0; k<e->npatches_g; k++) buffer0[k]=0.0;
        for (k=0; k<nonzeros; k++) {
          buffer0[vf_index[k]] = buffer1[k];
        }
        VF_LoadMatrixRow(e->row[j].local_index, buffer0, area);
        fread(&(e->row[j].raw_rowsum),sizeof(double),1,fp);
        fread(&(e->row[j].sym_rowsum),sizeof(double),1,fp);
        fread(&(e->row[j].smooth_rowsum),sizeof(double),1,fp);
      }
      if (VFLIB_Size>1) {
        for (j=0; j<VFLIB_Size; j++) {
          fread(&e->comm_plan.partners[j].proc,sizeof(int),1,fp);
          fread(&e->comm_plan.partners[j].cnt,sizeof(int),1,fp);
          fread(&e->comm_plan.partners[j].offset,sizeof(int),1,fp);
          e->comm_plan.partners[j].global_index = VF_Newi(e->comm_plan.partners[j].cnt);
          fread(e->comm_plan.partners[j].global_index,sizeof(int),e->comm_plan.partners[j].cnt,fp);
        }
      }
    }
    break;
#ifdef HAVE_XDR
  case VF_XDRFMT:
    xdrstdio_create(&xdr, fp, XDR_DECODE);
    xdr_string(&xdr, &magic_ptr, 8);
    xdr_u_int(&xdr, &magic_number);
    xdr_int(&xdr, &num_enclosures);
    xdr_int(&xdr, &max_surfaces);
    toc = newl(num_enclosures, "file toc");
    ReadXDRarray_int (&xdr, (char *)toc, num_enclosures);
    if (encl<0) {
      i0 = 0;
      i1 = num_enclosures;
    } else {
      i0 = encl;
      i1 = encl+1;
      xdr_setpos(&xdr,toc[encl]);
    }
    VF_Free(toc);
    for (i=i0; i<i1; i++) {
      e = VF_GetEnclosure(i);
      xdr_int   (&xdr, &e->nonblocking);
      xdr_int   (&xdr, &e->partial);
      xdr_double(&xdr, &e->asink);
      xdr_int   (&xdr, &e->npatches_g);
      xdr_int   (&xdr, &e->npatches_l);
      xdr_int   (&xdr, &e->host_npatches);
      xdr_int   (&xdr, &e->vf_method);
      xdr_int   (&xdr, &e->sym_method);
      xdr_int   (&xdr, &e->smoothed);
      VF_InitializeEnclosure(encl);
      ReadXDRarray_int(&xdr, (char*)(e->host2vflib_map),e->host_npatches);
      for (i=0; i<e->npatches_l; i++) {
        VF_InitializeSparseArray(&(e->row[i].array0));
        VF_InitializeSparseArray(&(e->row[i].array1));
        xdr_int   (&xdr, &e->row[i].local_index);
        xdr_int   (&xdr, &e->row[i].global_index);
        xdr_int   (&xdr, &e->row[i].host_gid);
        xdr_double(&xdr, &area);
        xdr_int(&xdr, &nonzeros);
        ReadXDRarray_float(&xdr, (char*)(buffer1), nonzeros);
        ReadXDRarray_int  (&xdr, (char*)(vf_index), nonzeros);
        for (k=0; k<e->npatches_g; k++) buffer0[k]=0.0;
        for (k=0; k<nonzeros; k++) {
          buffer0[vf_index[k]] = buffer1[k];
        }
        VF_LoadMatrixRow(e->row[i].local_index, buffer0, area);
        xdr_double(&xdr, &e->row[i].raw_rowsum);
        xdr_double(&xdr, &e->row[i].sym_rowsum);
        xdr_double(&xdr, &e->row[i].smooth_rowsum);
      }
      if (VFLIB_Size>1) {
        for (i=0; i<VFLIB_Size; i++) {
          xdr_int   (&xdr, &e->comm_plan.partners[i].proc);
          xdr_int   (&xdr, &e->comm_plan.partners[i].cnt);
          xdr_int   (&xdr, &e->comm_plan.partners[i].offset);
          e->comm_plan.partners[i].global_index = VF_Newi(e->comm_plan.partners[i].cnt);
          ReadXDRarray_int  (&xdr, (char*)(e->comm_plan.partners[i].global_index), e->comm_plan.partners[i].cnt);
        }
      }
    }
    xdr_destroy(&xdr);
    break;
#endif
  }
  fclose(fp);
  VF_StagedIO_End();  
}

#endif

void
VF_InfoRead(int encl, char *filename)
{
#ifdef HAVE_XDR
  XDR  xdr;
#endif
  FILE     *fp;
  char     magic_string[9];
#ifdef HAVE_XDR
  char *magic_ptr=magic_string;
#endif
  char     s[132],line[132],encl_id[132];
  unsigned int magic_number;
  int      j,ext_scheme;
  VFenclosure *e=VF_GetEnclosure(encl);

  if (VFLIB_Rank==0) {
    ext_scheme = -1;
    if (ext_scheme<0) {
      /*==========================================*/
      /* CHECK TO SEE IF THE FILE IS ASCII FORMAT */
      /*==========================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("ASCII   ",magic_string) == 0) {
        ext_scheme = VF_ASCII;
      }
      fclose(fp);
    }
    if (ext_scheme<0) {
      /*==================================================*/
      /* CHECK TO SEE IF THE FILE IS NATIVE BINARY FORMAT */
      /*==================================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("BINARY  ",magic_string) == 0) {
        fread(&magic_number,1,sizeof(int),fp);
        if (magic_number == 0xFF00FF00) ext_scheme = VF_BINARY;
      }
      fclose(fp);
    }
#ifdef HAVE_XDR
    if (ext_scheme<0) {
      /*========================================*/
      /* CHECK TO SEE IF THE FILE IS XDR FORMAT */
      /*========================================*/
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      xdrstdio_create(&xdr, fp, XDR_DECODE);
      xdr_string(&xdr, &magic_ptr, 8);
      if (strcmp("XDR     ",magic_string) == 0) ext_scheme = VF_XDRFMT;
      xdr_destroy(&xdr);
      fclose(fp);
    }
#endif
  }
  VF_BroadcastInt(&ext_scheme,1,0);
  if (ext_scheme<0) {
    printf("Viewfactor file format not supported\n");
    VF_Exit(1);
  }
  VF_StagedIO_Start();
  /*===================================*/
  /* CHECK TO SEE IF THE FILE IS VALID */
  /*===================================*/
  if ((fp=fopen(filename,"r")) == NULL) {
    perror("vf dbase");
    sprintf(s,"cannot open viewfactor file %s",filename);
    VF_Error(s);
  }
  switch (ext_scheme) {
  case VF_ASCII:
    fgets(line, 80, fp);
    fgets(encl_id, 80, fp);
    e->nonblocking   = GetIntValue(fp,":");
    e->partial       = GetIntValue(fp,":");
    e->asink         = GetDoubleValue(fp,":");
    e->npatches_g    = GetIntValue(fp,":");
    e->npatches_l    = GetIntValue(fp,":");
    e->host_npatches = GetIntValue(fp,":");
    e->vf_method     = GetIntValue(fp,":");
    e->sym_method    = GetIntValue(fp,":");
    e->smoothed      = GetIntValue(fp,":");
    break;
  case VF_BINARY:
    fread(magic_string,1,8,fp);
    fread(&magic_number,1,sizeof(int),fp);
    fread(&j,sizeof(int),1,fp);
    fread(encl_id,sizeof(char),j,fp);
    fread(&e->nonblocking,sizeof(int),1,fp);
    fread(&e->partial,sizeof(int),1,fp);
    fread(&e->asink,sizeof(double),1,fp);
    fread(&e->npatches_g,sizeof(int),1,fp);
    fread(&e->npatches_l,sizeof(int),1,fp);
    fread(&e->host_npatches,sizeof(int),1,fp);
    fread(&e->vf_method,sizeof(int),1,fp);
    fread(&e->sym_method,sizeof(int),1,fp);
    fread(&e->smoothed,sizeof(int),1,fp);
    break;
#ifdef HAVE_XDR
  case VF_XDRFMT:
    xdrstdio_create(&xdr, fp, XDR_DECODE);
    xdr_string(&xdr, &magic_ptr, 8);
    xdr_u_int (&xdr, &magic_number);
    xdr_int   (&xdr, &e->nonblocking);
    xdr_int   (&xdr, &e->partial);
    xdr_double(&xdr, &e->asink);
    xdr_int   (&xdr, &e->npatches_g);
    xdr_int   (&xdr, &e->npatches_l);
    xdr_int   (&xdr, &e->host_npatches);
    xdr_int   (&xdr, &e->vf_method);
    xdr_int   (&xdr, &e->sym_method);
    xdr_int   (&xdr, &e->smoothed);
    xdr_destroy(&xdr);
    break;
#endif
  }
  fclose(fp);
  VF_StagedIO_End();  
}

void
VF_MatrixRead(int encl, char *filename)
{
#ifdef HAVE_XDR
  XDR  xdr;
#endif
  FILE     *fp;
  char     magic_string[9];
#ifdef HAVE_XDR
  char *magic_ptr=magic_string;
#endif
  char     s[132],line[132],encl_id[132];
  unsigned int magic_number;
  int      i,j,k,ext_scheme;
  int      nonzeros,*vf_index;
  float    *buffer0, *buffer1;
  double   area;
  VFenclosure *e=VF_GetEnclosure(encl);

  VF_GetSPbuffer0_ptr(&buffer0);
  VF_GetSPbuffer1_ptr(&buffer1);
  VF_GetINTbuffer_ptr(&vf_index);
  if (VFLIB_Rank==0) {
    ext_scheme = -1;
    if (ext_scheme<0) {
      /*==========================================*/
      /* CHECK TO SEE IF THE FILE IS ASCII FORMAT */
      /*==========================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("ASCII   ",magic_string) == 0) {
        ext_scheme = VF_ASCII;
      }
      fclose(fp);
    }
    if (ext_scheme<0) {
      /*==================================================*/
      /* CHECK TO SEE IF THE FILE IS NATIVE BINARY FORMAT */
      /*==================================================*/
      magic_string[8] = (char)NULL;
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      fread(magic_string,1,8,fp);
      if (strcmp("BINARY  ",magic_string) == 0) {
        fread(&magic_number,1,sizeof(int),fp);
        if (magic_number == 0xFF00FF00) ext_scheme = VF_BINARY;
      }
      fclose(fp);
    }
#ifdef HAVE_XDR
    if (ext_scheme<0) {
      /*========================================*/
      /* CHECK TO SEE IF THE FILE IS XDR FORMAT */
      /*========================================*/
      if ((fp=fopen(filename,"r")) == NULL) {
        perror("vf dbase");
        sprintf(s,"cannot open viewfactor file %s",filename);
        VF_Error(s);
      }
      xdrstdio_create(&xdr, fp, XDR_DECODE);
      xdr_string(&xdr, &magic_ptr, 8);
      if (strcmp("XDR     ",magic_string) == 0) ext_scheme = VF_XDRFMT;
      xdr_destroy(&xdr);
      fclose(fp);
    }
#endif
  }
  VF_BroadcastInt(&ext_scheme,1,0);
  if (ext_scheme<0) {
    printf("Viewfactor file format not supported\n");
    VF_Exit(1);
  }
  VF_StagedIO_Start();
  /*===================================*/
  /* CHECK TO SEE IF THE FILE IS VALID */
  /*===================================*/
  if ((fp=fopen(filename,"r")) == NULL) {
    perror("vf dbase");
    sprintf(s,"cannot open viewfactor file %s",filename);
    VF_Error(s);
  }
  switch (ext_scheme) {
  case VF_ASCII:
    fgets(line, 80, fp);
    fgets(encl_id, 80, fp);
    e->nonblocking   = GetIntValue(fp,":");
    e->partial       = GetIntValue(fp,":");
    e->asink         = GetDoubleValue(fp,":");
    e->npatches_g    = GetIntValue(fp,":");
    e->npatches_l    = GetIntValue(fp,":");
    e->host_npatches = GetIntValue(fp,":");
    e->vf_method     = GetIntValue(fp,":");
    e->sym_method    = GetIntValue(fp,":");
    e->smoothed      = GetIntValue(fp,":");
    VF_InitializeEnclosure(encl);
    for (i=0; i<e->host_npatches; i++) {
      e->host2vflib_map[i] = GetIntValue(fp,":");
    }
    for (i=0; i<e->npatches_l; i++) {
      VF_InitializeSparseArray(&(e->row[i].array0));
      VF_InitializeSparseArray(&(e->row[i].array1));
      e->row[i].local_index   = GetIntValue(fp,":");
      e->row[i].global_index  = GetIntValue(fp,":");
      e->row[i].host_gid      = GetIntValue(fp,":");
      area                    = GetDoubleValue(fp,":");
      for (k=0; k<e->npatches_g; k++) {
        buffer0[k] = GetFloatValue(fp,":");
      }
      VF_LoadMatrixRow(e->row[i].local_index, buffer0, area);
      e->row[i].raw_rowsum    = GetDoubleValue(fp,":");
      e->row[i].sym_rowsum    = GetDoubleValue(fp,":");
      e->row[i].smooth_rowsum = GetDoubleValue(fp,":");
    }
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        e->comm_plan.partners[i].proc         = GetIntValue(fp,":");
        e->comm_plan.partners[i].cnt          = GetIntValue(fp,":");
        e->comm_plan.partners[i].offset       = GetIntValue(fp,":");
        e->comm_plan.partners[i].global_index = VF_Newi(e->comm_plan.partners[i].cnt);
        for (j=0; j<e->comm_plan.partners[i].cnt; j++) {
          e->comm_plan.partners[i].global_index[j] = GetIntValue(fp,":");
        }
      }
    }
    break;
  case VF_BINARY:
    fread(magic_string,1,8,fp);
    fread(&magic_number,1,sizeof(int),fp);
    fread(&j,sizeof(int),1,fp);
    fread(encl_id,sizeof(char),j,fp);
    fread(&e->nonblocking,sizeof(int),1,fp);
    fread(&e->partial,sizeof(int),1,fp);
    fread(&e->asink,sizeof(double),1,fp);
    fread(&e->npatches_g,sizeof(int),1,fp);
    fread(&e->npatches_l,sizeof(int),1,fp);
    fread(&e->host_npatches,sizeof(int),1,fp);
    fread(&e->vf_method,sizeof(int),1,fp);
    fread(&e->sym_method,sizeof(int),1,fp);
    fread(&e->smoothed,sizeof(int),1,fp);
    VF_InitializeEnclosure(encl);
    fread(e->host2vflib_map,sizeof(int),e->host_npatches,fp);
    for (i=0; i<e->npatches_l; i++) {
      VF_InitializeSparseArray(&(e->row[i].array0));
      VF_InitializeSparseArray(&(e->row[i].array1));
      fread(&(e->row[i].local_index),sizeof(int),1,fp);
      fread(&(e->row[i].global_index),sizeof(int),1,fp);
      fread(&(e->row[i].host_gid),sizeof(int),1,fp);
      fread(&area,sizeof(double),1,fp);
      fread(&nonzeros,sizeof(int),1,fp);
      fread(buffer1,sizeof(float),nonzeros,fp);
      fread(vf_index,sizeof(int),nonzeros,fp);
      for (k=0; k<e->npatches_g; k++) buffer0[k]=0.0;
      for (k=0; k<nonzeros; k++) {
        buffer0[vf_index[k]] = buffer1[k];
      }
      VF_LoadMatrixRow(e->row[i].local_index, buffer0, area);
      fread(&(e->row[i].raw_rowsum),sizeof(double),1,fp);
      fread(&(e->row[i].sym_rowsum),sizeof(double),1,fp);
      fread(&(e->row[i].smooth_rowsum),sizeof(double),1,fp);
    }
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        fread(&e->comm_plan.partners[i].proc,sizeof(int),1,fp);
        fread(&e->comm_plan.partners[i].cnt,sizeof(int),1,fp);
        fread(&e->comm_plan.partners[i].offset,sizeof(int),1,fp);
        e->comm_plan.partners[i].global_index = VF_Newi(e->comm_plan.partners[i].cnt);
        fread(e->comm_plan.partners[i].global_index,sizeof(int),e->comm_plan.partners[i].cnt,fp);
      }
    }
    break;
#ifdef HAVE_XDR
  case VF_XDRFMT:
    xdrstdio_create(&xdr, fp, XDR_DECODE);
    xdr_string(&xdr, &magic_ptr, 8);
    xdr_u_int (&xdr, &magic_number);
    xdr_int   (&xdr, &e->nonblocking);
    xdr_int   (&xdr, &e->partial);
    xdr_double(&xdr, &e->asink);
    xdr_int   (&xdr, &e->npatches_g);
    xdr_int   (&xdr, &e->npatches_l);
    xdr_int   (&xdr, &e->host_npatches);
    xdr_int   (&xdr, &e->vf_method);
    xdr_int   (&xdr, &e->sym_method);
    xdr_int   (&xdr, &e->smoothed);
    VF_InitializeEnclosure(encl);
    ReadXDRarray_int(&xdr, (char*)(e->host2vflib_map),e->host_npatches);
    for (i=0; i<e->npatches_l; i++) {
      VF_InitializeSparseArray(&(e->row[i].array0));
      VF_InitializeSparseArray(&(e->row[i].array1));
      xdr_int   (&xdr, &e->row[i].local_index);
      xdr_int   (&xdr, &e->row[i].global_index);
      xdr_int   (&xdr, &e->row[i].host_gid);
      xdr_double(&xdr, &area);
      xdr_int(&xdr, &nonzeros);
      ReadXDRarray_float(&xdr, (char*)(buffer1), nonzeros);
      ReadXDRarray_int  (&xdr, (char*)(vf_index), nonzeros);
      for (k=0; k<e->npatches_g; k++) buffer0[k]=0.0;
      for (k=0; k<nonzeros; k++) {
        buffer0[vf_index[k]] = buffer1[k];
      }
      VF_LoadMatrixRow(e->row[i].local_index, buffer0, area);
      xdr_double(&xdr, &e->row[i].raw_rowsum);
      xdr_double(&xdr, &e->row[i].sym_rowsum);
      xdr_double(&xdr, &e->row[i].smooth_rowsum);
    }
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        xdr_int   (&xdr, &e->comm_plan.partners[i].proc);
        xdr_int   (&xdr, &e->comm_plan.partners[i].cnt);
        xdr_int   (&xdr, &e->comm_plan.partners[i].offset);
        e->comm_plan.partners[i].global_index = VF_Newi(e->comm_plan.partners[i].cnt);
        ReadXDRarray_int  (&xdr, (char*)(e->comm_plan.partners[i].global_index), e->comm_plan.partners[i].cnt);
      }
    }
    xdr_destroy(&xdr);
    break;
#endif
  }
  fclose(fp);
  VF_StagedIO_End();  
}

void
VF_AllMatrixWrite(char *filename, int format) 
{
#ifdef HAVE_XDR
  XDR  xdr;
#endif
  FILE *fp;
  char s[80], *magic_string;
  int  i, j, k, n, nn, rownum;
  long *toc, toc_start;
  float *buffer0;
  int  max_surfaces, num_enclosures, *vf_index;
  unsigned int magic_number=0xFF00FF00;
  VFenclosure *e;

  max_surfaces   = VF_MaxPatches();
  num_enclosures = VF_NumEnclosures();
  toc            = VF_Newl(num_enclosures);
  VF_GetSPbuffer0_ptr(&buffer0);
  VF_GetINTbuffer_ptr(&vf_index);
  VF_StagedIO_Start();
  if ((fp=fopen(filename,"w")) == NULL) {
    perror("vf dbase");
    sprintf(s,"cannot open viewfactor file %s",filename);
    VF_Error(s);
  }
  switch (format) {
  case VF_ASCII:
    fprintf(fp,"ASCII   \n");
    fprintf(fp,"Number of Enclosures: %d\n",num_enclosures);
    fprintf(fp,"Max # of Surfaces:    %d\n",max_surfaces);
    for (i=0; i<num_enclosures; i++) {
      e = VF_GetEnclosure(i);
      fprintf(fp,"Enclosure ID:           %s\n",e->id);
      fprintf(fp,"Blocking Enclosure:     %d\n",e->nonblocking);
      fprintf(fp,"Partial Enclosure:      %d\n",e->partial);
      fprintf(fp,"Partial Enclosure Area: %g\n",e->asink);
      fprintf(fp,"# of Global Patches:    %d\n",e->npatches_g);
      fprintf(fp,"# of Local Patches:     %d\n",e->npatches_l);
      fprintf(fp,"# of Host Patches:      %d\n",e->host_npatches);
      fprintf(fp,"Calculation Method:     %d\n",e->vf_method);
      fprintf(fp,"Symmetry Method:        %d\n",e->sym_method);
      fprintf(fp,"Smoothed:               %d\n",e->smoothed);
      for (j=0; j<e->host_npatches; j++) {
        fprintf(fp,"host2vflib map:  %d\n",e->host2vflib_map[j]);
      }
      for (j=0; j<e->npatches_l; j++) {
        fprintf(fp,"row local index:  %d\n",e->row[j].local_index);
        fprintf(fp,"  global index:   %d\n",e->row[j].global_index);
        fprintf(fp,"  host gid:       %d\n",e->row[j].host_gid);
        fprintf(fp,"  area:           %f\n",e->row[j].area);
        VF_GetMatrixRow(j,&rownum,buffer0);
        for (k=0; k<e->npatches_g; k++) {
          fprintf(fp,"  VF(%d,%d): %f\n",e->row[j].global_index,k,buffer0[k]);
        }
        fprintf(fp,"  raw rowsum:     %f\n",e->row[j].raw_rowsum);
        fprintf(fp,"  sym rowsum:     %f\n",e->row[j].sym_rowsum);
        fprintf(fp,"  sm  rowsum:     %f\n",e->row[j].smooth_rowsum);
      }
      if (VFLIB_Size>1) {
        for (j=0; j<VFLIB_Size; j++) {
          fprintf(fp,"Comm partner proc     = %d",e->comm_plan.partners[j].proc);
          fprintf(fp,"  Comm partner count  = %d",e->comm_plan.partners[j].cnt);
          fprintf(fp,"  Comm partner offset = %d",e->comm_plan.partners[j].offset);
          for (k=0; k<e->comm_plan.partners[j].cnt; k++) {
            fprintf(fp,"  Comm partner global index = %d",e->comm_plan.partners[j].global_index[k]);
          }
        }
      }
    }
    break;
  case VF_BINARY:
    magic_string = "BINARY  ";
    fwrite(magic_string,1,8,fp);
    fwrite(&magic_number,1,sizeof(int),fp);
    fwrite(&num_enclosures,sizeof(int),1,fp);
    fwrite(&max_surfaces,sizeof(int),1,fp);
    toc_start = ftell(fp);
    fwrite(toc,sizeof(long),num_enclosures,fp);
    for (i=0; i<num_enclosures; i++) {
      toc[i] = ftell(fp);
      e      = VF_GetEnclosure(i);
      j      = strlen(e->id);
      fwrite(&j,sizeof(int),1,fp);
      fwrite(e->id,sizeof(char),j,fp);
      fwrite(&e->nonblocking,sizeof(int),1,fp);
      fwrite(&e->partial,sizeof(int),1,fp);
      fwrite(&e->asink,sizeof(double),1,fp);
      fwrite(&e->npatches_g,sizeof(int),1,fp);
      fwrite(&e->npatches_l,sizeof(int),1,fp);
      fwrite(&e->host_npatches,sizeof(int),1,fp);
      fwrite(&e->vf_method,sizeof(int),1,fp);
      fwrite(&e->sym_method,sizeof(int),1,fp);
      fwrite(&e->smoothed,sizeof(int),1,fp);
      fwrite(e->host2vflib_map,sizeof(int),e->host_npatches,fp);
      for (j=0; j<e->npatches_l; j++) {
        fwrite(&e->row[j].local_index,sizeof(int),1,fp);
        fwrite(&e->row[j].global_index,sizeof(int),1,fp);
        fwrite(&e->row[j].host_gid,sizeof(int),1,fp);
        fwrite(&e->row[j].area,sizeof(double),1,fp);
        nn = e->row[j].array0.cnt;
        for (n=0, k=0; k<nn; k++, n++) {
          vf_index[n] = e->row[j].array0.index[k];
          buffer0[n]  = e->row[j].array0.data[k];
        }
        nn = e->row[j].array1.cnt;
        for (k=0; k<nn; k++, n++) {
          vf_index[n] = e->row[j].array1.index[k];
          buffer0[n]  = e->row[j].array1.data[k];
        }
        VF_SortSparseArrayAux(vf_index,buffer0,n);
        fwrite(&n,sizeof(int),1,fp);
        fwrite(buffer0,sizeof(float),n,fp);
        fwrite(vf_index,sizeof(int),n,fp);
        fwrite(&e->row[j].raw_rowsum,sizeof(double),1,fp);
        fwrite(&e->row[j].sym_rowsum,sizeof(double),1,fp);
        fwrite(&e->row[j].smooth_rowsum,sizeof(double),1,fp);
      }
      if (VFLIB_Size>1) {
        for (j=0; j<VFLIB_Size; j++) {
          fwrite(&e->comm_plan.partners[j].proc,sizeof(int),1,fp);
          fwrite(&e->comm_plan.partners[j].cnt,sizeof(int),1,fp);
          fwrite(&e->comm_plan.partners[j].offset,sizeof(int),1,fp);
          fwrite(e->comm_plan.partners[j].global_index,sizeof(int),e->comm_plan.partners[j].cnt,fp);
        }
      }
    }
    fseek(fp, toc_start, SEEK_SET);
    fwrite(toc,sizeof(long),num_enclosures,fp);
    fseek(fp, 0, SEEK_END);
    break;
#ifdef HAVE_XDR
  case VF_XDRFMT:
    magic_string = "XDR     ";
    xdrstdio_create(&xdr, fp, XDR_ENCODE);
    xdr_string(&xdr, &magic_string, 8);
    xdr_u_int(&xdr, &magic_number);
    xdr_int(&xdr, &num_enclosures);
    xdr_int(&xdr, &max_surfaces);
    toc_start = xdr_getpos(&xdr);
    WriteXDRarray_int(&xdr, (char *)toc, num_enclosures);
    for (i=0; i<num_enclosures; i++) {
      toc[i] = xdr_getpos(&xdr);
      e      = VF_GetEnclosure(i);
      xdr_int   (&xdr, &e->nonblocking);
      xdr_int   (&xdr, &e->partial);
      xdr_double(&xdr, &e->asink);
      xdr_int   (&xdr, &e->npatches_g);
      xdr_int   (&xdr, &e->npatches_l);
      xdr_int   (&xdr, &e->host_npatches);
      xdr_int   (&xdr, &e->vf_method);
      xdr_int   (&xdr, &e->sym_method);
      xdr_int   (&xdr, &e->smoothed);
      WriteXDRarray_int(&xdr, (char*)(e->host2vflib_map),e->host_npatches);
      for (j=0; j<e->npatches_l; j++) {
        xdr_int   (&xdr, &e->row[j].local_index);
        xdr_int   (&xdr, &e->row[j].global_index);
        xdr_int   (&xdr, &e->row[j].host_gid);
        xdr_double(&xdr, &e->row[j].area);
        nn = e->row[j].array0.cnt;
        for (n=0, k=0; k<nn; k++, n++) {
          vf_index[n] = e->row[j].array0.index[k];
          buffer0[n]  = e->row[j].array0.data[k];
        }
        nn = e->row[j].array1.cnt;
        for (k=0; k<nn; k++, n++) {
          vf_index[n] = e->row[j].array1.index[k];
          buffer0[n]  = e->row[j].array1.data[k];
        }
        VF_SortSparseArrayAux(vf_index,buffer0,n);
        xdr_int(&xdr, &n);
        WriteXDRarray_float(&xdr, (char*)(buffer0), n);
        WriteXDRarray_int  (&xdr, (char*)(vf_index), n);
        xdr_double(&xdr, &e->row[j].raw_rowsum);
        xdr_double(&xdr, &e->row[j].sym_rowsum);
        xdr_double(&xdr, &e->row[j].smooth_rowsum);
      }
      if (VFLIB_Size>1) {
        for (j=0; j<VFLIB_Size; j++) {
          xdr_int   (&xdr, &e->comm_plan.partners[j].proc);
          xdr_int   (&xdr, &e->comm_plan.partners[j].cnt);
          xdr_int   (&xdr, &e->comm_plan.partners[j].offset);
          WriteXDRarray_int  (&xdr, (char*)(e->comm_plan.partners[j].global_index), e->comm_plan.partners[j].cnt);
        }
      }
    }
    xdr_setpos(&xdr, toc_start);
    WriteXDRarray_int(&xdr, (char*)toc, num_enclosures);
    xdr_destroy(&xdr);
    break;
#endif
  default:
    printf("Storage format %d unknown\n",format);
    break;
  }
  fclose(fp);
  VF_Free(toc);
  VF_StagedIO_End();
}

void
VF_MatrixWrite(int encl, char *filename, int format) 
{
#ifdef HAVE_XDR
  XDR  xdr;
#endif
  FILE  *fp;
  char  s[80], *magic_string;
  int   i, j, k, n, nn, rownum;
  float *buffer0;
  int   *vf_index;
  unsigned int magic_number=0xFF00FF00;
  VFenclosure *e=VF_GetEnclosure(encl);

  VF_StagedIO_Start();
  VF_GetSPbuffer0_ptr(&buffer0);
  VF_GetINTbuffer_ptr(&vf_index);
  if ((fp=fopen(filename,"w")) == NULL) {
    perror("vf dbase");
    sprintf(s,"cannot open viewfactor file %s",filename);
    VF_Error(s);
  }
  switch (format) {
  case VF_ASCII:
    fprintf(fp,"ASCII   \n");
    fprintf(fp,"Enclosure ID:           %s\n",e->id);
    fprintf(fp,"Blocking Enclosure:     %d\n",e->nonblocking);
    fprintf(fp,"Partial Enclosure:      %d\n",e->partial);
    fprintf(fp,"Partial Enclosure Area: %g\n",e->asink);
    fprintf(fp,"# of Global Patches:    %d\n",e->npatches_g);
    fprintf(fp,"# of Local Patches:     %d\n",e->npatches_l);
    fprintf(fp,"# of Host Patches:      %d\n",e->host_npatches);
    fprintf(fp,"Calculation Method:     %d\n",e->vf_method);
    fprintf(fp,"Symmetry Method:        %d\n",e->sym_method);
    fprintf(fp,"Smoothed:               %d\n",e->smoothed);
    for (i=0; i<e->host_npatches; i++) {
      fprintf(fp,"host2vflib map:  %d\n",e->host2vflib_map[i]);
    }
    for (j=0; j<e->npatches_l; j++) {
      fprintf(fp,"row local index:  %d\n",e->row[j].local_index);
      fprintf(fp,"  global index:   %d\n",e->row[j].global_index);
      fprintf(fp,"  host gid:       %d\n",e->row[j].host_gid);
      fprintf(fp,"  area:           %f\n",e->row[j].area);
      VF_GetMatrixRow(j,&rownum,buffer0);
      for (k=0; k<e->npatches_g; k++) {
        fprintf(fp,"  VF(%d,%d): %f\n",e->row[j].global_index,k,buffer0[k]);
      }
      fprintf(fp,"  raw rowsum:     %f\n",e->row[j].raw_rowsum);
      fprintf(fp,"  sym rowsum:     %f\n",e->row[j].sym_rowsum);
      fprintf(fp,"  sm  rowsum:     %f\n",e->row[j].smooth_rowsum);
    }
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        fprintf(fp,"Comm partner proc     = %d",e->comm_plan.partners[i].proc);
        fprintf(fp,"  Comm partner count  = %d",e->comm_plan.partners[i].cnt);
        fprintf(fp,"  Comm partner offset = %d",e->comm_plan.partners[i].offset);
        for (j=0; j<e->comm_plan.partners[i].cnt; j++) {
          fprintf(fp,"  Comm partner global index = %d",e->comm_plan.partners[i].global_index[j]);
        }
      }
    }
    break;
  case VF_BINARY:
    j            = strlen(e->id);
    magic_string = "BINARY  ";
    fwrite(magic_string,1,8,fp);
    fwrite(&magic_number,1,sizeof(int),fp);
    fwrite(&j,sizeof(int),1,fp);
    fwrite(e->id,sizeof(char),j,fp);
    fwrite(&e->nonblocking,sizeof(int),1,fp);
    fwrite(&e->partial,sizeof(int),1,fp);
    fwrite(&e->asink,sizeof(double),1,fp);
    fwrite(&e->npatches_g,sizeof(int),1,fp);
    fwrite(&e->npatches_l,sizeof(int),1,fp);
    fwrite(&e->host_npatches,sizeof(int),1,fp);
    fwrite(&e->vf_method,sizeof(int),1,fp);
    fwrite(&e->sym_method,sizeof(int),1,fp);
    fwrite(&e->smoothed,sizeof(int),1,fp);
    fwrite(e->host2vflib_map,sizeof(int),e->host_npatches,fp);
    for (j=0; j<e->npatches_l; j++) {
      fwrite(&e->row[j].local_index,sizeof(int),1,fp);
      fwrite(&e->row[j].global_index,sizeof(int),1,fp);
      fwrite(&e->row[j].host_gid,sizeof(int),1,fp);
      fwrite(&e->row[j].area,sizeof(double),1,fp);
      nn = e->row[j].array0.cnt;
      for (n=0, k=0; k<nn; k++, n++) {
        vf_index[n] = e->row[j].array0.index[k];
        buffer0[n]  = e->row[j].array0.data[k];
      }
      nn = e->row[j].array1.cnt;
      for (k=0; k<nn; k++, n++) {
        vf_index[n] = e->row[j].array1.index[k];
        buffer0[n]  = e->row[j].array1.data[k];
      }
      VF_SortSparseArrayAux(vf_index,buffer0,n);
      fwrite(&n,sizeof(int),1,fp);
      fwrite(buffer0,sizeof(float),n,fp);
      fwrite(vf_index,sizeof(int),n,fp);
      fwrite(&e->row[j].raw_rowsum,sizeof(double),1,fp);
      fwrite(&e->row[j].sym_rowsum,sizeof(double),1,fp);
      fwrite(&e->row[j].smooth_rowsum,sizeof(double),1,fp);
    }
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        fwrite(&e->comm_plan.partners[i].proc,sizeof(int),1,fp);
        fwrite(&e->comm_plan.partners[i].cnt,sizeof(int),1,fp);
        fwrite(&e->comm_plan.partners[i].offset,sizeof(int),1,fp);
        fwrite(e->comm_plan.partners[i].global_index,sizeof(int),e->comm_plan.partners[i].cnt,fp);
      }
    }
    break;
#ifdef HAVE_XDR
  case VF_XDRFMT:
    magic_string = "XDR     ";
    xdrstdio_create(&xdr, fp, XDR_ENCODE);
    xdr_string(&xdr, &magic_string, 8);
    xdr_u_int (&xdr, &magic_number);
    xdr_int   (&xdr, &e->nonblocking);
    xdr_int   (&xdr, &e->partial);
    xdr_double(&xdr, &e->asink);
    xdr_int   (&xdr, &e->npatches_g);
    xdr_int   (&xdr, &e->npatches_l);
    xdr_int   (&xdr, &e->host_npatches);
    xdr_int   (&xdr, &e->vf_method);
    xdr_int   (&xdr, &e->sym_method);
    xdr_int   (&xdr, &e->smoothed);
    WriteXDRarray_int(&xdr, (char*)(e->host2vflib_map),e->host_npatches);
    for (j=0; j<e->npatches_l; j++) {
      xdr_int   (&xdr, &e->row[j].local_index);
      xdr_int   (&xdr, &e->row[j].global_index);
      xdr_int   (&xdr, &e->row[j].host_gid);
      xdr_double(&xdr, &e->row[j].area);
      nn = e->row[j].array0.cnt;
      for (n=0, k=0; k<nn; k++, n++) {
        vf_index[n] = e->row[j].array0.index[k];
        buffer0[n]  = e->row[j].array0.data[k];
      }
      nn = e->row[j].array1.cnt;
      for (k=0; k<nn; k++, n++) {
        vf_index[n] = e->row[j].array1.index[k];
        buffer0[n]  = e->row[j].array1.data[k];
      }
      VF_SortSparseArrayAux(vf_index,buffer0,n);
      xdr_int(&xdr, &n);
      WriteXDRarray_float(&xdr, (char*)(buffer0), n);
      WriteXDRarray_int  (&xdr, (char*)(vf_index), n);
      xdr_double(&xdr, &e->row[j].raw_rowsum);
      xdr_double(&xdr, &e->row[j].sym_rowsum);
      xdr_double(&xdr, &e->row[j].smooth_rowsum);
    }
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        xdr_int   (&xdr, &e->comm_plan.partners[i].proc);
        xdr_int   (&xdr, &e->comm_plan.partners[i].cnt);
        xdr_int   (&xdr, &e->comm_plan.partners[i].offset);
        WriteXDRarray_int  (&xdr, (char*)(e->comm_plan.partners[i].global_index), e->comm_plan.partners[i].cnt);
      }
    }
    xdr_destroy(&xdr);
    break;
#endif
  default:
    printf("Storage format %d unknown\n",format);
    break;
  }
  fclose(fp);
  VF_StagedIO_End();
}

int
GetIntValue(FILE *fp, char *delimiter)
{
  char line[80], *token;
  int  length;
    
  fgets(line, 80, fp);
  length       = strlen(line)-1;
  line[length] = (char)NULL;
  token        = strtok(line,delimiter);
  token        = strtok(NULL," ");
  return (atoi(token));
}

float
GetFloatValue(FILE *fp, char *delimiter)
{
  char line[80], *token;
  int  length;
    
  fgets(line, 80, fp);
  length       = strlen(line)-1;
  line[length] = (char)NULL;
  token        = strtok(line,delimiter);
  token        = strtok(NULL," ");
  return ((float)atof(token));
}

double
GetDoubleValue(FILE *fp, char *delimiter)
{
  char line[80], *token;
  int  length;
    
  fgets(line, 80, fp);
  length       = strlen(line)-1;
  line[length] = (char)NULL;
  token        = strtok(line,delimiter);
  token        = strtok(NULL," ");
  return (atof(token));
}

#ifdef HAVE_XDR

void
ReadXDRarray_int(XDR *xdr, char *ptr, int cnt)
{
  unsigned int size=cnt;
  xdr_array(xdr, &ptr, &size, size, sizeof(int), (xdrproc_t)xdr_int);
}

void
ReadXDRarray_float(XDR *xdr, char *ptr, int cnt)
{
  unsigned int size=cnt;
  xdr_array(xdr, &ptr, &size, size, sizeof(float), (xdrproc_t)xdr_float);
}

void
WriteXDRarray_int(XDR *xdr, char *ptr, int cnt)
{
  unsigned int size=cnt;
  xdr_array(xdr, &ptr, &size, size, sizeof(int), (xdrproc_t)xdr_int);
}

void
WriteXDRarray_float(XDR *xdr, char *ptr, int cnt)
{
  unsigned int size=cnt;
  xdr_array(xdr, &ptr, &size, size, sizeof(float), (xdrproc_t)xdr_float);
}

#endif
