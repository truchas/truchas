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
@(#)    $RCSfile: VF_MatrixUtils.c,v $
@(#)    $Revision: 1.6.2.2 $  $Date: 2006/05/10 18:15:34 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_MatrixUtils.c,v $
@(#)
@(#)    DESCRIPTION:  Utilities to store/retrieve/manipulate viewfactors.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "vf.h"

static int    max_patches=0;
static int    *vf_index=NULL;
static float  *SPbuffer0=NULL;
static float  *SPbuffer1=NULL;
static float  *SPbuffer2=NULL;
static double *DPbuffer0=NULL;
static double *DPbuffer1=NULL;
static double *DPbuffer2=NULL;
static double *DPbuffer3=NULL;
static int    VF_buffers_initialized=FALSE;

int 
VF_MaxPatches(void)
{
  return max_patches;
}

void
VF_InitializeBuffers(int npatches)
{
  int doit=FALSE;

  if (!VF_buffers_initialized) {
    VF_buffers_initialized = TRUE;
    doit                   = TRUE;
  } else {
    if (npatches>max_patches) {
      VF_Free(SPbuffer0);
      VF_Free(SPbuffer1);
      VF_Free(SPbuffer2);
      VF_Free(DPbuffer0);
      VF_Free(DPbuffer1);
      VF_Free(DPbuffer2);
      VF_Free(DPbuffer3);
      VF_Free(vf_index);
      doit = TRUE;
    }
  }
  if (doit) {
    max_patches = MAX(npatches, 1000);
    SPbuffer0   = VF_Newf(max_patches);
    SPbuffer1   = VF_Newf(max_patches);
    SPbuffer2   = VF_Newf(max_patches);
    DPbuffer0   = VF_Newd(max_patches);
    DPbuffer1   = VF_Newd(max_patches);
    DPbuffer2   = VF_Newd(max_patches);
    DPbuffer3   = VF_Newd(max_patches);
    vf_index    = VF_Newi(max_patches);
  }
}

void
VF_FreeBuffers()
{
  VF_Free(SPbuffer0);
  VF_Free(SPbuffer1);
  VF_Free(SPbuffer2);
  VF_Free(DPbuffer0);
  VF_Free(DPbuffer1);
  VF_Free(DPbuffer2);
  VF_Free(DPbuffer3);
  VF_Free(vf_index);
  VF_buffers_initialized = FALSE;
  max_patches            = 0;
}

void
VF_LoadMatrixRow(int index, float vf_buffer[], double area)
{
  int i, ii, knt;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  VFrow       *row=&(enclosure->row[index]);

  for (knt=0, i=0; i<enclosure->npatches_g; i++) {
    if (vf_buffer[i] != 0.0) knt++;
  }
  if (knt>0) {
    VF_AllocateSparseArray(&(row->array0), knt);
    for (ii=0, i=0; i<enclosure->npatches_g; i++) {
      if (vf_buffer[i]!=0.0) {
        row->array0.data[ii]  = vf_buffer[i];
        row->array0.index[ii] = i;
        ii++;
      }
    }
    row->array0.cnt = knt;
  }
  row->diagonal = vf_buffer[row->global_index];
  row->area     = area;
}

void
VF_UpdateMatrixRow(int n, float vf_buffer[])
{
  int i, knt, nonzeros;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  for (nonzeros=0, i=0; i<enclosure->npatches_g; i++) {
    if (vf_buffer[i]!=0.0) {
      nonzeros++;
    }
  }
  
#define VF_USE_TWO_ARRAY 
#ifdef VF_USE_ONE_ARRAY 
  
  if (enclosure->row[n].array0.size==0) {
    VF_AllocateSparseArray(&(enclosure->row[n].array0), nonzeros);
  } else if (enclosure->row[n].array0.size<nonzeros) {
    VF_ReallocateSparseArray(&(enclosure->row[n].array0), nonzeros);
  }
  for (knt=0, i=0; i<enclosure->npatches_g; i++) {
    if (vf_buffer[i]!=0.0) {
      enclosure->row[n].array0.data[knt]  = vf_buffer[i];
      enclosure->row[n].array0.index[knt] = i;
      knt++;
    }
  }
  enclosure->row[n].array0.cnt = knt;

#else
    
  if (enclosure->row[n].array0.size==0) {
    VF_AllocateSparseArray(&(enclosure->row[n].array0), nonzeros);
    for (knt=0, i=0; i<enclosure->npatches_g; i++) {
      if (vf_buffer[i]!=0.0) {
        enclosure->row[n].array0.data[knt]  = vf_buffer[i];
        enclosure->row[n].array0.index[knt] = i;
        knt++;
      }
    }
    enclosure->row[n].array0.cnt = knt;
  } else if (nonzeros<=enclosure->row[n].array0.size) {
    enclosure->row[n].array0.cnt = nonzeros;
    for (knt=0, i=0; i<enclosure->npatches_g; i++) {
      if (vf_buffer[i]!=0.0) {
        enclosure->row[n].array0.data[knt]  = vf_buffer[i];
        enclosure->row[n].array0.index[knt] = i;
        knt++;
      }
    }
    enclosure->row[n].array0.cnt = knt;
    VF_FreeSparseArray(&(enclosure->row[n].array1));
  } else {
    nonzeros -= enclosure->row[n].array0.cnt;
    if (enclosure->row[n].array1.size == 0) {
      VF_AllocateSparseArray(&(enclosure->row[n].array1), nonzeros);
    } else if (nonzeros>enclosure->row[n].array1.size) {
      VF_ReallocateSparseArray(&(enclosure->row[n].array1), nonzeros);
    }
    for (knt=0, i=0; i<enclosure->npatches_g; i++) {
      if (vf_buffer[i]!=0.0) {
        if (knt==enclosure->row[n].array0.size) break;
        enclosure->row[n].array0.data[knt]  = vf_buffer[i];
        enclosure->row[n].array0.index[knt] = i;
        knt++;
      }
    }
    enclosure->row[n].array0.cnt = knt;
    for (knt=0; i<enclosure->npatches_g; i++) {
      if (vf_buffer[i]!=0.0) {
        enclosure->row[n].array1.data[knt]  = vf_buffer[i];
        enclosure->row[n].array1.index[knt] = i;
        knt++;
      }
    }
    enclosure->row[n].array1.cnt = knt;
  }

#endif

  enclosure->row[n].diagonal = vf_buffer[enclosure->row[n].global_index];
}

void
VF_SetRawRowsum(int index)
{
  double rowsum=0.0;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  rowsum += VF_SumSparseArray(&(enclosure->row[index].array0));
  rowsum += VF_SumSparseArray(&(enclosure->row[index].array1));
  enclosure->row[index].raw_rowsum    = rowsum;
  enclosure->row[index].sym_rowsum    = rowsum;
  enclosure->row[index].smooth_rowsum = rowsum;
}

void
VF_SetSymmetricRowsum(int index)
{
  double rowsum=0.0;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  rowsum += VF_SumSparseArray(&(enclosure->row[index].array0));
  rowsum += VF_SumSparseArray(&(enclosure->row[index].array1));
  enclosure->row[index].sym_rowsum    = rowsum;
  enclosure->row[index].smooth_rowsum = rowsum;
}

void
VF_SetSmoothedRowsum(int index)
{
  double rowsum=0.0;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  rowsum += VF_SumSparseArray(&(enclosure->row[index].array0));
  rowsum += VF_SumSparseArray(&(enclosure->row[index].array1));
  enclosure->row[index].smooth_rowsum = rowsum;
}

void
VF_GetMatrixRow(int index, int *rownum, float vf_buffer[])
{
  int         i;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  VFrow       *row=&(enclosure->row[index]);

  *rownum = row->global_index;
  for (i=0; i<enclosure->npatches_g; i++) vf_buffer[i] = 0.0;
  VF_ExpandSparseArray(&(row->array0), vf_buffer);
  VF_ExpandSparseArray(&(row->array1), vf_buffer);
}

void
VF_GetMatrixCol(int column, int root, float vf_buffer[])
{
  int   i;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  for (i=0; i<enclosure->npatches_g; i++) vf_buffer[i] = 0.0;
  for (i=0; i<enclosure->npatches_l; i++) {
    vf_buffer[i] = VF_GetMatrixEntry(i, column);
  }
  VF_GatherVFs(vf_buffer, root);
}

float
VF_GetMatrixEntry(int row, int col)
{
  int   i;
  float vf=0.0;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  i = VF_LocateIndexInSortedArray(col, enclosure->row[row].array0.index, 
                                       enclosure->row[row].array0.cnt);
  if (i!=VF_INVALID) {
    vf = enclosure->row[row].array0.data[i];
  } else {
    i = VF_LocateIndexInSortedArray(col, enclosure->row[row].array1.index, 
                                         enclosure->row[row].array1.cnt);
    if (i!=VF_INVALID) {
      vf = enclosure->row[row].array1.data[i];
    }
  }
  return vf;
}

int
VF_GetMatrixRowGID(int index)
{
  VFenclosure *enclosure=VF_CurrentEnclosure();
  return(enclosure->row[index].global_index);
}

float
VF_GetMatrixRowLastEntry(int n)
{
  int   last,index,found=0;
  float vf;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (!found && enclosure->row[n].array0.cnt>0) {
    last  = enclosure->row[n].array0.cnt-1;
    index = enclosure->row[n].array0.index[last];
    if (index==enclosure->npatches_g-1) {
      vf    = enclosure->row[n].array0.data[last];
      found = 1;
    }
  }
  if (!found && enclosure->row[n].array1.cnt>0) {
    last  = enclosure->row[n].array1.cnt-1;
    index = enclosure->row[n].array1.index[last];
    if (index==enclosure->npatches_g-1) {
      vf    = enclosure->row[n].array1.data[last];
      found = 1;
    }
  }
  if (!found) vf = 0.0;
  return (vf);
}

double 
VF_GetRawRowsum(int index)
{
  VFenclosure *enclosure=VF_CurrentEnclosure();
  return (enclosure->row[index].raw_rowsum);
}

double
VF_GetSmoothedRowsum(int index)
{
  VFenclosure *enclosure=VF_CurrentEnclosure();
  return (enclosure->row[index].smooth_rowsum);
}

double
VF_GetMatrixRowArea(int index)
{
  VFenclosure *enclosure=VF_CurrentEnclosure();
  return (enclosure->row[index].area);
}

void VF_GetMatrixAreas(double *areas)
{
  int   i, n;
  VFrow *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  row = enclosure->row; 

  n = enclosure->npatches_g; 
  for (i=0; i<n; i++) {
    areas[i] = 0.0;
  }
  n = enclosure->npatches_l;
  for (i=0; i<n; i++, row++) {
    areas[row->global_index] = row->area;
  }
  VF_ExchangeDouble(areas);
}

void 
VF_RowsumStats(double *mn0, double *mx0, double *m0, double *s0,
               double *mn1, double *mx1, double *m1, double *s1,
               double *mn2, double *mx2, double *m2, double *s2)
{
  int    i, nsurf_l;
  double err0, err1, err2, nsurf_g;
  double min0=MAXDOUBLE, max0=-MAXDOUBLE, sum0=0.0, sd0=0.0;
  double min1=MAXDOUBLE, max1=-MAXDOUBLE, sum1=0.0, sd1=0.0;
  double min2=MAXDOUBLE, max2=-MAXDOUBLE, sum2=0.0, sd2=0.0;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  VFrow  *row;

  nsurf_g   = (double)enclosure->npatches_g;
  nsurf_l   = enclosure->npatches_l;
  row       = enclosure->row;   
  for (i=0; i<nsurf_l; i++, row++) {
    err0  = row->raw_rowsum-1.0;
    min0  = MIN(min0,err0);
    max0  = MAX(max0,err0);
    sum0 += err0;
    sd0  += err0*err0;
    err1  = row->sym_rowsum-1.0;
    min1  = MIN(min1,err1);
    max1  = MAX(max1,err1);
    sum1 += err1;
    sd1  += err1*err1;
    err2  = row->smooth_rowsum-1.0;
    min2  = MIN(min2,err2);
    max2  = MAX(max2,err2);
    sum2 += err2;
    sd2  += err2*err2;
  }
  VF_GlobalMinDouble(&min0);
  VF_GlobalMaxDouble(&max0);
  VF_GlobalSumDouble(&sum0);
  VF_GlobalSumDouble(&sd0);
  VF_GlobalMinDouble(&min1);
  VF_GlobalMaxDouble(&max1);
  VF_GlobalSumDouble(&sum1);
  VF_GlobalSumDouble(&sd1);
  VF_GlobalMinDouble(&min2);
  VF_GlobalMaxDouble(&max2);
  VF_GlobalSumDouble(&sum2);
  VF_GlobalSumDouble(&sd2);
  sd0   = sqrt((nsurf_g*sd0-sum0*sum0)/(nsurf_g*(nsurf_g-1.0)));
  sd1   = sqrt((nsurf_g*sd1-sum1*sum1)/(nsurf_g*(nsurf_g-1.0)));
  sd1   = sqrt((nsurf_g*sd2-sum2*sum2)/(nsurf_g*(nsurf_g-1.0)));
  sum0 /= nsurf_g;
  sum1 /= nsurf_g;
  sum2 /= nsurf_g;
  *mn0  = min0;
  *mx0  = max0;
  *m0   = sum0;
  *s0   = sd0;
  *mn1  = min1;
  *mx1  = max1;
  *m1   = sum1;
  *s1   = sd1;
  *mn2  = min2;
  *mx2  = max2;
  *m2   = sum2;
  *s2   = sd2;
}

void VF_FillUpperDiagonal(void)
{
  int         i, j, master;
  int         nrows_g, nrows_l, icnt;
  float       *rowVF, *colVF, *VFptr;
  double      *areas, rowsum;
  VFrow       *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  
  VF_GetSPbuffer0_ptr(&colVF);
  VF_GetSPbuffer1_ptr(&rowVF);
  VF_GetDPbuffer0_ptr(&areas);
  VF_GetMatrixAreas(areas);
  row     = enclosure->row;
  nrows_g = enclosure->npatches_g-enclosure->partial;
  nrows_l = enclosure->npatches_l;
  for (icnt=0, i=0; i<nrows_g; i++) {
    master = -1;
    if (nrows_l>0 && icnt<nrows_l) {
      if (row->global_index==i) master = VFLIB_Rank;
    }
    VF_GlobalMaxInt(&master);
    VF_GetMatrixCol(i, master, colVF);
    if (master==VFLIB_Rank) {
      VF_GetMatrixRow(icnt, &j, rowVF);
      for (j=i+1; j<nrows_g; j++) {
        rowVF[j] = (float)(areas[j]*colVF[j]/areas[i]);
      }
      if (enclosure->partial) {
        for (rowsum=0.0, VFptr=rowVF, j=0; j<nrows_g; j++) {
          rowsum += *VFptr++;
        }
        rowVF[nrows_g] = (float)(MAX(0.0,1.0-rowsum));
      }
      VF_UpdateMatrixRow(icnt,rowVF);
      VF_SetRawRowsum(icnt);
      icnt++;
      row++;
    }
  }
}

void VF_ReciprocityStats(double *mean, double *sdev)
{
  int    i, j, icnt;
  uint64_t ncnt;
  int    master, nrows_g, nrows_l;
  float  *rowVF, *colVF;
  double err, std, cnt, tmp, *areas;
  VFrow  *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  
  VF_GetSPbuffer0_ptr(&colVF);
  VF_GetSPbuffer1_ptr(&rowVF);
  VF_GetDPbuffer0_ptr(&areas);
  nrows_g = enclosure->npatches_g;
  nrows_l = enclosure->npatches_l;
  ncnt    = 0;
  std     = 0.0;
  err     = 0.0;
  row     = enclosure->row;
  VF_GetMatrixAreas(areas);
  for (icnt=0, i=0; i<nrows_g; i++) {
    master = -1;
    if (nrows_l>0 && icnt<nrows_l) {
      if (row->global_index==i) master = VFLIB_Rank;
    }
    VF_GlobalMaxInt(&master);
    VF_GetMatrixCol(i, master, colVF);
    if (master==VFLIB_Rank) {
      VF_GetMatrixRow(icnt, &j, rowVF);
      for (j=i+1; j<nrows_g; j++) {
        ncnt++;
        tmp  = fabs(areas[j]*colVF[j] - areas[i]*rowVF[j]);
        err += tmp;
        std += tmp*tmp;
      }
      icnt++;
      row++;
    }
  }
  VF_GlobalSumUInt64(&ncnt);
  VF_GlobalSumDouble(&err);
  VF_GlobalSumDouble(&std);
  cnt   = (double)ncnt;
  std   = sqrt((cnt*std-err*err)/(cnt*(cnt-1.0)));
  err  /= cnt;
  *mean = err;
  *sdev = std;
}

void VF_GetSPbuffer0_ptr(float **ptr)
{
  *ptr = SPbuffer0;
}

void VF_GetSPbuffer1_ptr(float **ptr)
{
  *ptr = SPbuffer1;
}

void VF_GetSPbuffer2_ptr(float **ptr)
{
  *ptr = SPbuffer2;
}

void VF_GetDPbuffer0_ptr(double **ptr)
{
  *ptr = DPbuffer0;
}

void VF_GetDPbuffer1_ptr(double **ptr)
{
  *ptr = DPbuffer1;
}

void VF_GetDPbuffer2_ptr(double **ptr)
{
  *ptr = DPbuffer2;
}

void VF_GetDPbuffer3_ptr(double **ptr)
{
  *ptr = DPbuffer3;
}

void VF_GetINTbuffer_ptr(int **ptr)
{
  *ptr = vf_index;
}

int 
vf_get_nsurfaces_max(void)
{
    return(max_patches);
}

void VF_SortMatrixRows(void)
{
  int l,j,ir,i;
  int nelements;
  VFrow       rra;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  nelements = enclosure->npatches_l;
  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        rra = enclosure->row[--l-1];
      } else {
        rra                  = enclosure->row[ir-1];
        enclosure->row[ir-1] = enclosure->row[0];
        if (--ir == 1) {
          enclosure->row[0] = rra;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (enclosure->row[j-1].host_gid<enclosure->row[j].host_gid)) ++j;
        if (rra.host_gid<enclosure->row[j-1].host_gid) {
          enclosure->row[i-1] = enclosure->row[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      enclosure->row[i-1] = rra;
    }
  }
}

VFrow* VF_LocateLocalRow(int host_gid)
{
  int ascnd,bot,mid,top,nelements;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  nelements = enclosure->npatches_l;
  if (nelements>0) {
    bot   = 0;
    top   = nelements-1;
    ascnd = enclosure->row[top].host_gid > enclosure->row[bot].host_gid;
    if (host_gid == enclosure->row[bot].host_gid) return(&(enclosure->row[bot]));
    if (host_gid == enclosure->row[top].host_gid) return(&(enclosure->row[top]));
    while (top-bot > 1) {
      mid = (top+bot) >> 1;
      if (host_gid == enclosure->row[mid].host_gid) return(&(enclosure->row[mid]));
      if ((host_gid>enclosure->row[mid].host_gid) == ascnd) {
        bot = mid;
      } else {
        top = mid;
      }
    }
  }
  return(NULL);
}

double Ddot1(VFenclosure *e, double *v1, double *v2)
{
  int    i, ii;
  double rc = 0.0;
  VFrow  *row;

  for (row=e->row, i=0; i<e->npatches_l; i++, row++) {
    ii  = row->global_index;
    rc += v1[ii] * v2[i];
  }
  VF_GlobalSumDouble(&rc);
  return(rc);
}

double Ddot0(int n, double *v1, double *v2)
{
  int    i, fourx, rest;
  double rc = 0.0;

  fourx = n / 4;
  rest  = n % 4;
  for (i=0; i<fourx; i++) {
    rc += *v1 * *v2;
    v1++; v2++;
    rc += *v1 * *v2;
    v1++; v2++;
    rc += *v1 * *v2;
    v1++; v2++;
    rc += *v1 * *v2;
    v1++; v2++;
  }
  for (i=0; i<rest; i++) {
    rc += *v1 * *v2;
    v1++; v2++;
  }
  return(rc);
}

double Ddot(int n, double *v1, double *v2)
{
  int    i, fourx, rest;
  double rc = 0.0;

  fourx = n / 4;
  rest  = n % 4;
  for (i=0; i<fourx; i++) {
    rc += *v1 * *v2;
    v1++; v2++;
    rc += *v1 * *v2;
    v1++; v2++;
    rc += *v1 * *v2;
    v1++; v2++;
    rc += *v1 * *v2;
    v1++; v2++;
  }
  for (i=0; i<rest; i++) {
    rc += *v1 * *v2;
    v1++; v2++;
  }
  VF_GlobalSumDouble(&rc);
  return(rc);
}

void Dscal(int n, double a, double r[])
{
  int    i;
  double *ptr;

  for (ptr=r, i=0; i<n; i++) {
    *ptr++ *= a;
  }
}

void Daxpy1(VFenclosure *e, double a, double *x, double *y)
{
  int   i, ii;
  VFrow *row;

  for (row=e->row, i=0; i<e->npatches_l; i++, row++) {
    ii     = row->global_index;
    y[ii] += a*x[ii];
  }
}

void Daxpy2(VFenclosure *e, double a, double *x, double *y)
{
  int   i, ii;
  VFrow *row;

  for (row=e->row, i=0; i<e->npatches_l; i++, row++) {
    ii     = row->global_index;
    y[ii] += a*x[i];
  }
}

void Daxpy(int n, double a, double *x, double *y)
{
  int i, rest, fourx;

  fourx = n / 4;
  rest  = n % 4;
  for (i=0; i<fourx; i++) {
    *y = (a * *x) + *y;
    x++;
    y++;
    *y = (a * *x) + *y;
    x++;
    y++;
    *y = (a * *x) + *y;
    x++;
    y++;
    *y = (a * *x) + *y;
    x++;
    y++;
  }
  for (i=0; i<rest; i++) {
    *y = (a * *x) + *y;
    x++;
    y++;
  }
}

void Dcopy(int n, double *x, double *y)
{
  int i, rest, fourx;

  fourx = n / 4;
  rest  = n % 4;
  for (i=0; i<fourx; i++) {
    *y = *x;
    x++;
    y++;
    *y = *x;
    x++;
    y++;
    *y = *x;
    x++;
    y++;
    *y = *x;
    x++;
    y++;
  }
  for (i=0; i<rest; i++) {
    *y = *x;
    x++;
    y++;
  }
}

