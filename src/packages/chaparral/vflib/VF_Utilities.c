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
@(#)    $RCSfile: VF_Utilities.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_Utilities.c,v $
@(#)
@(#)    DESCRIPTION:  Miscellanious utilities.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(sun) || defined(sgi) || defined(dec) || defined(aix) || defined(linux)
#  include <unistd.h>
#endif

#ifdef VF_NO_MPI
#  if defined(sun) || defined(sgi) || defined(dec) || defined(aix) || defined(linux)
#    define USE_GETRUSAGE
#    include <sys/resource.h>
#  else
#    if defined(hp) || defined (hpux) || defined (WIN32)
#      define USE_CLOCK
#      include <time.h>
#    else
#      define USE_TIMEOFDAY
#      include <sys/time.h>
#    endif
#  endif
#endif

#include "vf.h"

void VF_Error(char *s)
{
  fprintf(stdout,"Error: %s\n",s);
  VF_Exit(1);
}

double VF_Clock(void)
{
#ifndef VF_NO_MPI
  return (MPI_Wtime());
#else
#ifdef USE_GETRUSAGE
  struct rusage rusage;
  long   secs,mics;

  getrusage(RUSAGE_SELF,&rusage);
  secs = rusage.ru_utime.tv_sec  + rusage.ru_stime.tv_sec;
  mics = rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec;
  return ((double)secs+(double)mics*1.0e-6);
#else
#ifdef USE_CLOCK
  static clock_t last_num_ticks     = 0;
  static int     clock_rollovers    = 0;
  static double  inv_clocks_per_sec = 1.0/(double)CLOCKS_PER_SEC;
  static double  clock_width        = (double)(1L<<((int)sizeof(clock_t)*8-2))
                                        *4.0/(double)CLOCKS_PER_SEC;    
  double  value;
  clock_t num_ticks;
    
  num_ticks      = clock();
  value          = num_ticks*inv_clocks_per_sec;
  last_num_ticks = num_ticks;
  if (num_ticks<last_num_ticks) clock_rollovers++;
  if (clock_rollovers) value += clock_rollovers*clock_width;
  return (value);
#else
  struct timeval tv;
    
  (void) gettimeofday(&tv, (struct timezone *)NULL);
  return ((double)tv.tv_sec + ((double)tv.tv_usec/1000000.0));
#endif
#endif
#endif
}

#ifdef WIN32
void sleep(clock_t wait)
{
  clock_t goal;
  goal = wait*CLOCKS_PER_SECOND+clock();
  while (goal>clock());
}
#endif

/*+++++++++++++++++++++++++++++++++++++++++++++++++++
  use Heapsort algorithm to sort in ascending order
---------------------------------------------------*/
void VF_SortIntegerArray(int array[], int nelements)
{
  int l,j,ir,i,rra;

  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        rra = array[--l-1];
      } else {
        rra         = array[ir-1];
        array[ir-1] = array[0];
        if (--ir == 1) {
          array[0] = rra;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (array[j-1]<array[j])) ++j;
        if (rra<array[j-1]) {
          array[i-1] = array[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      array[i-1] = rra;
    }
  }
}

void VF_SortIntegerPairArray(int *index, int *data, int nelements)
{
  int   l,j,ir,i,k,tmp_index,tmp_data;

  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        k         = --l-1;
        tmp_index = index[k];
        tmp_data  = data[k];
      } else {
        tmp_index   = index[ir-1];
        index[ir-1] = index[0];
        tmp_data    = data[ir-1];
        data[ir-1]  = data[0];
        if (--ir == 1) {
          index[0] = tmp_index;
          data[0]  = tmp_data;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (index[j-1]<index[j])) ++j;
        if (tmp_index<index[j-1]) {
          index[i-1] = index[j-1];
          data[i-1]  = data[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      index[i-1] = tmp_index;
      data[i-1]  = tmp_data;
    }
  }
}

int VF_LocateIndexInArray(int element, int list[], int nelements)
{
  int i;

  if (nelements>0) {
    for (i=0; i<nelements; i++) {
      if (element == list[i]) return(i);
    }
  }
  return(-1);
}

int VF_LocateIndexInSortedArray(int element, int list[], int nelements)
{
  int ascnd,bot,mid,top;

  if (nelements>0) {
    bot   = 0;
    top   = nelements-1;
    ascnd = list[top] > list[bot];
    if (element == list[bot]) return(bot);
    if (element == list[top]) return(top);
    while (top-bot > 1) {
      mid = (top+bot) >> 1;
      if (element == list[mid]) return(mid);
      if ((element>list[mid]) == ascnd) {
        bot = mid;
      } else {
        top = mid;
      }
    }
  }
  return(-1);
}

void VF_SortSparseArray(VFsparse_array *array)
{
  int   nelements=array->cnt;
  int   *index=array->index;
  int   l,j,ir,i,k,tmp_index;
  float *data=array->data;
  float tmp_data;

  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        k         = --l-1;
        tmp_index = index[k];
        tmp_data  = data[k];
      } else {
        tmp_index   = index[ir-1];
        index[ir-1] = index[0];
        tmp_data    = data[ir-1];
        data[ir-1]  = data[0];
        if (--ir == 1) {
          index[0] = tmp_index;
          data[0]  = tmp_data;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (index[j-1]<index[j])) ++j;
        if (tmp_index<index[j-1]) {
          index[i-1] = index[j-1];
          data[i-1]  = data[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      index[i-1] = tmp_index;
      data[i-1]  = tmp_data;
    }
  }
}

void VF_SortSparseArrayAux(int *index, float *data, int nelements)
{
  int   l,j,ir,i,k,tmp_index;
  float tmp_data;

  if (nelements>1) {
    l  = (nelements>>1)+1;
    ir = nelements;
    for (;;) {
      if (l>1) {
        k         = --l-1;
        tmp_index = index[k];
        tmp_data  = data[k];
      } else {
        tmp_index   = index[ir-1];
        index[ir-1] = index[0];
        tmp_data    = data[ir-1];
        data[ir-1]  = data[0];
        if (--ir == 1) {
          index[0] = tmp_index;
          data[0]  = tmp_data;
          return;
        }
      }
      i = l;
      j = l<<1;
      while (j<=ir) {
        if ((j<ir) && (index[j-1]<index[j])) ++j;
        if (tmp_index<index[j-1]) {
          index[i-1] = index[j-1];
          data[i-1]  = data[j-1];
          j += (i=j);
        } else {
          j = ir+1;
        }
      }
      index[i-1] = tmp_index;
      data[i-1]  = tmp_data;
    }
  }
}


double
VF_SumSparseArray(VFsparse_array *array)
{
  int    i, cnt;
  float  *data;
  double sum=0.0;

  cnt  = array->cnt;
  data = array->data;
  for (i=0; i<cnt; i++) {
    sum += *data++;
  }
  return sum;
}

void
VF_ExpandSparseArray(VFsparse_array *array, float *buffer)
{
  int   i, cnt;
  int   *index;
  float *data;
  
  cnt   = array->cnt;
  index = array->index;
  data  = array->data;
  for (i=0; i<cnt; i++) {
    buffer[index[i]] = data[i];
  }
}

void
VF_InitializeSparseArray(VFsparse_array *array)
{
  array->cnt   = 0;
  array->size  = 0;
  array->data  = NULL;
  array->index = NULL;
}

void
VF_AllocateSparseArray(VFsparse_array *array, int n)
{
  array->cnt   = n;
  array->size  = n;
  array->data  = VF_Newf(n);
  array->index = VF_Newi(n);
}

void
VF_ReallocateSparseArray(VFsparse_array *array, int n)
{
  array->cnt   = n;
  array->size  = n;
  array->data  = VF_ReNewf(array->data, n);
  array->index = VF_ReNewi(array->index, n);
}

void
VF_FreeSparseArray(VFsparse_array *array)
{
  array->cnt  = 0;
  array->size = 0;
  VF_Free(array->data);
  VF_Free(array->index);
}

#define M  714025
#define IA 1366
#define IC 150889

double ran2(long *idum)
{
  static long iy,ir[100];
  static int  iff=0;
  int    j;
  void   nrerror();

  if (*idum<0 || iff==0) {
    iff = 1;
    if ((*idum=(IC-(*idum))%M)<0) *idum = -(*idum);
    for (j=1; j<=97; j++) {
      *idum = (IA*(*idum)+IC) % M;
      ir[j] = (*idum);
    }
    *idum = (IA*(*idum)+IC) % M;
    iy    = (*idum);
  }
  j = 1 + 97.0*iy/M;
  if (j>97 || j<1) {
    fprintf(stderr,"RAN2: This cannot happen.");
    VF_Exit(0);
  }
  iy    = ir[j];
  *idum = (IA*(*idum)+IC) % M;
  ir[j] = (*idum);
  return (double) iy/M;
}

#undef M
#undef IA
#undef IC
