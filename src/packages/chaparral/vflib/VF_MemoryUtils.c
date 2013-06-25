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
@(#)    $RCSfile: VF_MemoryUtils.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_MemoryUtils.c,v $
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

#include "vf.h"

void VF_MemoryError(int size, char* type, char* function, char* file, int lineno);

void *VF_Malloc_void(int size, char *file, int lineno)
{
  void* m;

  if (size==0) return(NULL);
  m = (void *)malloc(size);
  if (m==NULL) {
    VF_MemoryError(size, "bytes", "allocate", file, lineno);
  } else {
    memset(m,0,size);
  }
  return(m);
}

char *VF_Malloc_char(int size, char *file, int lineno)
{
  char *m;

  if (size==0) return(NULL);
  m = (char *)malloc(size*sizeof(char));
  if (m==NULL) {
    VF_MemoryError(size, "chars", "allocate", file, lineno);
  } else {
    memset(m,0,size);
  }
  return(m);
}

int *VF_Malloc_int(int size, char *file, int lineno)
{
  int  *m;

  if (size==0) return(NULL);
  m = (int *)malloc(size*sizeof(int));
  if (m==NULL) {
    VF_MemoryError(size, "ints", "allocate", file, lineno);
  } else {
    memset(m,0,size*sizeof(int));
  }
  return(m);
}

long *VF_Malloc_long(int size, char *file, int lineno)
{
  long  *m;

  if (size==0) return(NULL);
  m = (long *)malloc(size*sizeof(long));
  if (m==NULL) {
    VF_MemoryError(size, "longs", "allocate", file, lineno);
  } else {
    memset(m,0,size*sizeof(long));
  }
  return(m);
}

float *VF_Malloc_float(int size, char *file, int lineno)
{
  float *m;

  if (size==0) return(NULL);
  m = (float *)malloc(size*sizeof(float));
  if (m==NULL) {
    VF_MemoryError(size, "floats", "allocate", file, lineno);
  } else {
    memset(m,0,size*sizeof(float));
  }
  return(m);
}

double *VF_Malloc_double(int size, char *file, int lineno)
{
  double *m;

  if (size==0) return(NULL);
  m = (double *)malloc(size*sizeof(double));
  if (m==NULL) {
    VF_MemoryError(size, "doubles", "allocate", file, lineno);
  } else {
    memset(m,0,size*sizeof(double));
  }
  return(m);
}

void **VF_MallocPtr_void(int size, char *file, int lineno)
{
  void **m;

  if (size==0) return(NULL);
  m = (void **)malloc(size*sizeof(void*));
  if (m==NULL) VF_MemoryError(size, "void ptrs", "allocate", file, lineno);
  return(m);
}

char **VF_MallocPtr_char(int size, char *file, int lineno)
{
  char **m;

  if (size==0) return(NULL);
  m = (char **)malloc(size*sizeof(char*));
  if (m==NULL) VF_MemoryError(size, "char ptrs", "allocate", file, lineno);
  return(m);
}

int **VF_MallocPtr_int(int size, char *file, int lineno)
{
  int  **m;

  if (size==0) return(NULL);
  m = (int **)malloc(size*sizeof(int*));
  if (m==NULL) VF_MemoryError(size, "int ptrs", "allocate", file, lineno);
  return(m);
}

long **VF_MallocPtr_long(int size, char *file, int lineno)
{
  long **m;

  if (size==0) return(NULL);
  m = (long **)malloc(size*sizeof(long*));
  if (m==NULL) VF_MemoryError(size, "long ptrs", "allocate", file, lineno);
  return(m);
}

float **VF_MallocPtr_float(int size, char *file, int lineno)
{
  float **m;

  if (size==0) return(NULL);
  m = (float **)malloc(size*sizeof(float*));
  if (m==NULL) VF_MemoryError(size, "float ptrs", "allocate", file, lineno);
  return(m);
}

double **VF_MallocPtr_double(int size, char *file, int lineno)
{
  double **m;

  if (size==0) return(NULL);
  m = (double **)malloc(size*sizeof(double*));
  if (m==NULL) VF_MemoryError(size, "double ptrs", "allocate", file, lineno);
  return(m);
}

void *VF_Realloc_void(void *p, int size, char *file, int lineno)
{
  void* m;

  m = (void *)realloc(p,size);
  if (m==NULL) VF_MemoryError(size, "bytes", "reallocate", file, lineno);
  return(m);
}

char *VF_Realloc_char(char *p, int size, char *file, int lineno)
{
  char *m;

  m = (char *)realloc(p,size*sizeof(char));
  if (m==NULL) VF_MemoryError(size, "chars", "reallocate", file, lineno);
  return(m);
}

int *VF_Realloc_int(int *p, int size, char *file, int lineno)
{
  int  *m;

  m = (int *)realloc(p,size*sizeof(int));
  if (m==NULL) VF_MemoryError(size, "ints", "reallocate", file, lineno);
  return(m);
}

long *VF_Reallc_long(long *p, int size, char *file, int lineno)
{
  long *m;

  m = (long *)realloc(p,size*sizeof(long));
  if (m==NULL) VF_MemoryError(size, "longs", "reallocate", file, lineno);
  return(m);
}

float *VF_Realloc_float(float *p, int size, char *file, int lineno)
{
  float *m;

  m = (float *)realloc(p,size*sizeof(float));
  if (m==NULL) VF_MemoryError(size, "floats", "reallocate", file, lineno);
  return(m);
}

double *VF_Realloc_double(double *p, int size, char *file, int lineno)
{
  double *m;

  m = (double *)realloc(p,size*sizeof(double));
  if (m==NULL) VF_MemoryError(size, "doubles", "reallocate", file, lineno);
  return(m);
}

void VF_MemoryError(int size, char* type, char* func, char *file, int lineno) 
{
  printf("Error: Cannot %s memory (%d %s) in %s, line %d\n",func,size,type,file,lineno);
  printf("       Current memory usage is at %.2f Mb\n",VF_GetMemoryUsage());
  VF_Exit(1);
}

void VF_MemoryInfo(char *s) 
{
#ifdef __PUMAGON__
  int fragments;    /* total number of links */
  int total_free;   /* total free memory (bytes) */
  int largest_free; /* largest free link (bytes) */
  int total_used;   /* total currently malloced memory (bytes) */

  heap_info(&fragments,&total_free,&largest_free,&total_used);
  VF_PrintSyncStart();
  if (VFLIB_Rank==0) printf("%s\n",s);
  printf("P%4.4d> Memory layout info:\n",VFLIB_Rank);
  printf("       Total memory usage = %.2f Mb\n",(double)total_used/1048576.0);
  printf("       Total free memory  = %.2f Mb\n",(double)total_free/1048576.0);
  printf("       Largest free link  = %.2f Mb\n",(double)largest_free/1048576.0);
  printf("       Number of links    = %d\n",fragments);
  VF_PrintSyncEnd();
#else
  VF_PrintSyncStart();
  if (VFLIB_Rank==0) printf("%s\n",s);
  printf("P%4.4d> Memory layout info:\n",VFLIB_Rank);
  printf("       Total memory usage = %.2f Mb\n",VF_GetMemoryUsage());
  VF_PrintSyncEnd();
#endif
}

double
VF_GetMemoryUsage(void)
{
  int size=0;
#ifdef __PUMAGON__
  int fragments;    /* total number of links     */
  int total_free;   /* total free memory (bytes) */
  int largest_free; /* largest free link (bytes) */

  heap_info(&fragments,&total_free,&largest_free,&size);
#elif !defined (WIN32)
  char s[20];
  void *addr;
    
  addr = (void*)sbrk(0);
  sprintf(s,"%u",addr);
  sscanf(s,"%d",&size);
#endif
  /* CONVERT TO MEGEBYTES */
  return (double)size/1048576.0;
}

