/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C Transpose arrays
*
*/
#include <hdf5.h>

#include <danu_fort_transpose.h>

void transpose_array_int(int len, const int *src, int * dest)
{
  int i,j;
  j=len-1;
  for(i=0;i<len;i++) {
    dest[i]=src[j];
    j--;
  }
}

void transpose_array_hsize(int len, const hsize_t *src, int * dest)
{
  int i,j;
  j=len-1;
  for(i=0;i<len;i++) {
    dest[i]=(int)src[j];
    j--;
  }
}

