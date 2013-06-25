/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */

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

