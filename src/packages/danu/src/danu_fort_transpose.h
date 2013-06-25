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
* danu_fort_transpose.h
*
*  DANU FORTRAN/C transpose size arrays
*
*/

#ifndef DANU_FORT_TRANS_H
#define DANU_FORT_TRANS_H

#include <hdf5.h>

void transpose_array_int(int len, const int *src, int *dest);
void transpose_array_hsize(int len, const hsize_t *src, int *dest);

#endif /* DANU_FORT_TRANS_H */

