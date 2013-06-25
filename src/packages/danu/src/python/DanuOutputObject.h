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
 * Danu Output Objects
 *
 * Python Interface Ouput Objects
 *
 */

#ifndef DANU_PY_OUTPUT_OBJ_H
#define DANU_PY_OUTPUT_OBJ_H

#include "DanuH5Object.h"

/* HARD CODED from HDF5 public include file H5Fpublic.h 
 * FILE_ACC_RDONLY H5F_ACC_RDONLY
 * FILE_ACC_RDWR   H5F_ACC_RDWR
*/
#define FILE_ACC_RDONLY 0x0000u
#define FILE_ACC_RDWR   0x0001u
#define FILE_ACC_APPEND FILE_ACC_RDWR

typedef struct {

  H5Obj *h5;
  int  access;

} Output;  



#endif

