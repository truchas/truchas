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
 * Danu HDF5 Objects
 *
 * Python Interface Ouput Objects
 *
 */

#ifndef DANU_PY_H5OBJ_H
#define DANU_PY_H5OBJ_H

#include <hdf5.h>

/* Useful macros */
#define GET_H5OBJECT_ID(ptr)    ( (ptr)->hid )

typedef struct {
  hid_t hid;           /* HDF5 indentifier */
  char  *name;         /* Name of the object */
} H5Obj;

#endif

