/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

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

