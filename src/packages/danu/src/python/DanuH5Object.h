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

#ifndef DANU_PY_HDF5_OBJECT_H
#define DANU_PY_HDF5_OBJECT_H

#include <hdf5.h>
#include <Python.h>

#include "DanuH5Obj.h"

H5Obj * allocate_h5_object(hid_t hid, const char *name);
void    deallocate_h5_object(H5Obj *h5);


herr_t     writePythonAttribute(H5Obj *h5, const char *name, PyObject *value); 
PyObject * readPythonAttribute(H5Obj *h5, const char *name); 
PyObject * readPythonAllAttributes(H5Obj *h5);



#endif

