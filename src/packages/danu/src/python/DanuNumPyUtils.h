/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 *  
 * General NumPy Interface Utilities
 *
 */

#ifndef DANU_NUMPY_UTIL_H
#define DANU_NUMPY_UTIL_H

#include <hdf5.h>

#include <Python.h>

int numpy_typecode(PyObject *t);
int numpy_dataset_typecode(hid_t loc, const char *name);
int numpy_convert_typecode(int danu_typecode);

#endif
