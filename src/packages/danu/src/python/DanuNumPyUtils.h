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
