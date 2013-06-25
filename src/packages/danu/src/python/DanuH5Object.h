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

