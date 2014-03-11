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
 * General Python Interface Utilities
 *
 */

#ifndef DANU_PY_UTIL_H
#define DANU_PY_UTIL_H

#include <Python.h>

PyObject * convertCharListToPyList(const char * const * names, int num);
PyObject * convertIntListToPyList(int *int_list, int num);
PyObject * buildPyList(PyObject **objects, int num);

#endif
