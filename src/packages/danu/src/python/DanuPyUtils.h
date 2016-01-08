/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

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
