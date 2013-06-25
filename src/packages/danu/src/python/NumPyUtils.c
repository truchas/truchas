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

#include <Python.h>
#include <numpy/arrayobject.h>

#include <danu_dataset_types.h>
#include <danu_dataset.h>

#include "DanuPyException.h"
#include "DanuNumPyUtils.h"

/* External functions built from the numpy.i SWIG interface file*/
extern const char * pytype_string(PyObject *t);


int numpy_typecode(PyObject *input) 
{
  const char *type_string = pytype_string(input);
  int code = NPY_NOTYPE;

  if ( strcmp(type_string,"double") == 0 ) {
    code = NPY_DOUBLE;
  }
  else if ( strcmp(type_string, "int") == 0 ) {
    code = NPY_INT;
  }
  else if ( strcmp(type_string, "double") == 0 ) {
    code = NPY_FLOAT;
  }
  
  return code; 

}

int numpy_dataset_typecode(hid_t loc, const char *name)
{
  int danu_code;
  int ret_code;

  if ( danu_dataset_type(loc,name,&danu_code) >= 0 ) {
    ret_code=numpy_convert_typecode(danu_code);
  }
  else {
    throw_exception("Failed to stat dataset type.");
  }

  return ret_code;

}

int numpy_convert_typecode(int danu_typecode)
{
  int ret_code = NPY_NOTYPE;

  if (danu_typecode == DANU_DATASET_INT ) {
    ret_code=NPY_INT;
  }
  else if (danu_typecode == DANU_DATASET_FLOAT ) {
    ret_code=NPY_FLOAT;
  }
  else if (danu_typecode == DANU_DATASET_DOUBLE ) {
    ret_code=NPY_DOUBLE;
  }
  else if (danu_typecode == DANU_DATASET_STRING ) {
    ret_code=NPY_STRING;
  }
  else {
    throw_exception("Unsupported Danu datatype.");
  }

  return ret_code;

}
