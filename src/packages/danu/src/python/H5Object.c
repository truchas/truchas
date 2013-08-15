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
 * Python Interface for HDF5 Objects
 *
 */

/* System includes */
#include<string.h>

/* External Package Includes */
#include <hdf5.h>
#include <Python.h>

/* Danu includes */
#include <danu.h>

#include "DanuH5Object.h" 


H5Obj * allocate_h5_object(hid_t hid, const char *name)
{
    H5Obj *h5 = DANU_MALLOC(H5Obj,1);
    if ( h5 != NULL ) {
	h5->name = DANU_MALLOC(char,strlen(name)+1);
	strcpy(h5->name, name);
	h5->hid = hid;
    }

    return h5;
}

void deallocate_h5_object(H5Obj *h5)
{
    DANU_FREE(h5->name);
    DANU_FREE(h5);
}

herr_t writePythonAttribute(H5Obj *h5, const char *name, PyObject *value)
{
    hid_t hid = GET_H5OBJECT_ID(h5);
    herr_t err;

    if ( PyInt_Check(value) ) {
	err = danu_attr_write_int(hid,name,(int)PyInt_AsLong(value));
    }
    else if ( PyFloat_Check(value) ) {
	err = danu_attr_write_double(hid,name,PyFloat_AsDouble(value));
    }
    else if ( PyString_Check(value) ) {
	err = danu_attr_write_string(hid,name,PyString_AsString(value));
    }
    else {
	DANU_ERROR_MESS("Unsupported Python data type");
	err = DANU_FAILURE;
    }

    return err;
}

PyObject * readPythonAttribute(H5Obj *h5, const char *name)
{
  hid_t hid = GET_H5OBJECT_ID(h5);
  int exists;
  herr_t err = danu_attr_exists(hid,name,&exists);
  char type;
  int ivalue;
  double dvalue;
  float fvalue;
  char * string;
  size_t str_size;

  if ( H5_RETURN_OK(err) && exists ) {

      type = danu_attr_get_type(hid,name);

      switch(type) {
	  case 'i':
	      err = danu_attr_read_int(hid,name,&ivalue);
	      return PyInt_FromLong((long)ivalue);
	      break;
	  case 'd': 
	      err = danu_attr_read_double(hid,name,&dvalue);
	      return PyFloat_FromDouble(dvalue);
	      break;
	  case 'f':
	      err = danu_attr_read_float(hid,name,&fvalue);
	      return PyFloat_FromDouble((double)fvalue);
	      break;
	  case 'c':
	      str_size = danu_attr_get_size(hid,name);
	      if ( str_size ) {
		  string = DANU_MALLOC(char,str_size+1);
		  err = danu_attr_read_string(hid,name,string,str_size+1);
		  if ( H5_RETURN_OK(err) ) {
		      return PyString_FromString(string);
		  }
		  else {
		      return NULL;
		  }
	      }
	      else {
		  return NULL;
	      }
	      break;
	  case '\0':
	      return NULL;
	      break;
	  default:
	      DANU_ERROR_MESS("Unsupported Python data type");
	      return NULL;
	      break;
      }

  }
  else {

      if ( !exists ) {
	  DANU_ERROR_MESS("Attribute does not exist");
      }

      if ( H5_RETURN_FAIL(err) ) {
	  DANU_ERROR_MESS("Failed to status the existence of attribute");
      }
  }

  return NULL;
}

PyObject * readPythonAllAttributes(H5Obj *h5)
{
    PyObject * dict = NULL;
    hid_t hid = GET_H5OBJECT_ID(h5);
    char **attr_names;
    size_t *size;
    int num, i, num_found, ret;
    herr_t err = danu_attr_count(hid,&num);
    PyObject *value;

    const int max_name_len = 256;

    if ( H5_RETURN_OK(err) && (num > 0) ) {

	attr_names = DANU_MALLOC(char *,num);
	size = DANU_MALLOC(size_t, num);
	for(i=0; i<num; i++) {
	    size[i] = max_name_len;
	    attr_names[i] = DANU_MALLOC(char,max_name_len);
	}

	err = danu_attr_names(hid,num,attr_names,&num_found);
        if ( H5_RETURN_OK(err) && ( num_found > 0 ) ) {
            dict = PyDict_New();
	    for(i=0; i<num; i++ ) {
		value = readPythonAttribute(h5,attr_names[i]);
                ret = PyDict_SetItemString(dict,attr_names[i],value);
	    }
	}
	else {
	    DANU_ERROR_MESS("Failed to retrieve attribute names");
	}

	DANU_FREE(size);
    }
    else {
	if ( H5_RETURN_FAIL(err) ) {
	    DANU_ERROR_MESS("Failed to count number of attributes");
	}
#if 0
	if ( num == 0 ) {
	    DANU_ERROR_MESS("No attributes present");
	}
#endif
    }

    return dict;
}

















      




