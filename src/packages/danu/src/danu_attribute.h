/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_attribute.h
*
*  DANU scalar attributes
*
*/

#ifndef DANU_ATTR_H
#define DANU_ATTR_H

#include <hdf5.h>

#include <danu_error.h>

herr_t danu_attr_exists(hid_t loc, const char *attr_name, int *flag);
herr_t danu_attr_count(hid_t loc, int *num_found);
herr_t danu_attr_names(hid_t loc, int num, char **names, int *num_found);

char danu_attr_get_type(hid_t loc, const char *attr_name);
size_t danu_attr_get_size(hid_t loc, const char *name);

herr_t danu_attr_write(hid_t loc,const char * name, void * value, hid_t type);
herr_t danu_attr_read(hid_t loc,const char * name, void * buffer, hid_t type);

herr_t danu_attr_write_int(hid_t loc,const char * name, int value);
herr_t danu_attr_read_int(hid_t loc,const char * name, int * buffer);

herr_t danu_attr_write_uint(hid_t loc,const char * name, unsigned int value);
herr_t danu_attr_read_uint(hid_t loc,const char * name, unsigned int * buffer);

herr_t danu_attr_write_float(hid_t loc,const char * name, float value);
herr_t danu_attr_read_float(hid_t loc,const char * name, float * buffer);

herr_t danu_attr_write_double(hid_t loc,const char * name, double value);
herr_t danu_attr_read_double(hid_t loc,const char * name, double * buffer);

herr_t danu_attr_write_string(hid_t loc,const char * name, const char * string);
herr_t danu_attr_read_string(hid_t loc,const char * name, char * buffer, size_t buf_len);

#endif
