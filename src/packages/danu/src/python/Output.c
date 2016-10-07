/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu Output Objects
 *
 * Python Interface for Output Objects
 *
 */

/* System includes */

#include <hdf5.h>
#include <Python.h>

/* Danu includes */
#include <danu.h>

#include "DanuH5Object.h" 
#include "DanuFile.h"

#include "DanuOutput.h"

Output * allocate_output_object(hid_t hid, const char *name, int access)
{
    Output *f = DANU_MALLOC(Output,1);
    if ( f != NULL ) {
        f->h5     = allocate_h5_object(hid,name);
	f->access = access;
    }
    return f;
}

void deallocate_output_object(Output *f)
{
    deallocate_h5_object(f->h5);
    DANU_FREE(f);
}

Output * construct_output_object(const char *outputname, const char mode)
{
  hid_t fid;
  herr_t err;
  Output *fh = NULL;
  int access = 0;

  if ( FILE_EXISTS(outputname) ) {
    switch(mode) {
      case 'r': case 'R':
        err = output_file_open_rdonly(outputname,&fid);
        access = DANU_FILE_ACC_RDONLY;
        break;
      case 'w': case 'W':
        err = output_file_open_rdwr(outputname,&fid);
        access = DANU_FILE_ACC_RDWR;
        break;
      case 'a': case 'A':
        err = output_file_open_append(outputname,&fid);
        access = DANU_FILE_ACC_APPEND;
        break;
      default:
        err = output_file_open_rdonly(outputname,&fid);
        access = DANU_FILE_ACC_RDONLY;
        break;
    }
  }
  else {
    err = output_file_create(outputname,&fid);
    access = DANU_FILE_ACC_RDWR;
  }
  if ( H5_ISA_VALID_ID(fid) ) {
    fh = allocate_output_object(fid,outputname,access);
  }
  return fh;
}

void deconstruct_output_object(Output *f)
{
    hid_t fid = GET_H5OBJECT_ID(f->h5);
    //printf("Destroying output object\n\n");
    if ( H5_ISA_VALID_ID(fid) ) {
        output_file_close(&fid);
    }
    deallocate_output_object(f);
}

PyObject * meshPythonList(Output *f)
{
  hid_t fid = GET_OUTPUT_OBJ_HID(f);
  herr_t err;

  PyObject * list = NULL;
  PyObject * item;
  char **names;
  int num,i,ret, num_found;
  Py_ssize_t py_num, py_idx;

  err = mesh_count(fid,&num);
  printf("\n\n\n%s Found %d meshes \n",__FILE__,num);
  if ( H5_RETURN_OK(err) && num > 0  ) {

    if ( sizeof(Py_ssize_t) < sizeof(int) ) {
      DANU_WARN_MESS("Warning size of Py_ssize_t < size of int");
    }
    py_num = (Py_ssize_t) num;

    names = DANU_MALLOC(char *, num);
    if ( H5_RETURN_OK(mesh_list(fid, num,names,&num_found) ) ) {
      printf("mesh_list found %d meshes\n", num_found);
      list = PyList_New(py_num);
      for(i=0; i < num; i++) {
        item = PyString_FromString(names[i]);
	  py_idx = (Py_ssize_t) i;
	  ret = PyList_SetItem(list,py_idx,item);
      }
    }
    else {
      DANU_ERROR_MESS("Failed to find mesh names");
    }
  }
  else {
    if ( H5_RETURN_FAIL(err) ) {
      DANU_ERROR_MESS("Failed to count meshes in output file");
    }

    if ( !num ) {
      DANU_ERROR_MESS("No meshes found in output file");
    }

  }

  return list;
}




        
        



