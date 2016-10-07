/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu File Objects
 *
 * Python Interface for File Objects
 *
 */

/* System includes */

#include <hdf5.h>
#include <Python.h>

/* Danu */
#include <danu.h>

#include "DanuH5Object.h" 
#include "DanuFile.h"



File * allocate_file_object(hid_t hid, const char *name, int access)
{
    File *f = DANU_MALLOC(File,1);
    if ( f != NULL ) {
        f->h5 = allocate_h5_object(hid,name);
        f->access = access;
    }
    return f;
}

void deallocate_file_object(File *f)
{
    deallocate_h5_object(f->h5);
    DANU_FREE(f);
}

File * construct_file_object(const char *filename, const char mode)
{
  hid_t fid;
  File *fh = NULL;
  int access = 0;
 
  if ( FILE_EXISTS(filename) ) {
 
    switch(mode) {
      case 'r': case 'R':
        fid = danu_file_open_rdonly(filename);
        access = DANU_FILE_ACC_RDONLY;
        break;
      case 'w': case 'W':
        fid = danu_file_open_rdwr(filename);
        access = DANU_FILE_ACC_RDWR;
        break;
      case 'a': case 'A':
        fid = danu_file_open_append(filename);
        access = DANU_FILE_ACC_APPEND;
        break;
      default:
        fid = danu_file_open_rdonly(filename);
        access = DANU_FILE_ACC_RDONLY;
        break;
    }
  
  }
  else {
    fid = danu_file_create(filename);
    access = DANU_FILE_ACC_RDWR;
  }
  if ( H5_ISA_VALID_ID(fid) ) {
    fh = allocate_file_object(fid,filename,access);
  }
  return fh;
}

void deconstruct_file_object(File *fh)
{
    hid_t fid = GET_H5OBJECT_ID(fh->h5);
    if ( H5_ISA_VALID_ID(fid) ) {
        danu_file_close(fid);
    }
    deallocate_file_object(fh);
}


