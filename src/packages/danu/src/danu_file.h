/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * danu_file.h
 *
 *  DANU Files 
 *
 *
 *  Purpose:
 *
 *          The HDF5 library
 *
 *
 */

#ifndef DANU_FILE_H
#define DANU_FILE_H

#include <hdf5.h>

/* useful defines */
#define DANU_FILE_ACT_OPEN   0x1000 
#define DANU_FILE_ACT_CREATE 0x2000 

#define DANU_FILE_ACC_RDONLY H5F_ACC_RDONLY
#define DANU_FILE_ACC_RDWR   H5F_ACC_RDWR
#define DANU_FILE_ACC_APPEND H5F_ACC_RDWR


#define ACCESS_IS_VALID(a)   (    (a) == DANU_FILE_ACC_RDONLY \
                               || (a) == DANU_FILE_ACC_RDWR \
                               || (a) == DANU_FILE_ACC_APPEND )
#define ACCESS_IS_INVALID(a) ( ! ACCESS_IS_VALID(a) )

#define ACTION_IS_VALID(a)    ( (a) == DANU_FILE_ACT_OPEN || (a) == DANU_FILE_ACT_CREATE )
#define ACTION_IS_INVALID(a)  ( ! ACTION_IS_VALID(a) )

hid_t danu_file_create(const char *name);
hid_t danu_file_open(const char *name, unsigned access, unsigned action);

hid_t danu_file_open_rdonly(const char *name);
hid_t danu_file_open_append(const char *name);
hid_t danu_file_open_rdwr(const char *name);

herr_t danu_file_flush(hid_t fid);
herr_t danu_file_flush_global(hid_t fid);
herr_t danu_file_close(hid_t fid);


#endif

