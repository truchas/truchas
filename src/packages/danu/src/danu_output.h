/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* 
*
*
*/

#ifndef TOUT_OUTPUT_H
#define TOUT_OUTPUT_H

#include <hdf5.h>

/* Public defines */
#define DANU_VERSION_MAJOR 1
#define DANU_VERSION_MINOR 0

/* Public Functions */
hbool_t output_file_is_valid(hid_t fid);

hid_t  output_file_open(const char *filename, unsigned access, unsigned action);
herr_t output_file_create(const char *filename, hid_t *fid);
herr_t output_file_open_rdonly(const char *filename, hid_t *fid);
herr_t output_file_open_append(const char *filename, hid_t *fid);
herr_t output_file_open_rdwr(const char * filename, hid_t *fid);
herr_t output_file_close(hid_t *fid);

#endif
