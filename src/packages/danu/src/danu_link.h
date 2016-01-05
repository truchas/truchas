/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * group.c
 *
 *  DANU Groups 
 *
 *
 *  Purpose:
 *
 *          The HDF5 library
 *
 *
 */

#ifndef DANU_LINK_H
#define DANU_LINK_H

#include <hdf5.h>

#include <danu_types.h>

/* Link types */
typedef enum d_link_t {
    DANU_LINK_SOFT = 1,
    DANU_LINK_HARD = 2,
    DANU_LINK_EXTERNAL =3
} d_link_t;


herr_t danu_link_create_soft(hid_t link_loc, const char *link_name, const char * target);
herr_t danu_link_create_hard(hid_t link_loc, const char *link_name, hid_t obj_loc, const char * obj_name);
herr_t danu_link_create_external(hid_t link_loc, const char *link_name, const char *file, const char *obj_name);

herr_t danu_link_exists(hid_t loc_id, const char * name, hbool_t *flag);

#endif
