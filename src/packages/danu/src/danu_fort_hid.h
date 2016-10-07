/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_fort_hid.h
*
*  DANU FORTRAN Interface for HID types
*
*/

#ifndef DANU_FORT_HID_H
#define DANU_FORT_HID_H

#include <hdf5.h>

/* Simple struct to hold the hid_t identifier */
typedef struct f_hid_t_struct {
  hid_t id;
} *hid_t_ptr;

#define GET_HID_VALUE(a)        ( (*a)->id )
#define SET_HID_VALUE(a,b)      ( (*a)->id = b )


hid_t_ptr create_hid_struct(hid_t id);
void      destroy_hid_struct(hid_t_ptr p);


#endif /* DANU_FORT_HID_H */

