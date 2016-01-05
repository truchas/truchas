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

#ifndef TOUT_OFFSET_H
#define TOUT_OFFSET_H

#define DANU_OFFSET_ATTR_NAME "Offset"

#include <hdf5.h>

herr_t danu_set_offset(hid_t loc, int offset);
herr_t danu_get_offset(hid_t loc, int *offset);

#endif

