/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_fort_hid.h
*
*  DANU Fortran HID_T Struct 
*
*/
#include <stdlib.h>
#include<string.h>

#include <hdf5.h>

#include <danu_memory.h>
#include <danu_fort_hid.h>

/* Struct creation, destruction and manipulation routines */
hid_t_ptr create_hid_struct(hid_t id)
{
    hid_t_ptr ptr = DANU_MALLOC(struct f_hid_t_struct,1);

    ptr->id = id;
#if 0
    danu_debug_printf("%s %d ptr=%lx id=%d\n", __FILE__, __LINE__,
	                     ptr, ptr->id);
#endif

    return ptr;
}

void destroy_hid_struct( hid_t_ptr ptr)
{
    DANU_FREE(ptr);
}
