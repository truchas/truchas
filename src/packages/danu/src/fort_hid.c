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
