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

