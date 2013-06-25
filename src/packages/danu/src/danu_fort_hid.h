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

