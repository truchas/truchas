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
 * Danu File Objects
 *
 * Python Interface File Objects
 *
 */

#ifndef DANU_PY_FILE_H
#define DANU_PY_FILE_H

#include <hdf5.h>

#include "DanuFileObject.h"

/* MACROS */
#define GET_FOBJECT_HID(a)           ( (a)->h5->hid )
#define GET_FOBJECT_HOBJECT(a)       ( (a)->h5 )
#define GET_FOBJECT_ACCESS(a)        ( (a)->access )

/* Prototypes */
File * allocate_file_object(hid_t, const char *, int);
void   deallocate_file_object(File *f);

File * construct_file_object(const char *, const char);
void   deconstruct_file_object(File *fh);

#endif


