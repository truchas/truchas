/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu Output Objects
 *
 * Python Interface Ouput Objects
 *
 */

#ifndef DANU_PY_OUTPUT_H
#define DANU_PY_OUTPUT_H

#include <hdf5.h>

#include "DanuOutputObject.h"

/* Object Macros */
#define GET_OUTPUT_OBJ_HID(a)           ( (a)->h5->hid )
#define GET_OUTPUT_OBJ_HOBJECT(a)       ( (a)->h5 )
#define GET_OUTPUT_OBJ_ACCESS(a)        ( (a)->access )

/* Prototypes */
Output * allocate_output_object(hid_t,  const char *, int access);
void     deallocate_output_object(Output *f);

Output * construct_output_object(const char *, const char);
void     deconstruct_output_object(Output *fh);


PyObject * meshPythonList(Output *file);





#endif

