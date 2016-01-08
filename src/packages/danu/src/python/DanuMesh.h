/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu Mesh Objects
 *
 * Python Interface Ouput Objects
 *
 */

#ifndef DANU_PY_MESH_H
#define DANU_PY_MESH_H

#include <hdf5.h>

#include <danu.h>

#include "DanuOutputObject.h"
#include "DanuMeshObject.h"



#define GET_MOBJECT_HID(a)           ( (a)->h5->hid )
#define GET_MOBJECT_H5OBJECT(a)      ( (a)->h5 )
#define GET_MOBJECT_DIM(a)           ( (a)->dim )
#define GET_MOBJECT_TYPE(a)          ( (a)->mesh_type )
#define GET_MOBJECT_ELEM(a)          ( (a)->elem_type )


/* Prototypes */
Mesh * allocate_mesh_object(hid_t mid, const char *meshname,tmesh_t mtype,telem_t etype);
void   deallocate_mesh_object(Mesh *m);

Mesh * construct_mesh_object(Output *fh, const char *meshname, tmesh_t mtype, telem_t etype);
Mesh * construct_mesh_object_by_id(hid_t mid);
void   deconstruct_mesh_object(Mesh * m);

#endif

